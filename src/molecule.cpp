#include "molecule.h"


template <class C, bool xyzless>
void Molecule::constructor(C input_object, py::module_& nx, py::module_& pyutils) {
    checkpoint("Entered SDF constructor");
    this->nx = &nx;
    this->pyutils = &pyutils;
    #ifdef KDMOL_LOG
        if constexpr(!xyzless) {
            Logger::clear_log();
            Logger::log("[CTOR] NEW CALL");
            Logger::log("INIT:: filename = {}", input_object);
        }
    #endif
    
    checkpoint("Setting nx");
    Utils::PyModules::set_nx(this->nx);
    Utils::PyModules::set_pyutils(this->pyutils);
    checkpoint("Done setting nx");
    this->filename = ""; // Some initialization
    
    if constexpr(xyzless) {
        read_graph(py::handle(input_object["graph"]).cast<py::object>());
        mol_data.set_provided_dihedrals(py::handle(input_object["fixed_dihedrals"]).cast<py::dict>());
    } else {
        read_sdf(input_object);
    }
    gen_cyclic_parts();
    gen_dofs();
    gen_fragment_graph();
}

template void Molecule::constructor<py::dict, true>(py::dict input_object, py::module_& nx, py::module_& pyutils);
template void Molecule::constructor<const std::string, false>(const std::string input_object, py::module_& nx, py::module_& pyutils);

void Molecule::read_sdf(const std::string& filename) {
    auto mylines = Utils::readlines(filename);
    
    assertm(mylines.size() > 4, fmt::format("Provided file {} seems to be corrupted", filename));
    auto natoms_str = mylines[3].substr(0, 3),
         nbonds_str = mylines[3].substr(3, 3);
    boost::algorithm::trim(natoms_str);
    boost::algorithm::trim(nbonds_str);
    auto natoms = boost::lexical_cast<unsigned int>(natoms_str), 
         nbonds = boost::lexical_cast<unsigned int>(nbonds_str);
    mol_data = mol_data_t(natoms, nbonds, dofs, discrete_dofs);
    mol_data.molgr = nx->attr("Graph")();
    mol_data.xyzless_init = false;

    std::vector<std::string> parts;
    int idx = 0;
    for (unsigned int i = 4; i < natoms + 4; ++i) {
        assertm(mylines.size() > i, fmt::format("Provided file {} seems to be corrupted", filename));
        boost::algorithm::trim(mylines[i]);
        boost::split(parts, mylines[i], boost::is_any_of(" "), boost::token_compress_on);

        auto xyz = mol_data.xyz[idx];
        xyz[0] = boost::lexical_cast<double>(parts[0]);
        xyz[1] = boost::lexical_cast<double>(parts[1]);
        xyz[2] = boost::lexical_cast<double>(parts[2]);

        mol_data.molgr.attr("add_node")(idx);
        mol_data.atom_symbols[idx] = parts[3];
        parts.empty();
        ++idx;
    }
    
    int edge_count = 0;
    std::vector<std::pair<int, int>> edge_list(nbonds);
    for (unsigned int i = natoms + 4; i < natoms + nbonds + 4; ++i) {
        assertm(mylines.size() > i, fmt::format("Provided file {} seems to be corrupted", filename));
        auto atA_str  = mylines[i].substr(0, 3),
             atB_str  = mylines[i].substr(3, 3),
             bondtype = mylines[i].substr(7, 3);
        boost::algorithm::trim(atA_str);
        boost::algorithm::trim(atB_str);
        boost::algorithm::trim(bondtype);
        auto atA = boost::lexical_cast<unsigned int>(atA_str) - 1,
             atB = boost::lexical_cast<unsigned int>(atB_str) - 1;
        
        mol_data.bondtypes[edge_count] = boost::lexical_cast<unsigned int>(bondtype);
        mol_data.bondlengths[edge_count] = mol_data.xyz.get_distance(atA, atB);

        mol_data.molgr.attr("add_edge")(atA, atB);
        edge_list[edge_count] = std::pair(atA, atB);
        edge_count++;
    }
    mol_data.edge_indexer = Utils::EdgeIndexer(natoms, edge_list);
    mol_data.atom_in_frag = mol_data_t::atom_in_frag_t(natoms);

    for (const auto& node : mol_data.molgr.attr("nodes")())
        mol_data.polys[node.cast<int>()] = poly_t(mol_data.molgr, mol_data.xyz, node.cast<int>());
    
    this->filename = filename;
    checkpoint("readsdf");
    #ifdef KDMOL_LOG
        py::list edge_data, node_data;
        for (const auto& idx_h : mol_data.molgr.attr("nodes")()) {
            py::list new_item;
            idx = idx_h.cast<int>();
            auto xyz = mol_data.xyz[idx];
            new_item.append(idx);
            new_item.append(xyz[0]);
            new_item.append(xyz[1]);
            new_item.append(xyz[2]);
            new_item.append(mol_data.atom_symbols[idx]);
            node_data.append(new_item);
        }
         
        for (const auto& py_edge : mol_data.molgr.attr("edges")()) {
            py::list new_item;
            auto edge = py_edge.cast<std::pair<int, int>>();
            auto edge_idx = mol_data.bond_data(edge).index;
            new_item.append(edge.first);
            new_item.append(edge.second);
            new_item.append(mol_data.bondlengths[edge_idx]);
            new_item.append(mol_data.bondtypes[edge_idx]);
            edge_data.append(new_item);
        }

        log_routine("readsdf");
        log_in("file", this->filename);
        log_check("node_data_SET1", node_data);
        log_check("edge_data_SET2", edge_data);
        log_finish();
    #endif
}

void Molecule::read_graph(py::object input_graph) {
    auto natoms = py::handle(input_graph.attr("number_of_nodes"))().cast<int>(), 
         nbonds = py::handle(input_graph.attr("number_of_edges"))().cast<int>();
    mol_data = mol_data_t(natoms, nbonds, dofs, discrete_dofs);
    mol_data.molgr = nx->attr("Graph")();
    mol_data.xyzless_init = true;
    
    for (const auto& iteridx : input_graph.attr("nodes")()) {
        auto py_idx = py::handle(iteridx).cast<py::int_>();
        mol_data.molgr.attr("add_node")(py_idx);
        mol_data.atom_symbols[py_idx.cast<int>()] = py::handle(input_graph.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("symbol")).cast<std::string>();
    }
    
    int edge_count = 0;
    std::vector<std::pair<int, int>> edge_list(nbonds);
    for (const auto& iteredge : input_graph.attr("edges")()) {
        auto atA = py::handle(iteredge.attr("__getitem__")(0)).cast<int>(),
             atB = py::handle(iteredge.attr("__getitem__")(1)).cast<int>();
        
        mol_data.bondtypes[edge_count] = py::handle(input_graph.attr("__getitem__")(atA).attr("__getitem__")(atB).attr("__getitem__")("bond_type")).cast<int>();
        mol_data.bondlengths[edge_count] = py::handle(input_graph.attr("__getitem__")(atA).attr("__getitem__")(atB).attr("__getitem__")("length")).cast<double>();

        mol_data.molgr.attr("add_edge")(atA, atB);
        edge_list[edge_count] = std::pair(atA, atB);
        edge_count++;
    }
    mol_data.edge_indexer = Utils::EdgeIndexer(natoms, edge_list);
    mol_data.atom_in_frag = mol_data_t::atom_in_frag_t(natoms);

    for (const auto& iteridx : mol_data.molgr.attr("nodes")()) {
        auto py_idx = py::handle(iteridx).cast<py::int_>();
        auto idx = py::handle(iteridx).cast<int>();
        auto passed_poly = py::handle(input_graph.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("poly")).cast<py::object>();
        mol_data.polys[idx] = poly_t(passed_poly, idx);
    }
}

void Molecule::gen_cyclic_parts() {
    py::list py_edges, bridges, nonbridges;
    for (const auto& py_edge : mol_data.molgr.attr("edges")()) {
        auto edge = py_edge.cast<std::pair<int, int>>();
        if (edge.first < edge.second)
            py_edges.append(py::make_tuple(edge.first, edge.second));
        else
            py_edges.append(py::make_tuple(edge.second, edge.first));
    }
    
    auto graph = nx->attr("Graph")();
    graph.attr("add_edges_from")(py_edges);
    for(auto bridge : nx->attr("bridges")(graph)) {
        auto bridge_obj = bridge.cast<py::tuple>();
        auto atA = bridge_obj[0].cast<py::int_>(),
             atB = bridge_obj[1].cast<py::int_>();
        if (atA < atB)
            bridges.append(py::make_tuple(atA, atB));
        else
            bridges.append(py::make_tuple(atB, atA));
    }
    
    for (const auto& py_edge : mol_data.molgr.attr("edges")()) {
        auto edge = py_edge.cast<std::pair<int, int>>();
        py::tuple fixed_edge;
        if (edge.first < edge.second)
            fixed_edge = py::make_tuple(edge.first, edge.second);
        else
            fixed_edge = py::make_tuple(edge.second, edge.first);
        
        if(!bridges.attr("__contains__")(fixed_edge).cast<bool>())
            nonbridges.append(fixed_edge);
    }

    auto cyclgraph = nx->attr("Graph")();
    cyclgraph.attr("add_edges_from")(nonbridges);
    auto conn_comps_raw = py::list(nx->attr("connected_components")(cyclgraph));
    py::list conn_comps; // Needed for reproducible CP order
    for (const auto& item : conn_comps_raw) {
        auto newitem = py::list(py::handle(item).cast<py::set>());
        newitem.attr("sort")();
        conn_comps.append(newitem);

    }
    conn_comps.attr("sort")();
    checkpoint("conncomp_order");
    #ifdef KDMOL_LOG
        int cc_idx = 0;
        for (const auto& item : conn_comps) {
            log_routine("conncomp_order");
            log_in("file", this->filename);
            log_in("cc_idx", cc_idx++);
            log_check("nodes", item.cast<py::list>());
            log_finish();
        }
    #endif

    cyclic_parts.reserve(conn_comps.size());
    for(auto& comp : conn_comps) {
        // auto cp_indices = comp.cast<std::unordered_set<int>>();
        subgraph_t cp_subgraph = mol_data.molgr.attr("subgraph")(comp);
        cyclic_parts.push_back(cyclic_part_t(&mol_data, &dofs, &discrete_dofs, cp_subgraph, this->filename));
    }

    checkpoint("gencyclparts_A");
    #ifdef KDMOL_LOG
        log_routine("gencyclparts_A");
        log_in("file", this->filename);
        log_check("bridges_SET2", bridges);
        log_check("nonbridges_SET2", nonbridges);
        log_finish();
    #endif
}

void Molecule::gen_dofs() {
    #ifdef KDMOL_LOG
        py::list added_bonds;
    #endif

    // Noncyclic DOFs
    for(auto bridge : nx->attr("bridges")(mol_data.molgr)) {
        auto bridge_obj = bridge.cast<py::tuple>();
        auto atA = bridge_obj[0].cast<int>(),
             atB = bridge_obj[1].cast<int>();
        if ((mol_data.number_of_neighbors(atA) > 1) && 
            (mol_data.number_of_neighbors(atB) > 1)) {
            #ifdef KDMOL_DEBUG
                if (atA > atB) { std::swap(atA, atB); }
            #endif
            dofs_t::bondatoms_t bond_atoms {atA, atB};
            dofs_t::sideatoms_t side_atoms {mol_data.get_some_neighbor(atA, atB),
                                            mol_data.get_some_neighbor(atB, atA)};
            
            // enum dof_type : bool { fixed = false, free = true }
            auto doftype = static_cast<dofs_t::dof_type>(mol_data.get_bondtype(bond_atoms) == 1);

            auto dihedral = 0.0;
            if (!mol_data.xyzless_init)
                dihedral = mol_data.xyz.get_dihedral(side_atoms.first, bond_atoms.first,
                                                     bond_atoms.second, side_atoms.second);
            else
                dihedral = mol_data.get_provided_dihedral(side_atoms.first, bond_atoms.first,
                                                          bond_atoms.second, side_atoms.second);
            dofs.new_dof(bond_atoms, side_atoms, doftype, dihedral);

            checkpoint("noncyclicDOF");
            #ifdef KDMOL_LOG
                py::list bridge_log = bridge.cast<py::list>(),
                         sides_log = py::cast(side_atoms).cast<py::list>();

                if (atA > atB) {
                    py::int_ x = bridge_log[0];
                    bridge_log[0] = bridge_log[1];
                    bridge_log[1] = x;
                    
                    x = sides_log[0];
                    sides_log[0] = sides_log[1];
                    sides_log[1] = x;
                }

                std::string type_str;
                if(doftype == dofs_t::fixed)
                    type_str = "fixed";
                else if (doftype == dofs_t::free)
                    type_str = "free";
                else
                    throw std::runtime_error("unknown dof type");

                log_routine("noncyclicDOF");
                log_in("file", this->filename);
                log_in("bond", bridge_log);
                log_check("type", type_str);
                log_check("sides", sides_log);
                log_check("dihedral", as_list(dihedral));
                log_finish();

                added_bonds.append(py::cast(bond_atoms));
            #endif
        }
    }

    checkpoint("noncyclicDOF_total");
    #ifdef KDMOL_LOG
        log_routine("noncyclicDOF_total");
        log_in("file", this->filename);
        log_check("added_bonds_SET2", added_bonds);
        log_finish();
    #endif
    dofs.noncyclic_size = dofs.size();

    // Include cyclic DOFs
    dofs.cpdof_starts.reserve(cyclic_parts.size());
    dofs.cpdof_sizes.reserve(cyclic_parts.size());
    for (const auto& cyclic_part : cyclic_parts) {
        auto cp_dofs = cyclic_part.get_dofs();
        auto cp_dofs_size = cp_dofs.size();
        dofs.cpdof_starts.push_back(dofs.size());

        for (const auto [dof, doftype, dihedral] : cp_dofs) {
            if (!dofs.contain(dof)) {
                dofs_t::bondatoms_t bond_atoms = {std::get<1>(dof), std::get<2>(dof)};
                dofs_t::sideatoms_t side_atoms = {std::get<0>(dof), std::get<3>(dof)};
                dofs.new_dof(bond_atoms, side_atoms, doftype, dihedral);
            } else
                cp_dofs_size--;
        }
        dofs.cpdof_sizes.push_back(cp_dofs_size);
    }
    dofs.finalize();

    // Distribute cyclic DOFs
    for (int i = 0; i < cyclic_parts.size(); ++i)
        cyclic_parts[i].map_dofvalues(dofs.cpdof_starts[i], dofs.cpdof_sizes[i]);

    checkpoint("cp_checkFreeDof");
    #ifdef KDMOL_LOG
        {
        py::list log_freedih;
        for (int i = 0; i < cyclic_parts.size(); ++i) {
            py::list newitem, dih_item;
            newitem.append(i);
            for (int j = dofs.cpdof_starts.at(i); j < dofs.cpdof_starts.at(i) + dofs.cpdof_sizes.at(i); ++j) {
                if (dofs.doftype_container.at(j) != dofs_t::dof_type::free)
                    continue;
                auto dih = std::make_tuple(
                    dofs.side_container.at(j).first,
                    dofs.bond_container.at(j).first,
                    dofs.bond_container.at(j).second,
                    dofs.side_container.at(j).second
                );
                dih_item.append(py::cast(dih));
            }
            dih_item.attr("sort")();
            newitem.append(dih_item);
            log_freedih.append(newitem);
        }
        log_routine("cp_checkFreeDof");
        log_in("file", this->filename);
        log_check("freedih_SET1", log_freedih);
        log_finish();
        }
    #endif

    // Include discrete DOFs of cyclic parts
    discrete_dofs.cpdof_starts.reserve(cyclic_parts.size());
    discrete_dofs.cpdof_sizes.reserve(cyclic_parts.size());
    for (const auto& cyclic_part : cyclic_parts) {
        auto cp_ddofs = cyclic_part.get_discrete_dofs();
        auto cp_ddofs_size = cp_ddofs.size();
        discrete_dofs.cpdof_starts.push_back(discrete_dofs.size());
        discrete_dofs.cpdof_sizes.push_back(cp_ddofs_size);

        for (const auto& owner_atoms : cp_ddofs)
            discrete_dofs.new_dof(owner_atoms);
    }
    discrete_dofs.finalize();

    // Distribute discrete DOFs
    for (int i = 0; i < cyclic_parts.size(); ++i)
        cyclic_parts[i].map_discrete_dofvalues(discrete_dofs.cpdof_starts[i], discrete_dofs.cpdof_sizes[i]);
    
    checkpoint("mol_discrete_dofs");
    #ifdef KDMOL_LOG
        {
        py::list log_ddof;
        for (int cp_idx = 0; cp_idx < cyclic_parts.size(); ++cp_idx) {
            py::list cur_ddofs;
            for (int i = discrete_dofs.cpdof_starts[cp_idx]; 
                 i < discrete_dofs.cpdof_starts[cp_idx] + discrete_dofs.cpdof_sizes[cp_idx]; ++i)
            {
                auto cur_atoms = py::handle(py::cast(discrete_dofs.owner_container[i])).cast<py::list>();
                cur_atoms.attr("sort")();
                cur_ddofs.append(cur_atoms);
            }
            if (cur_ddofs.size() == 0)
                continue;
            cur_ddofs.attr("sort")();

            py::list newitem;
            newitem.append(cp_idx);
            newitem.append(cur_ddofs);
            log_ddof.append(newitem);
        }

        log_routine("mol_discrete_dofs");
        log_in("file", this->filename);
        log_check("discretedofs_SET1", log_ddof);
        log_finish();
        }
    #endif
}

void Molecule::gen_fragment_graph() {
    auto fragments = nx->attr("Graph")();
    fragments.attr("add_edges_from")(mol_data.molgr.attr("edges")());
    for (int i = 0; i < dofs.noncyclic_size; ++i) {
        // if (dofs.doftype_container[i] == dofs_t::free)
        fragments.attr("remove_edge")(dofs.bond_container[i].first, dofs.bond_container[i].second);
    }
    auto frag_comp = py::list(nx->attr("connected_components")(fragments));
    #ifdef KDMOL_LOG
        for (int i = 0; i < frag_comp.size(); ++i) {
            py::list x(frag_comp[i]);
            x.attr("sort")();
            frag_comp[i] = x;
        }
        frag_comp.attr("sort")();
    #endif
    
    mol_data.fg_container.reserve(frag_comp.size());
    std::vector<std::vector<int>> fg_edges(dofs.noncyclic_size);
    std::unordered_map<int, int> edgeidx_param_map;
    #ifdef KDMOL_LOG
        py::list log_fg_nodes;
    #endif
    for (int i = 0; i < frag_comp.size(); ++i) {
        const auto& comp = py::set(frag_comp[i]);
        auto atomidx_set = comp.cast<mol_data_t::fragatoms_t>();

        for (const auto atom_idx : atomidx_set)
            mol_data.atom_in_frag[atom_idx] = i;

        int k = 0;
        for (int j = 0; j < dofs.noncyclic_size; ++j) {
            // if (dofs.doftype_container[j] == dofs_t::fixed)
            //     continue;
            if(atomidx_set.find(dofs.bond_container[j].first) != atomidx_set.cend() ||
               atomidx_set.find(dofs.bond_container[j].second) != atomidx_set.cend()) {
                fg_edges[k].push_back(i);
                edgeidx_param_map[k] = j;
            }
            k++;
        }

        mol_data.fg_container.push_back(std::move(atomidx_set));
        #ifdef KDMOL_LOG
            auto nodes = comp.cast<py::list>();
            nodes.attr("sort")();
            py::list item;
            item.append(py::cast(i));
            item.append(nodes);
            log_fg_nodes.append(item);
        #endif
    }
    checkpoint("gen_frag_nodes");
    #ifdef KDMOL_LOG
        log_routine("gen_frag_nodes");
        log_in("file", this->filename);
        log_check("nodes_SET1", log_fg_nodes);
        log_finish();
    #endif
    mol_data.fg_container.shrink_to_fit();

    // --------------------------------
    // This might be redundant
    std::map<std::pair<int, int>, int> edge_param_map;
    for (int i = 0; i < fg_edges.size(); ++i) {
        if (fg_edges[i][0] < fg_edges[i][1])
            edge_param_map[std::make_pair(fg_edges[i][0], fg_edges[i][1])] = edgeidx_param_map[i];
        else
            edge_param_map[std::make_pair(fg_edges[i][1], fg_edges[i][0])] = edgeidx_param_map[i];
    }
    // --------------------------------
    
    auto fg = nx->attr("Graph")();
    if (frag_comp.size() == 1)
        fg.attr("add_node")(0);
    else
        fg.attr("add_edges_from")(py::cast(fg_edges));
    
    checkpoint("gen_frag_grA");
    #ifdef KDMOL_LOG
        auto nodes = py::list(fg.attr("nodes")());
        nodes.attr("sort")();
        log_routine("gen_frag_grA");
        log_in("file", this->filename);
        log_check("node_data", nodes);
        log_check("edge_data_SET2", py::list(fg.attr("edges")()));
        log_finish();
    #endif

    headnode = 0;
    py::int_ py_headnode = py::cast(headnode);
    std::vector<int> endnodes;
    for (const auto& node : fg.attr("nodes")()) {
        py::int_ node_idx = node.cast<py::int_>();
        if ((py_headnode.cast<int>() != node_idx.cast<int>()) && (py::list(fg.attr("neighbors")(node_idx)).size() == 1))
            endnodes.push_back(node.cast<py::int_>());
    }
    const int number_of_fragments = fg.attr("number_of_nodes")().cast<int>();
    dfg = dfg_t(number_of_fragments);

    for (const auto& enode : endnodes) {
        auto solvepath = py::list(nx->attr("shortest_path")(fg, py::cast(enode), py_headnode));
        for (int i = 0; i < solvepath.size() - 1; ++i) {
            auto vA = vertex(solvepath[i].cast<int>(), dfg),
                 vB = vertex(solvepath[i + 1].cast<int>(), dfg);
            if (!edge(vA, vB, dfg).second) {
                dfg_t::edge_descriptor edge_descr = add_edge(vA, vB, dfg).first;
                if (vA < vB)
                    dfg[edge_descr].param_idx = edge_param_map[{vA, vB}];
                else
                    dfg[edge_descr].param_idx = edge_param_map[{vB, vA}];
            }
        }
    }

    checkpoint("gen_dfg");
    #ifdef KDMOL_LOG
        py::list dfg_nodes, dfg_edges;
        boost::graph_traits <dfg_t>::vertex_iterator i, end;
        boost::property_map <dfg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dfg);

        for (boost::tie(i, end) = vertices(dfg); i != end; ++i) {
            dfg_nodes.append(get(index_map, *i));
        }
        
        boost::graph_traits<dfg_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(dfg); ei != ei_end; ++ei)
            dfg_edges.append(py::make_tuple(get(index_map, source(*ei, dfg)),
                                            get(index_map, target(*ei, dfg))));

        log_routine("gen_dfg");
        log_in("file", this->filename);
        log_check("node_data", dfg_nodes);
        log_check("edge_data_SET2", dfg_edges);
        log_finish();
    #endif

    frag_refs = fragref_container_t(number_of_fragments, nullptr);

    for (auto& cp : cyclic_parts) {
        auto someatom = cp.get_some_atom();
        for (int j = 0; j < mol_data.fg_container.size(); j++) {
            if (mol_data.belongs_to_fragment(someatom, j)) {
                frag_refs[j] = &cp;
                break;
            }
        }
    }

    frag_container = fragment_container_t(number_of_fragments - cyclic_parts.size());
    int cur_idx = 0;
    for (int i = 0; i < number_of_fragments; ++i) {
        if (frag_refs[i] != nullptr)
            continue; // Skip if corresponds to a cyclic part
        frag_container[cur_idx] = fragment_t(&mol_data, &dofs, &discrete_dofs, i);
        frag_refs[i] = &frag_container[cur_idx];
        cur_idx++;
    }

    for (int i = 0; i < frag_refs.size(); ++i)
        frag_refs[i]->configure_assembly(mol_data.molgr, dfg, i);
    
    #ifdef KDMOL_DEBUG
        mol_data.edge_indexer.validate_frag_indices();
    #endif

    checkpoint("frags_done");
    #ifdef KDMOL_LOG
        py::list nodetypes, nodeatoms, edgeparams;
        for (const auto& node : fg.attr("nodes")()) {
            py::int_ node_idx = node.cast<py::int_>();
            py::list new_item, myatoms = py::cast(mol_data.fg_container[node_idx]);
            new_item.append(node_idx);
            myatoms.attr("sort")();
            new_item.append(myatoms);
            nodeatoms.append(new_item);
        }

        int node_idx = 0;
        for (auto fragref : frag_refs) {
            py::list new_item;
            new_item.append(node_idx);
            if (dynamic_cast<cyclic_part_t*>(fragref) != nullptr)
                new_item.append(1);
            else
                new_item.append(0);
            nodetypes.append(new_item);
            node_idx++;
        }

        for (boost::tie(ei, ei_end) = edges(dfg); ei != ei_end; ++ei) {
            auto curidx = dfg[*ei].param_idx;
            edgeparams.append(py::make_tuple(get(index_map, source(*ei, dfg)),
                                             get(index_map, target(*ei, dfg)),
                                             py::make_tuple(dofs.side_container[curidx].first,
                                                            dofs.bond_container[curidx].first,
                                                            dofs.bond_container[curidx].second,
                                                            dofs.side_container[curidx].second)));
        }

        log_routine("frags_done");
        log_in("file", this->filename);
        log_check("nodetypes_SET1", nodetypes);
        log_check("nodeatoms_SET1", nodeatoms);
        log_check("edgeparams_SET2", edgeparams);
        log_finish();
    #endif

    int free_count = 0;
    for (int i = 0; i < dofs.size(); ++i) {
        if (dofs.doftype_container[i] == dofs_t::free)
            dofs.free_indexer[i] = free_count++;
    }

    todo_queue.reserve(number_of_fragments);
    mol_data.frame_container = mol_data_t::frame_container_t(frag_comp.size());
    mol_data.frame_container[this->headnode].assign(boost::numeric::ublas::identity_matrix
                            <typename mol_data_t::frame_container_t::value_type::value_type>
                                (mol_data.frame_container[this->headnode].size1()));
    
    boost::graph_traits <dfg_t>::vertex_iterator vert_i, vert_end;
    boost::graph_traits<dfg_t>::in_edge_iterator in, in_end;
    boost::property_map <dfg_t, boost::vertex_index_t>::type dfg_index_map = get(boost::vertex_index, dfg);
    mol_data.frag_neighbors = mol_data_t::frag_neighbors_t(number_of_fragments);
    for (boost::tie(vert_i, vert_end) = vertices(dfg); vert_i != vert_end; ++vert_i) {
        const auto frag_idx = get(dfg_index_map, *vert_i); 
        for (boost::tie(in, in_end) = in_edges(*vert_i, dfg); in != in_end; ++in) {
            mol_data.frag_neighbors[frag_idx].push_back(get(dfg_index_map, source(*in, dfg)));
        }
        mol_data.frag_neighbors[frag_idx].shrink_to_fit();
    }

    checkpoint("mol_checkfree");
    #ifdef KDMOL_LOG
        {
        py::list log_free_dofs;
        for(const auto [dof_idx, _] : dofs.free_indexer) {
            auto dof_tuple = py::cast(std::make_tuple(
                dofs.side_container[dof_idx].first,
                dofs.bond_container[dof_idx].first,
                dofs.bond_container[dof_idx].second,
                dofs.side_container[dof_idx].second
            ));
            log_free_dofs.append(dof_tuple);
        }

        log_routine("mol_checkfree");
        log_in("file", this->filename);
        log_check("freedofs_SET4", log_free_dofs);
        log_finish();
        }
    #endif

    for (auto& cp : cyclic_parts)
        cp.index_owned_dofs();

    noncyclic_dof_idxs = active_idxs_t();    
    for (int i = 0; i < dofs.noncyclic_size; ++i)
        if (dofs.doftype_container[i] == dofs_t::free)
            noncyclic_dof_idxs.push_back(dofs.free_indexer[i]);
    noncyclic_dof_idxs.shrink_to_fit();

    #ifdef MOL_VALIDATE
        main_validator = validator_t(mol_data, dofs);
    #endif

    #ifdef MOL_OVERLAP_DETECTION
        overlap_detector = overlap_detector_t(mol_data.molgr, mol_data.xyz, mol_data.atom_symbols, "molecule");
    #endif
}

std::pair<std::vector<Utils::fulldof_t>, double*> Molecule::get_ps() {
    dofs.dof_vector = dofs_t::dof_vector_t(dofs.free_indexer.size());
    if (dofs.free_indexer.size() > 0)
        dofs.dof_rawpointer = &dofs.dof_vector[0];
    else
        dofs.dof_rawpointer = nullptr;

    for (const auto [ full_idx, free_idx ] : dofs.free_indexer)
        dofs.dof_vector.at(free_idx) = dofs.dofvalue_container[full_idx];
    
    std::vector<Utils::fulldof_t> dof_list;
    dof_list.reserve(dofs.free_indexer.size());
    for (const auto [ i, _ ] : dofs.free_indexer)
        dof_list.push_back({dofs.side_container[i].first,
                            dofs.bond_container[i].first,
                            dofs.bond_container[i].second,
                            dofs.side_container[i].second});
    return {dof_list, dofs.dof_rawpointer};
}

py::tuple Molecule::get_ps_py() {
    py::list dof_list;
    for (const auto& p : dofs.free_indexer) {
        dof_list.append(py::make_tuple(dofs.side_container[p.first].first,
                                       dofs.bond_container[p.first].first,
                                       dofs.bond_container[p.first].second,
                                       dofs.side_container[p.first].second));
    }

    const size_t size = dof_list.size();
    dofs.dof_rawpointer = new std::remove_pointer<Dofs::dof_rawpointer_t>::type[size];
    for (const auto [ full_idx, free_idx ] : dofs.free_indexer)
        dofs.dof_rawpointer[free_idx] = dofs.dofvalue_container[full_idx];

    py::capsule free_when_done(dofs.dof_rawpointer, [](void *f) {
        double *foo = reinterpret_cast<dofs_t::dof_rawpointer_t>(f);
        delete[] foo;
    });

    dofs.dof_nparray = dofs_t::dof_nparray_t(
        {static_cast<int>(size)}, // shape
        {sizeof(double)}, // C-style contiguous strides for double
        dofs.dof_rawpointer, // the data pointer
        free_when_done);
    return py::make_tuple(dof_list, dofs.dof_nparray);
}

std::pair<Molecule::owner_container_t, int*> Molecule::get_discrete_ps() {
    discrete_dofs.dof_vector = discrete_dofs_t::dof_vector_t(discrete_dofs.size());
    if (discrete_dofs.size() > 0)
        discrete_dofs.dof_rawpointer = &discrete_dofs.dof_vector[0];
    else
        discrete_dofs.dof_rawpointer = nullptr;

    for (int i = 0; i < discrete_dofs.size(); ++i)
        discrete_dofs.dof_rawpointer[i] = discrete_dofs.dofvalue_container[i];
    
    owner_container_t owner_container = discrete_dofs.owner_container;
    return {owner_container, discrete_dofs.dof_rawpointer};
}

py::array_t<double> Molecule::get_xyz_py() const {
    auto& xyz = mol_data.xyz.to_boost_format();

    const size_t size = xyz.size1() * xyz.size2();
    assertm(xyz.size2() == 3, "Unexpected behavior");
    double *coord_array = new double[size];
    
    for (int i = 0; i < xyz.size1(); ++i)
        for (int j = 0; j < 3; ++j)
            coord_array[i * 3 + j] = xyz(i, j);

    py::capsule free_when_done(coord_array, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        delete[] foo;
    });
    
    return py::array_t<double>(
        {static_cast<int>(xyz.size1()), static_cast<int>(xyz.size2())}, // shape
        {static_cast<int>(xyz.size2()) * 8, 8}, // C-style contiguous strides for double
        coord_array, // the data pointer
        free_when_done);
}

py::tuple Molecule::get_discrete_ps_py() {
    const size_t size = discrete_dofs.size();
    discrete_dofs.dof_rawpointer = new std::remove_pointer<discrete_dofs_t::dof_rawpointer_t>::type[size];

    py::capsule free_when_done(discrete_dofs.dof_rawpointer, [](void *f) {
        int *foo = reinterpret_cast<discrete_dofs_t::dof_rawpointer_t>(f);
        delete[] foo;
    });

    discrete_dofs.dof_nparray = discrete_dofs_t::dof_nparray_t(
        {static_cast<int>(size)}, // shape
        {sizeof(int)}, // C-style contiguous strides for int
        discrete_dofs.dof_rawpointer, // the data pointer
        free_when_done);
    return py::make_tuple(py::cast(discrete_dofs.owner_container), discrete_dofs.dof_nparray);
}

aps_return_t Molecule::apply_ps() {
    checkpoint("APS copy input values");
    assertm(dofs.pointer_is_good(), "Array of continuous DOFs was not initialized. Call get_ps() first.");
    for (const auto& [dof_idx, free_idx] : dofs.free_indexer)
        dofs.dofvalue_container[dof_idx] = dofs.dof_rawpointer[free_idx];

    checkpoint("APS copy discrete values");
    assertm(discrete_dofs.pointer_is_good(), "Array of discrete DOFs was not initialized. Call get_discrete_ps() first.");
    for (int i = 0; i < discrete_dofs.size(); ++i)
        discrete_dofs.dofvalue_container[i] = discrete_dofs.dof_rawpointer[i];

    checkpoint("aps_enter");
    #ifdef KDMOL_LOG
        Logger::aps_count++;

        py::list input_data, ddof_data;
        for (const auto& [dof_idx, free_idx] : dofs.free_indexer)
            input_data.append(py::make_tuple(dofs.side_container[dof_idx].first,
                                             dofs.bond_container[dof_idx].first,
                                             dofs.bond_container[dof_idx].second,
                                             dofs.side_container[dof_idx].second,
                                             dofs.dofvalue_container[dof_idx]));
        
        for (int i = 0; i < discrete_dofs.size(); ++i) {
            py::list newitem, cur_atoms;
            cur_atoms = py::handle(py::cast(discrete_dofs.owner_container[i])).cast<py::list>();
            cur_atoms.attr("sort")();
            newitem.append(cur_atoms);
            newitem.append(discrete_dofs.dofvalue_container[i]);
            ddof_data.append(newitem);
        }
        ddof_data.attr("sort")();

        log_routine("aps_enter");
        log_in("call_idx", Logger::aps_count);
        log_check("param_data_SET4", input_data);
        log_check("ddof_data", ddof_data);
        log_finish();
    #endif
    
    // assertm(mol_data.refiners.get_num_active() != 0, "apply_ps was called but there are no DOFs to refine");

    for (auto& cp : cyclic_parts) {
        auto result = cp.apply_ps();
        if (result != aps_return_t::success)
            return result;
    }
        
    for (auto& cp : cyclic_parts)
        cp.update_configuration(dfg);
    
    #ifdef KDMOL_LOG
        py::list nbs;
        for (int i = 0; i < mol_data.frag_neighbors.size(); ++i) {
            py::list newitem;
            newitem.append(i);
            newitem.append(py::cast(mol_data.frag_neighbors[i]));
            newitem[1].attr("sort")();
            nbs.append(newitem);
        }

        // log_routine("aps_nb_check");
        // log_in("call_idx", Logger::aps_count);
        // log_check("nbs_SET1", nbs);
        // log_finish();
    #endif

    using namespace boost::numeric::ublas;
    todo_queue.clear();
    todo_queue.push_back(this->headnode); // OPT: Precompute the queue
    Utils::uvector4d_t temp_coord;
    se_matrix_t temp_transform, temp_tmat;

    while (todo_queue.size() != 0) {
        auto& cur_frag = todo_queue.back();

        for (const auto& atom : mol_data.fg_container[cur_frag]) {
            noalias(temp_coord) = prod(mol_data.frame_container[cur_frag], mol_data.local_xyz.get_4d(atom));
            for (int i = 0; i < 3; ++i)
                mol_data.xyz[atom](i) = temp_coord(i);
            
            #ifdef KDMOL_LOG
                int fragtype;
                if (dynamic_cast<cyclic_part_t*>(frag_refs[cur_frag]) != nullptr)
                    fragtype = 1;
                else
                    fragtype = 0;

                // log_routine("aps_xyzcheck");
                // log_in("call_idx", Logger::aps_count);
                // log_in("curfrag", cur_frag);
                // log_in("atom", atom);
                // log_check("fragtype", as_list(fragtype));
                // log_check("frame", mol_data.frame_container[cur_frag]);
                // log_check("local_xyz", mol_data.local_xyz.get_4d(atom));
                // log_check("res_xyz", mol_data.xyz[atom]);
                // log_finish();
            #endif
        }

        const auto& dep_frags = mol_data.frag_neighbors[cur_frag];
        for (const auto& nb_frag : dep_frags) {
            #ifdef KDMOL_LOG
                if (!edge(nb_frag, cur_frag, dfg).second)
                    throw std::runtime_error(fmt::format("Edge {}-{} not found", cur_frag, nb_frag));
            #endif

            const auto param_idx = dfg[edge(nb_frag, cur_frag, dfg).first].param_idx;
            noalias(temp_transform) = prod(mol_data.frame_container[cur_frag], mol_data.transition_map[param_idx]);
            auto dof_value = dofs.get_tmat(param_idx, temp_tmat, 1);
            noalias(mol_data.frame_container[nb_frag]) = prod(temp_transform, temp_tmat);

            #ifdef KDMOL_LOG
                int fragtype, nb_fragtype;
                if (dynamic_cast<cyclic_part_t*>(frag_refs[cur_frag]) != nullptr)
                    fragtype = 1;
                else
                    fragtype = 0;
                if (dynamic_cast<cyclic_part_t*>(frag_refs[nb_frag]) != nullptr)
                    nb_fragtype = 1;
                else
                    nb_fragtype = 0;

                // log_routine("aps_framecheck");
                // log_in("call_idx", Logger::aps_count);
                // log_in("curfrag", cur_frag);
                // log_in("nb_frag", nb_frag);
                // log_check("fragtype", as_list(fragtype));
                // log_check("nb_fragtype", as_list(nb_fragtype));
                // log_check("start_frame", mol_data.frame_container[cur_frag]);
                // log_check("transition", mol_data.transition_map[param_idx]);
                // log_check("dof_value", as_list(dof_value));
                // log_check("t_mat", temp_tmat);
                // log_check("res_frame", mol_data.frame_container[nb_frag]);
                // log_finish();
            #endif
        }

        #ifdef KDMOL_LOG
            if (todo_queue.capacity() != mol_data.fg_container.size())
                throw std::runtime_error(fmt::format("Unexpected capacity of 'todo_queue'. Actual={} but should be {}", 
                                                                    todo_queue.capacity(), mol_data.fg_container.size()));
            
            py::list log_xyzs, log_frames;
            for (const auto& atom : mol_data.fg_container[cur_frag]) {
                py::list newitem;
                newitem.append(atom);
                newitem.append(Utils::vector_to_list(mol_data.xyz[atom]));
                log_xyzs.append(newitem);
            }

            for (const auto& nb_frag : dep_frags) {
                py::list newitem;
                newitem.append(nb_frag);
                newitem.append(Utils::matrix_to_list(mol_data.frame_container[nb_frag]));
                log_frames.append(newitem);
            }

            // log_routine("aps_reconstr_iter");
            // log_in("call_idx", Logger::aps_count);
            // log_in("curfrag", cur_frag);
            // log_check("xyzs_SET1", log_xyzs);
            // log_check("depframes_SET1", log_frames);
            // log_finish();
        #endif

        todo_queue.pop_back();
        todo_queue.insert(todo_queue.end(), dep_frags.cbegin(), dep_frags.cend());
    }

    checkpoint("aps_finish");
    #ifdef KDMOL_LOG
        {
        py::list input_data, log_xyzs, ddof_data;
        for (const auto& [dof_idx, free_idx] : dofs.free_indexer)
            input_data.append(py::make_tuple(dofs.side_container[dof_idx].first,
                                             dofs.bond_container[dof_idx].first,
                                             dofs.bond_container[dof_idx].second,
                                             dofs.side_container[dof_idx].second,
                                             dofs.dofvalue_container[dof_idx]));
        
        for (int atom = 0; atom < mol_data.xyz.size(); ++atom) {
            auto xyz = mol_data.xyz[atom];
            py::list newitem;
            newitem.append(atom);
            newitem.append(Utils::vector_to_list(xyz));
            log_xyzs.append(newitem);
        }

        for (int i = 0; i < discrete_dofs.size(); ++i) {
            py::list newitem, cur_atoms;
            cur_atoms = py::handle(py::cast(discrete_dofs.owner_container[i])).cast<py::list>();
            cur_atoms.attr("sort")();
            newitem.append(cur_atoms);
            newitem.append(discrete_dofs.dofvalue_container[i]);
            ddof_data.append(newitem);
        }
        ddof_data.attr("sort")();

        log_routine("aps_finish");
        log_in("file", this->filename);
        log_in("call_idx", Logger::aps_count);
        log_check("param_data_SET4", input_data);
        log_check("ddof_data", ddof_data);
        log_check("coords_SET1", log_xyzs);
        log_finish();
        }
    #endif

    #ifdef MOL_VALIDATE
        if (!main_validator.validate(mol_data.xyz, mol_data, dofs)) {
            mol_data.refiners.set_active(noncyclic_dof_idxs, {});
            return aps_return_t::validationfail;
            // throw Utils::GeomFailure(fmt::format("Validation failed"));
        }
    #endif

    #ifdef MOL_OVERLAP_DETECTION
        if (!overlap_detector(mol_data.xyz)) {
            mol_data.refiners.set_active(noncyclic_dof_idxs, {});
            return aps_return_t::overlap;
            // throw Utils::OverlapDetected(OVERLAP_MOLECULE_MESSAGE);
        }
    #endif

    // Actions on successful termination
    mol_data.refiners.set_successful();
    return aps_return_t::success;
}

void Molecule::reconfigure(const py::dict& upd_data) {
    std::set<int> reconfigure_frags;

    if (upd_data.attr("__contains__")("bonds").cast<bool>()) {
        const auto bond_items = upd_data["bonds"]["items"].cast<std::vector<std::pair<int, int>>>();
        const auto bond_values = upd_data["bonds"]["values"].cast<std::vector<double>>();
        if (bond_items.size() != bond_values.size())
            throw std::runtime_error(fmt::format("Lists of bond items and values have different sizes: {} vs. {}", 
                                            bond_items.size(), bond_values.size()));
        
        const auto nitems = bond_items.size();
        for (int i = 0; i < nitems; ++i) {
            const auto& bond_data = mol_data.bond_data(bond_items[i]);
            mol_data.bondlengths[bond_data.index] = bond_values[i];
            if ((reconfigure_frags.find(bond_data.parent_frag) == reconfigure_frags.end()) &&
                    (dynamic_cast<cyclic_part_t*>(frag_refs[bond_data.parent_frag]) == nullptr)) { // TODO: Should we really ignore cyc. parts?
                reconfigure_frags.insert(bond_data.parent_frag);
            }
        }
    }
    
    if (upd_data.attr("__contains__")("vangles").cast<bool>()) {
        const auto vangle_items = upd_data["vangles"]["items"].cast<std::vector<std::tuple<int, int, int>>>();
        const auto vangle_values = upd_data["vangles"]["values"].cast<std::vector<double>>();
        if (vangle_items.size() != vangle_values.size())
            throw std::runtime_error(fmt::format("Lists of vangle items and values have different sizes: {} vs. {}", 
                                            vangle_items.size(), vangle_values.size()));
        
        const auto nitems = vangle_items.size();
        std::unordered_map<int, std::vector<std::tuple<int, int, int, double>> > map_items;
        for (int i = 0; i < nitems; ++i) {
            const auto [ left_idx, central_idx, right_idx ] = vangle_items[i];
            map_items[central_idx].push_back(std::make_tuple(left_idx, central_idx, right_idx, vangle_values[i]));
        }

        for (const auto& [central_idx, center_data] : map_items) {
            std::vector<std::pair<int, int>> cur_vangle_items;
            std::vector<double> cur_vangle_values;
            for (auto [left_idx, central_idx, right_idx, vangle_value] : center_data) {
                cur_vangle_items.push_back(std::make_pair(left_idx, right_idx));
                cur_vangle_values.push_back(vangle_value);
            }

            mol_data.polys[central_idx].change_vangles(cur_vangle_items, cur_vangle_values);

            const auto parent_frag = mol_data.atom_in_frag[central_idx];
            if ((reconfigure_frags.find(parent_frag) == reconfigure_frags.end()) &&
                    (dynamic_cast<cyclic_part_t*>(frag_refs[parent_frag]) == nullptr)) { // TODO: Should we really ignore cyc. parts?
                reconfigure_frags.insert(parent_frag);
            }
        }
    }

    for (const auto& frag_idx : reconfigure_frags)
        frag_refs[frag_idx]->configure_assembly(mol_data.molgr, dfg);
}

py::list Molecule::get_biggest_ringfrag_atoms() const {
    if (mol_data.fg_container.size() == 0)
        return py::list();
    auto max_atoms = mol_data.fg_container[0].size();
    auto res_idx = 0;
    for (int i = 1; i < mol_data.fg_container.size(); ++i) {
        const auto& atoms = mol_data.fg_container[i];
        if (atoms.size() > max_atoms) {
            res_idx = i;
            max_atoms = atoms.size();
        }
    }
    return py::cast(std::vector<int>(mol_data.fg_container[res_idx].cbegin(), mol_data.fg_container[res_idx].cend()));
}

int Molecule::get_num_flexible_rings() const {
    int result = 0;
    for (const auto& cp : cyclic_parts)
        result += cp.get_num_tlcsolvers();
    return result;
}

int Molecule::get_num_rigid_rings() const {
    int result = 0;
    for (const auto& cp : cyclic_parts)
        result += cp.get_num_idsolvers();
    return result;
}
