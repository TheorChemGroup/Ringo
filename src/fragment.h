#ifndef FRAGMENT_H
#define FRAGMENT_H

#include "generic_fragment.h"
#include "utils.h"


template <class G, class MD, class DF, class DD, class D>
class Fragment : public GenericFragment<G, MD, DF, DD, D> {
    private:
        using parent_t = GenericFragment<G, MD, DF, DD, D>;
        int index;
        #ifdef KDMOL_LOG
            int call_idx;
        #endif

    public:
        using main_graph_t = G;
        using mol_data_t = MD;
        using dofs_t = DF;
        using discrete_dofs_t = DD;
        using dfg_t = D;

        using dir_g_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
        using point_t = typename mol_data_t::poly_t::point_t;

        Fragment() { }

        Fragment(mol_data_t* mol_data, dofs_t* dofs, discrete_dofs_t* discrete_dofs, const int i) : 
                    parent_t(mol_data, dofs, discrete_dofs), 
                    index(i)
        {
            #ifdef KDMOL_LOG
                call_idx = 0;
            #endif
        }

        void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg, const int INDEX) override;

        void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg) override
        { this->configure_assembly(molgr, dfg, index); }
};


template <class G, class MD, class DF, class DD, class D>
void Fragment<G, MD, DF, DD, D>::configure_assembly(const main_graph_t& molgr, const dfg_t& dfg, const int INDEX) {
    #ifdef KDMOL_LOG
        call_idx++;
    #endif
    if (INDEX != this->index)
        throw std::runtime_error(fmt::format("Unexpected index {} vs. {}", INDEX, this->index));
    auto myvertex = vertex(this->index, dfg);
    int out_count = 0;
    typename boost::graph_traits<dfg_t>::out_edge_iterator out, out_end;
    typename dfg_t::edge_descriptor out_edge;
    for (boost::tie(out, out_end) = out_edges(myvertex, dfg); out != out_end; ++out) {
        out_count++;
        out_edge = *out;
    }

    int any_atom, start_atom;
    if (out_count == 1) {
        const auto bond = this->dofs->bond_container[dfg[out_edge].param_idx];
        if (this->mol_data->belongs_to_fragment(bond.second, this->index))
            start_atom = bond.second;
        else if (this->mol_data->belongs_to_fragment(bond.first, this->index))
            start_atom = bond.first;
        else
            throw std::runtime_error(fmt::format("Not either atom of the bond {}-{} belongs to atoms of the fragment {}", bond.first, bond.second, this->index));
        any_atom = 0;
    } else { // if out_count == 0
        const auto& fragment_atoms = this->mol_data->fg_container[this->index];
        start_atom = *(fragment_atoms.begin());
        any_atom = 1;
    }

    std::set<int> end_atoms;
    typename boost::property_map <dfg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dfg);
    typename boost::graph_traits<dfg_t>::in_edge_iterator in, in_end;
    for (boost::tie(in, in_end) = in_edges(myvertex, dfg); in != in_end; ++in) {
        auto in_edge = *in;
        auto in_node_idx = get(index_map, source(in_edge, dfg));
        const auto bond = this->dofs->bond_container[dfg[in_edge].param_idx];
        if (this->mol_data->belongs_to_fragment(bond.second, this->index))
            end_atoms.insert(bond.first);
        else if (this->mol_data->belongs_to_fragment(bond.first, this->index))
            end_atoms.insert(bond.second);
        else
            throw std::runtime_error(fmt::format("Not either atom of the bond {}-{} belongs to atoms of the fragment {}", bond.first, bond.second, this->index));
    }

    // #ifdef KDMOL_LOG
    //     {
    //     py::list end_atoms_list;
    //     for (const auto& atom : end_atoms)
    //         end_atoms_list.append(atom);

    //     log_routine("frag_intconfigA");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("start_atom", as_list(start_atom));
    //     log_check("any_atom", as_list(any_atom));
    //     log_check("end_atoms_SET1", end_atoms_list);
    //     log_check("startatom", as_list(start_atom));
    //     log_finish();
    //     }
    // #endif

    int dir_g_size;
    std::map<int, int> dirg_idxs;// normal -> dir_g
    std::vector<int> dirg_atoms; // dir_g -> normal
    { // This part only adds some elements to end_atoms and computes the size of dir_g
        py::list subgraph_nodes;
        for (const auto& atom : this->mol_data->fg_container[this->index])
            subgraph_nodes.append(atom);
        for (const auto& atom : end_atoms)
            subgraph_nodes.append(atom);
        auto undir_g = molgr.attr("subgraph")(subgraph_nodes);
        for (const auto& atom : undir_g.attr("nodes")()) {
            py::int_ atom_idx = py::handle(atom).cast<py::int_>();
            if ((py::list(undir_g.attr("neighbors")(atom_idx)).size() == 1) && (atom_idx.cast<int>() != start_atom))
                end_atoms.insert(atom_idx.cast<int>());
        }

        dir_g_size = py::handle(undir_g).attr("number_of_nodes")().cast<int>();
        int k = 0;
        for (const auto& atom : py::list(undir_g.attr("nodes")())) {
            dirg_atoms.push_back(py::handle(atom).cast<int>());
            dirg_idxs[py::handle(atom).cast<int>()] = k++;
        }
        dirg_atoms.shrink_to_fit();

        // #ifdef KDMOL_LOG
        //     {
        //     py::list end_atoms_list;
        //     for (const auto& atom : end_atoms)
        //         end_atoms_list.append(atom);

        //     log_routine("frag_intconfigB");
        //     log_in("idx", this->index);
        //     log_in("callidx", this->call_idx);
        //     log_check("subgraph_nodes_SET1", subgraph_nodes);
        //     log_check("end_atoms_SET1", end_atoms_list);
        //     log_finish();
        //     }
        // #endif
    }

    // Create in undirected Python prototype of 
    py::list subgraph_nodes;
    for (const auto& atom : this->mol_data->fg_container[this->index])
        subgraph_nodes.append(atom);
    for (const auto& atom : end_atoms)
        subgraph_nodes.append(atom);
    for (const auto& [atom, _] : this->mol_data->polys[start_atom])
        subgraph_nodes.append(atom);
    auto undir_g = molgr.attr("subgraph")(subgraph_nodes); // caution: subgraph_nodes contains duplicates
                                                            // OPT: check if undir_g is not the same as undir_g generated above 
    // #ifdef KDMOL_LOG
    //     log_routine("frag_intconfigC");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("subgraph_nodes_SET1", py::list(undir_g.attr("nodes")()));
    //     log_finish();
    // #endif

    dir_g_t dir_g(dirg_atoms.size());
    py::int_ start_atom_py = py::cast(start_atom);
    for (const auto& enode : end_atoms) {
        auto solvepath = py::list(Utils::PyModules::get_nx()->attr("shortest_path")(undir_g, start_atom_py, py::cast(enode)));
        for (int i = 0; i < solvepath.size() - 1; ++i) {
            auto vA = vertex(dirg_idxs[solvepath[i].cast<int>()], dir_g),
                    vB = vertex(dirg_idxs[solvepath[i + 1].cast<int>()], dir_g);
            if (!edge(vA, vB, dir_g).second)
                add_edge(vA, vB, dir_g);
        }
    }

    // #ifdef KDMOL_LOG
    //     {
    //     py::list dfg_nodes, dfg_edges;
    //     boost::graph_traits <dir_g_t>::vertex_iterator i, end;
    //     boost::property_map <dir_g_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dir_g);

    //     for (boost::tie(i, end) = vertices(dir_g); i != end; ++i) {
    //         dfg_nodes.append(dirg_atoms[get(index_map, *i)]);
    //     }
        
    //     boost::graph_traits<dir_g_t>::edge_iterator ei, ei_end;
    //     for (boost::tie(ei, ei_end) = edges(dir_g); ei != ei_end; ++ei)
    //         dfg_edges.append(py::make_tuple(dirg_atoms[get(index_map, source(*ei, dir_g))],
    //                                         dirg_atoms[get(index_map, target(*ei, dir_g))]));

    //     log_routine("frag_intconfigD");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("undir_size", as_list(py::handle(undir_g).attr("number_of_nodes")().cast<int>()));
    //     log_check("dir_size", as_list(dirg_atoms.size()));
    //     log_check("node_data_SET1", dfg_nodes);
    //     log_check("edge_data_SET2", dfg_edges);
    //     log_finish();
    //     }
    // #endif

    Utils::MapXyzContainer<std::map<int, point_t>> frag_xyz;
    #ifdef KDMOL_LOG
        std::set<int> done_ones;
        done_ones.insert(start_atom);
        py::list bondlengths;
    #endif
    frag_xyz.set_atom(start_atom, {0.0, 0.0, 0.0});
    for (const auto& [atom, xyz] : this->mol_data->polys[start_atom]) { // TODO Check if I didn't use str.bindings somewhere
        if (atom == start_atom)
            continue;
        auto bond_length = this->mol_data->bondlengths[this->mol_data->bond_data({start_atom, atom}).index];
        frag_xyz[atom];
        for (int i = 0; i < 3; ++i)
            frag_xyz[atom](i) = xyz(i) * bond_length;
        #ifdef KDMOL_LOG
            done_ones.insert(atom);
            py::list lengthitem;
            lengthitem.append(atom);
            lengthitem.append(bond_length);
            bondlengths.append(lengthitem);
        #endif
    }
    
    // #ifdef KDMOL_LOG
    //     {
    //     py::list done_ones_py, coords_log;
    //     for (const auto& atom : done_ones)
    //         done_ones_py.append(py::cast(atom));
        
    //     for (const auto& [atom, xyz] : frag_xyz) {
    //         py::list cur;
    //         cur.append(atom);
    //         cur.append(xyz(0));
    //         cur.append(xyz(1));
    //         cur.append(xyz(2));
    //         coords_log.append(cur);
    //     }

    //     log_routine("frag_intconfigE");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("bondlengths_SET1", bondlengths);
    //     log_check("done_ones_SET1", done_ones_py);
    //     log_check("coords_SET1", coords_log);
    //     log_finish();
    //     }
    // #endif

    #ifdef KDMOL_LOG
        int step = 0;
    #endif
    std::vector<std::pair<int, int>> todo_atoms;
    todo_atoms.reserve(dirg_atoms.size());
    typename boost::property_map <dir_g_t, boost::vertex_index_t>::type dirg_index_map = get(boost::vertex_index, dir_g);
    typename boost::graph_traits<dir_g_t>::out_edge_iterator dirg_out, dirg_out_end;
    for (boost::tie(dirg_out, dirg_out_end) = out_edges(vertex(dirg_idxs[start_atom], dir_g), dir_g); dirg_out != dirg_out_end; ++dirg_out)
        todo_atoms.push_back({start_atom, dirg_atoms[get(dirg_index_map, target(*dirg_out, dir_g))]});
    py::list todo_atoms_py = py::cast(todo_atoms);
    todo_atoms_py.attr("sort")(); // OPT: Why sort??
    todo_atoms = todo_atoms_py.cast<decltype(todo_atoms)>();
    // #ifdef KDMOL_LOG
    //     log_routine("frag_intconfigX");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_in("step", step);
    //     log_check("todoatoms", todo_atoms);
    //     log_finish();
    // #endif

    #ifdef KDMOL_LOG
        int iter = 0;
    #endif
    while (todo_atoms.size() != 0) {
        auto [prev_atom, cur_atom] = todo_atoms.back();

        // #ifdef KDMOL_LOG
        //     iter++;
        //     step++;
        //     log_routine("frag_intIter");
        //     log_in("idx", this->index);
        //     log_in("callidx", this->call_idx);
        //     log_in("iter", iter);
        //     log_check("prev_atom", as_list(prev_atom));
        //     log_check("cur_atom", as_list(cur_atom));
        //     log_check("todo", todo_atoms);
        //     log_finish();
        // #endif
        
        todo_atoms.pop_back();
        if (end_atoms.find(cur_atom) != end_atoms.end())
            continue;

        if (this->mol_data->number_of_neighbors(prev_atom) > 1) {
            const auto param_idx = this->dofs->param_from_bond({prev_atom, cur_atom});
            const auto& bond_atoms = this->dofs->bond_container.at(param_idx);
            const auto& side_atoms = this->dofs->side_container.at(param_idx);
            int past_side, next_side;
            if (bond_atoms.second == cur_atom) {
                past_side = side_atoms.first;
                next_side = side_atoms.second;
            } else { // if bond_atoms.first == cur_atom
                past_side = side_atoms.second;
                next_side = side_atoms.first;
            }

            auto startframe = this->mol_data->polys[prev_atom].get_frame(frag_xyz,
                                            {cur_atom, past_side}, cur_atom);
            for (int i = 0; i < 3; ++i) {
                startframe(i, 0) = -startframe(i, 0);
                startframe(i, 2) = -startframe(i, 2);
            }
            
            se_matrix_t temp_tmat;
            auto ang = this->dofs->get_tmat(param_idx, temp_tmat, -1);
            startframe = prod(startframe, temp_tmat);
            this->mol_data->polys[cur_atom].restore_from_frame(frag_xyz, {prev_atom, next_side}, this->mol_data, startframe);

            // #ifdef KDMOL_LOG
            //     {
            //     log_routine("frag_intconfigG");
            //     log_in("idx", this->index);
            //     log_in("callidx", this->call_idx);
            //     log_in("step", step);
            //     log_in("iter", iter);
            //     log_check("cur_atom", as_list(cur_atom));
            //     log_check("past_side", as_list(past_side));
            //     log_check("next_side", as_list(next_side));
            //     log_check("startframe", Utils::repr_matrix(startframe));
            //     log_check("ang", as_list(ang));
            //     log_check("tmat", Utils::repr_matrix(temp_tmat));
            //     log_finish();
            //     }
            // #endif
        } else {
            auto startframe = this->mol_data->polys[prev_atom].get_frame(frag_xyz, cur_atom, cur_atom);
            for (int i = 0; i < 3; ++i) {
                startframe(i, 0) = -startframe(i, 0);
                startframe(i, 2) = -startframe(i, 2);
            }
            this->mol_data->polys[cur_atom].restore_from_frame(frag_xyz, prev_atom, this->mol_data, startframe);

            // #ifdef KDMOL_LOG
            //     {
            //     log_routine("frag_intconfigG");
            //     log_in("idx", this->index);
            //     log_in("callidx", this->call_idx);
            //     log_in("step", step);
            //     log_in("iter", iter);
            //     log_check("prev_atom", as_list(prev_atom));
            //     log_check("cur_atom", as_list(cur_atom));
            //     log_check("startframe", startframe);
            //     log_finish();
            //     }
            // #endif
        }

        for (boost::tie(dirg_out, dirg_out_end) = out_edges(vertex(dirg_idxs[cur_atom], dir_g), dir_g); 
                    dirg_out != dirg_out_end; ++dirg_out) {
            auto next_atom = dirg_atoms[get(dirg_index_map, target(*dirg_out, dir_g))];
            todo_atoms.push_back({cur_atom, next_atom});
            #ifdef KDMOL_LOG
                if (frag_xyz.find(next_atom) == frag_xyz.end())
                    throw std::runtime_error(fmt::format("The next atom={} was not computed. curr={}", next_atom, cur_atom));
                done_ones.insert(dirg_atoms[get(dirg_index_map, target(*dirg_out, dir_g))]);
            #endif
        }
        // #ifdef KDMOL_LOG
        //     {
        //     py::list coords_log;
        //     for (const auto& atom : done_ones) {
        //         py::list cur;
        //         cur.append(atom);
        //         cur.append(frag_xyz[atom](0));
        //         cur.append(frag_xyz[atom](1));
        //         cur.append(frag_xyz[atom](2));
        //         coords_log.append(cur);
        //     }
        //     log_routine("frag_intconfigF");
        //     log_in("idx", this->index);
        //     log_in("callidx", this->call_idx);
        //     log_in("step", step);
        //     log_in("iter", iter);
        //     log_check("prev_atom", as_list(prev_atom));
        //     log_check("cur_atom", as_list(cur_atom));
        //     log_check("coords_SET1", coords_log);
        //     log_finish();
        //     }
        // #endif
    }

    se_matrix_t startframe;
    if (out_count == 1) {
        auto param_idx = dfg[out_edge].param_idx;
        if (this->mol_data->belongs_to_fragment(this->dofs->bond_container[param_idx].second, this->index))
            startframe = this->mol_data->polys[this->dofs->bond_container[param_idx].second].get_frame(frag_xyz,
                                                {this->dofs->bond_container[param_idx].first,
                                                this->dofs->side_container[param_idx].second},
                                                this->dofs->bond_container[param_idx].second);
        else // if this->mol_data->belongs_to_fragment(this->dofs->bond_container[param_idx].first, this->index)
            startframe = this->mol_data->polys[this->dofs->bond_container[param_idx].first].get_frame(frag_xyz,
                                                {this->dofs->bond_container[param_idx].second,
                                                this->dofs->side_container[param_idx].first},
                                                this->dofs->bond_container[param_idx].first);
        for (int i = 0; i < 3; ++i) {
            startframe(i, 0) = -startframe(i, 0);
            startframe(i, 2) = -startframe(i, 2);
        }
    } else
        startframe.assign(boost::numeric::ublas::identity_matrix<se_matrix_t::value_type>(startframe.size1()));
    
    // #ifdef KDMOL_LOG
    //     {
    //     log_routine("frag_intFrame");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("frame", startframe);
    //     log_finish();
    //     }
    // #endif

    Utils::uvector4d_t temp_vec;
    temp_vec(3) = 1.0;
    se_matrix_t start_inv;
    Utils::invert_matrix(startframe, start_inv);
    auto& local_xyz = this->mol_data->local_xyz;
    #ifdef KDMOL_LOG
        py::list localcoords_log;
    #endif
    for (const auto& atom : this->mol_data->fg_container[this->index]) {
        for (int i = 0; i < 3; ++i)
            temp_vec(i) = frag_xyz[atom](i);
        noalias(local_xyz.get_4d(atom)) = prod(start_inv, temp_vec);
        #ifdef KDMOL_LOG
            py::list py_coords;
            py_coords.append(py::cast(atom));
            py_coords.append(py::cast(local_xyz.get_4d(atom)(0)));
            py_coords.append(py::cast(local_xyz.get_4d(atom)(1)));
            py_coords.append(py::cast(local_xyz.get_4d(atom)(2)));
            localcoords_log.append(py_coords);
        #endif
    }
    // #ifdef KDMOL_LOG
    //     {
    //     log_routine("frag_intLocalCoords");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("localcoords_SET1", localcoords_log);
    //     log_finish();
    //     }
    // #endif

    #ifdef KDMOL_LOG
        py::list trans_log;
    #endif
    for (boost::tie(in, in_end) = in_edges(myvertex, dfg); in != in_end; ++in) {
        auto in_edge = *in;
        auto in_node_idx = get(index_map, source(in_edge, dfg));
        int param_idx = dfg[in_edge].param_idx;
        se_matrix_t newframe;
        if (this->mol_data->belongs_to_fragment(this->dofs->bond_container[param_idx].second, this->index))
            newframe = this->mol_data->polys[this->dofs->bond_container[param_idx].second].get_frame(frag_xyz,
                                            {this->dofs->bond_container[param_idx].first,
                                                this->dofs->side_container[param_idx].second},
                                                this->dofs->bond_container[param_idx].first);
        else
            newframe = this->mol_data->polys[this->dofs->bond_container[param_idx].first].get_frame(frag_xyz,
                                            {this->dofs->bond_container[param_idx].second,
                                                this->dofs->side_container[param_idx].first},
                                                this->dofs->bond_container[param_idx].second);

        noalias(this->mol_data->transition_map[param_idx]) = prod(start_inv, newframe);
        #ifdef KDMOL_LOG
            py::list temp_matr;
            temp_matr.append(py::cast(target(in_edge, dfg)));
            temp_matr.append(py::cast(source(in_edge, dfg)));
            temp_matr.append(Utils::matrix_to_list(this->mol_data->transition_map[param_idx]));
            trans_log.append(temp_matr);
        #endif
    }
    // #ifdef KDMOL_LOG
    //     log_routine("frag_intTransitions");
    //     log_in("idx", this->index);
    //     log_in("callidx", this->call_idx);
    //     log_check("transitions_SET2", trans_log);
    //     log_finish();
    // #endif

    { // Assign bonds that belong to this fragment 
        boost::property_map <dir_g_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dir_g);
        boost::graph_traits<dir_g_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(dir_g); ei != ei_end; ++ei) {
            int vA = dirg_atoms[get(index_map, source(*ei, dir_g))],
                vB = dirg_atoms[get(index_map, target(*ei, dir_g))];
            this->mol_data->assign_parent_fragment({vA, vB}, this->index);
        }
    }
}

#endif // FRAGMENT_H
