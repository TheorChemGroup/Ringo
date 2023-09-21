#ifndef POLY_H
#define POLY_H

#include "utils.h"

using Utils::umatrix_t;
using Utils::uvector3d_t;
using Utils::uvector4d_t;


inline void get_dmat_sph(const double ang, so_matrix_t& res) {
    res.assign(boost::numeric::ublas::identity_matrix<typename so_matrix_t::value_type>(res.size1()));
    const double cos_ang = cos(ang);
    const double sin_ang = sin(ang);
    res(0, 0) = cos_ang; res(0, 1) = -sin_ang;
    res(1, 0) = sin_ang; res(1, 1) = cos_ang;
}

inline void get_vmat_sph(const double ang, so_matrix_t& res) {
    res.assign(boost::numeric::ublas::identity_matrix<typename so_matrix_t::value_type>(res.size1()));
    const double cos_ang = cos(ang);
    const double sin_ang = sin(ang);
    res(1, 1) = -cos_ang; res(1, 2) = -sin_ang;
    res(2, 1) = sin_ang; res(2, 2) = -cos_ang;
}

template <class T>
class Polyhedron {
    public:
        using point_t = T;
        using coordinates_t = Utils::MapXyzContainer<std::map<int, point_t>>;
        using dgr_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS>;
        
        Polyhedron() { }

        template<class main_graph_t, class coord_container_t>
        Polyhedron(const main_graph_t& molgr, coord_container_t& xyz_data, const int& central_idx);

        template<class passed_poly_t>
        Polyhedron(const passed_poly_t& molgr, const int& central_idx);

        template <class C>
        se_matrix_t get_frame(C& xyz_data, const std::pair<int, int>& nb_pair, const int& origin_idx);
        
        template <class C>
        se_matrix_t get_frame(C& xyz_data, const int& nb_atom, const int& origin_idx);

        template <class C>
        so_matrix_t get_frame_3d(C& xyz_data, const std::pair<int, int>& nb_pair);
        
        inline so_matrix_t get_localframe_3d(const std::pair<int, int>& nb_pair);

        template <class C, class M>
        void restore_neighbors(C& global_coords, const std::pair<int, int>& frame_atoms, const M* mol_data);
        
        template <class C, class M>
        void restore_from_frame(C& global_coords, const std::pair<int, int>& frame_atoms, const M* mol_data, const se_matrix_t& frame);
        
        template <class C, class M>
        void restore_from_frame(C& global_coords, const int& frame_atom, const M* mol_data, const se_matrix_t& frame);

        template < template<class ...> class PairContainer, class ... PairArgs,
                   template<class ...> class ValueContainer, class ... ValueArgs >
        void change_vangles(const PairContainer<std::pair<int, int>, PairArgs...>& vangle_items,
                            const ValueContainer<double, ValueArgs...>& vangle_values);

        template <class G, class C>
        inline void resolve_triangle(G& subgr, C& new_coords);
        
        template <class G, class C>
        inline void resolve_branched(G& subgr, C& new_coords);

        template <class C, template<class ...> class NodeContainer, class ... NodeArgs>
        inline void find_best_fit(C& new_coords, const NodeContainer<int, NodeArgs...>& nodes);

        inline double get_vangle(const int atA, const int atC) const noexcept
        { return local_coords.get_vangle(atA, central_idx, atC); }

        inline double get_vangle_rad(const int atA, const int atC) const noexcept
        { return get_vangle(atA, atC) * DEG2RAD; }

        inline double get_original_vangle(const int atA, const int atC) const noexcept
        { return coords_backup.get_vangle(atA, central_idx, atC); }

        inline double get_sph_angle(const int a, const int b, const int c)  noexcept;

        inline int get_center() const noexcept { return central_idx; }

        inline typename coordinates_t::iterator begin() noexcept { return local_coords.begin(); }
        inline typename coordinates_t::iterator end() noexcept { return local_coords.end(); }
        inline typename coordinates_t::const_iterator cbegin() const noexcept { return local_coords.cbegin(); }
        inline typename coordinates_t::const_iterator cend() const noexcept { return local_coords.cend(); }

    private:
        int central_idx;
        coordinates_t local_coords, coords_backup;
        #ifdef KDMOL_LOG
            std::unordered_map<std::string, int> call_idxs;
            std::string vangle_items, vangle_values;
        #endif
};

template <class T>
template<class main_graph_t, class coord_container_t>
Polyhedron<T>::Polyhedron(const main_graph_t& molgr, coord_container_t& xyz_data, const int& central_idx)
{
    this->central_idx = central_idx;
    point_t central_xyz = xyz_data.at(central_idx);

    for (const auto& py_idx : molgr.attr("neighbors")(central_idx)) {
        auto idx = py::handle(py_idx).cast<int>();
        point_t& nb_xyz = local_coords[idx];
        nb_xyz = xyz_data.at(idx);
        nb_xyz -= central_xyz;
        nb_xyz /= boost::numeric::ublas::norm_2(nb_xyz);
    }
    local_coords.set_atom(this->central_idx, {0.0, 0.0, 0.0});
    coords_backup = local_coords;

    #ifdef KDMOL_LOG
        call_idxs["polyrestore_singleA"] = 0;
        call_idxs["poly_frameS"] = 0;
        call_idxs["polyrestore_input"] = 0;

        py::list log_atoms, log_coords;
        for(const auto& [atom_idx, _] : local_coords)
            if (atom_idx != this->central_idx)
                log_atoms.append(atom_idx);
        
        log_atoms.attr("sort")();
        for(const auto& atom : log_atoms)
            log_coords.append(Utils::vector_to_list(local_coords[atom.cast<int>()]));

        // log_routine("poly_construct");
        // log_in("central_atom", central_idx);
        // log_check("atoms", log_atoms);
        // log_check("coords", log_coords);
        // log_finish();
    #endif
}

template <class T>
template<class passed_poly_t>
Polyhedron<T>::Polyhedron(const passed_poly_t& passed_poly, const int& central_idx) {
    this->central_idx = central_idx;

    auto idxs = py::handle(passed_poly.attr("idxs")).cast<std::vector<int>>();
    auto coords = py::handle(passed_poly.attr("coords")).cast<std::vector<std::vector<double>>>();
    assertm(idxs.size() == coords.size(), "List sizes are not equal");
    for (int i = 0; i < idxs.size(); ++i) {
        local_coords.set_atom(idxs[i], {coords[i][0], coords[i][1], coords[i][2]});
    }
    local_coords.set_atom(this->central_idx, {0.0, 0.0, 0.0});
    coords_backup = local_coords;
}

template <class T>
inline double Polyhedron<T>::get_sph_angle(const int a, const int b, const int c)  noexcept
{
    auto res = local_coords.get_dihedral(a, this->central_idx, b, c);
    checkpoint("get_sph_angle");
    #ifdef KDMOL_LOG
        log_routine("get_sph_angle");
        log_in("a", a);
        log_in("b", b);
        log_in("c", c);
        log_in("central_atom", this->central_idx);
        log_check("a_coord", local_coords.at(a));
        log_check("b_coord", local_coords.at(b));
        log_check("c_coord", local_coords.at(c));
        log_check("central_atom_coord", local_coords.at(this->central_idx));
        log_check("res", as_list(res));
        log_finish();
    #endif
    return res;
}

template <class T>
template <class C>
se_matrix_t Polyhedron<T>::get_frame(C& xyz_data, const std::pair<int, int>& nb_pair, const int& origin_idx) {
    const auto at0 = xyz_data.at(nb_pair.first);
    const auto at1 = xyz_data.at(central_idx);
    const auto at2 = xyz_data.at(nb_pair.second);
    
    uvector3d_t xv = at0 - at1;
    uvector3d_t yv = at2 - at1;
    uvector3d_t zv = Utils::gs_rand(xv, yv);
    
    se_matrix_t res;
    auto origin = xyz_data.at(origin_idx);
    for(int i = 0; i < 3; ++i) {
        res(i, 0) = xv(i);
        res(i, 1) = yv(i);
        res(i, 2) = zv(i);
        res(i, 3) = origin(i);
    }
    res(3, 3) = 1.0;

    // #ifdef KDMOL_LOG
    //     log_routine("poly_frame");
    //     log_in("central_atom", central_idx);
    //     log_in("atA", nb_pair.first);
    //     log_in("atB", nb_pair.second);
    //     log_in("at0", at0);
    //     log_in("at1", at1);
    //     log_in("at2", at2);
    //     log_check("xv", xv);
    //     log_check("yv", yv);
    //     log_check("zv", zv);
    //     log_finish();
    // #endif
    return res;
}

template <class T>
template <class C>
se_matrix_t Polyhedron<T>::get_frame(C& xyz_data, const int& nb_atom, const int& origin_idx) {
    const auto at0 = xyz_data.at(nb_atom);
    const auto at1 = xyz_data.at(central_idx);
    uvector3d_t at2;
    at2(0) = 1.0;
    at2(1) = 0.0;
    at2(2) = 0.0;
    
    uvector3d_t xv = at0 - at1;
    uvector3d_t yv = at2 - at1;
    uvector3d_t zv = Utils::gs_rand(xv, yv);
    
    se_matrix_t res;
    auto origin = xyz_data.at(origin_idx);
    for(int i = 0; i < 3; ++i) {
        res(i, 0) = xv(i);
        res(i, 1) = yv(i);
        res(i, 2) = zv(i);
        res(i, 3) = origin(i);
    }
    res(3, 3) = 1.0;
    // #ifdef KDMOL_LOG
    //     call_idxs["poly_frameS"] += 1;

    //     log_routine("poly_frameS");
    //     log_in("central_atom", central_idx);
    //     log_in("nbatom", nb_atom);
    //     log_in("call_idx", call_idxs["poly_frameS"]);
    //     log_check("at0", at0);
    //     log_check("at1", at1);
    //     log_check("at2", at2);
    //     log_check("xv", xv);
    //     log_check("yv", yv);
    //     log_check("zv", zv);
    //     log_finish();
    // #endif
    return res;
}

template <class T>
template <class C>
so_matrix_t Polyhedron<T>::get_frame_3d(C& xyz_data, const std::pair<int, int>& nb_pair) {
    // const auto at0 = xyz_data.at(nb_pair.first);
    // const auto at1 = xyz_data.at(central_idx);
    // const auto at2 = xyz_data.at(nb_pair.second);
    const auto at0 = Utils::at(xyz_data, nb_pair.first);
    const auto at1 = Utils::at(xyz_data, central_idx);
    const auto at2 = Utils::at(xyz_data, nb_pair.second);
    
    uvector3d_t xv = at0 - at1;
    uvector3d_t yv = at2 - at1;
    uvector3d_t zv = Utils::gs_rand(xv, yv);
    
    so_matrix_t res;
    for(int i = 0; i < 3; ++i) {
        res(i, 0) = xv(i);
        res(i, 1) = yv(i);
        res(i, 2) = zv(i);
    }
    return res;
}

template <class T>
inline so_matrix_t Polyhedron<T>::get_localframe_3d(const std::pair<int, int>& nb_pair) {
    return get_frame_3d(this->local_coords, nb_pair);
}

template <class T>
template <class C, class M>
void Polyhedron<T>::restore_neighbors(C& global_coords, const std::pair<int, int>& frame_atoms, const M* mol_data) {
    auto global_frame = get_frame_3d(global_coords, frame_atoms);
    auto local_frame = get_localframe_3d(frame_atoms);

    so_matrix_t transition_m, inv_local_frame;
    Utils::invert_matrix(local_frame, inv_local_frame);
    noalias(transition_m) = prod(global_frame, inv_local_frame); 
    
    // #ifdef KDMOL_LOG
    //     log_routine("poly_nbrestoreA");
    //     log_in("central_atom", central_idx);
    //     log_in("idxs", frame_atoms);
    //     log_check("transition", transition_m);
    //     log_finish();
    // #endif

    // #ifdef KDMOL_LOG
    //     py::list newcoords;
    // #endif
    for (const auto& [atom_idx, atom_xyz] : local_coords) {
        if ((atom_idx == this->central_idx) || (atom_idx == frame_atoms.first) || (atom_idx == frame_atoms.second))
            continue;
        noalias(global_coords.at(atom_idx)) = 
                mol_data->bondlengths[mol_data->bond_data({this->central_idx, atom_idx}).index] * 
                prod(transition_m, atom_xyz) + global_coords.at(this->central_idx);
        // #ifdef KDMOL_LOG
        //     py::list newitem;
        //     newitem.append(pair.first);
        //     newitem.append(Utils::vector_to_list(global_coords.at(pair.first)));
        //     newcoords.append(newitem);
        // #endif
    }
    
    // #ifdef KDMOL_LOG
    //     log_routine("poly_nbrestoreB");
    //     log_in("central_atom", central_idx);
    //     log_in("idxs", frame_atoms);
    //     log_check("newcoords_SET1", newcoords);
    //     log_finish();
    // #endif
}

template <class T>
template <class C, class M>
void Polyhedron<T>::restore_from_frame(C& global_coords, const std::pair<int, int>& frame_atoms, const M* mol_data, const se_matrix_t& frame) {
    auto local_frame = get_frame(this->local_coords, frame_atoms, central_idx);
    // #ifdef KDMOL_LOG
    //     log_routine("polyrestore_pairA");
    //     log_in("central_atom", central_idx);
    //     log_in("frame_atoms", as_list(std::make_pair(frame_atoms.first, frame_atoms.second)));
    //     log_check("frame", frame);
    //     log_check("localframe", local_frame);
    //     log_finish();
    // #endif
    se_matrix_t transition_m, inv_local_frame;
    Utils::invert_matrix(local_frame, inv_local_frame);
    noalias(transition_m) = prod(frame, inv_local_frame);

    uvector4d_t coord_temp, local_temp;
    for (const auto& [atom_idx, atom_xyz] : local_coords) {
        if ((atom_idx == this->central_idx) || (atom_idx == frame_atoms.first))
            continue;
        for (int i = 0; i < 3; ++i)
            local_temp(i) = mol_data->bondlengths[mol_data->bond_data({this->central_idx, atom_idx}).index] * atom_xyz(i);
        local_temp(3) = 1.0;
        noalias(coord_temp) = prod(transition_m, local_temp);
        global_coords[atom_idx]; // Remove it!!
        for (int i = 0; i < 3; ++i)
            global_coords.at(atom_idx)(i) = coord_temp(i);
    }
}

template <class T>
template <class C, class M>
void Polyhedron<T>::restore_from_frame(C& global_coords, const int& frame_atom, const M* mol_data, const se_matrix_t& frame) {
    checkpoint("polyrestore_input");
    #ifdef KDMOL_LOG
        call_idxs["polyrestore_input"] += 1;
        log_routine("polyrestore_input");
        log_in("central_atom", central_idx);
        log_in("call_idx", call_idxs["polyrestore_input"]);
        log_check("frame_atom", as_list(frame_atom));
        log_check("frame", frame);
        log_finish();
    #endif

    auto local_frame = get_frame(this->local_coords, frame_atom, central_idx);
    
    checkpoint("polyrestore_singleA");
    #ifdef KDMOL_LOG
        call_idxs["polyrestore_singleA"] += 1;
        log_routine("polyrestore_singleA");
        log_in("central_atom", central_idx);
        log_in("frame_atom", frame_atom);
        log_in("call_idx", call_idxs["polyrestore_singleA"]);
        log_check("frame", frame);
        log_check("localframe", local_frame);
        log_finish();
    #endif
    se_matrix_t transition_m, inv_local_frame;
    Utils::invert_matrix(local_frame, inv_local_frame);
    noalias(transition_m) = prod(frame, inv_local_frame);

    uvector4d_t coord_temp, local_temp;
    for (const auto& [atom_idx, atom_xyz] : local_coords) {
        if ((atom_idx == this->central_idx) || (atom_idx == frame_atom)) {
            // fmt::print("----------> Skipped {} atom\n", atom_idx);
            continue;
        }
        for (int i = 0; i < 3; ++i)
            local_temp(i) = mol_data->bondlengths[mol_data->bond_data({this->central_idx, atom_idx}).index] * atom_xyz(i);
        local_temp(3) = 1.0;
        noalias(coord_temp) = prod(transition_m, local_temp);
        global_coords[atom_idx]; // Remove it!!!!
        // fmt::print("----------> Added {} atom\n", atom_idx);
        for (int i = 0; i < 3; ++i)
            global_coords.at(atom_idx)(i) = coord_temp(i);
    }
}

template <class T>
template < template<class ...> class PairContainer, class ... PairArgs,
           template<class ...> class ValueContainer, class ... ValueArgs >
void Polyhedron<T>::change_vangles(const PairContainer<std::pair<int, int>, PairArgs...>& vangle_items,
                    const ValueContainer<double, ValueArgs...>& vangle_values)
{
    checkpoint("changevangEnter");
    #ifdef KDMOL_LOG
        this->vangle_items = py::repr(py::handle(py::cast(vangle_items))).cast<std::string>();
        this->vangle_values = py::repr(py::handle(py::cast(vangle_values))).cast<std::string>();
        log_routine("changevangEnter");
        log_in("central_atom", central_idx);
        log_check("items", vangle_items);
        log_check("values", vangle_values);
        log_finish();
    #endif

    local_coords = coords_backup;

    const auto nx = Utils::PyModules::get_nx();
    auto gr = nx->attr("Graph")();
    gr.attr("add_edges_from")(vangle_items);
    auto conn_comps = py::list(nx->attr("connected_components")(gr));

    std::vector<py::object> subgraphs;
    subgraphs.reserve(conn_comps.size());
    #ifdef KDMOL_LOG
        py::list subgraph_nodes;
    #endif
    for(auto& comp : conn_comps) {
        #ifdef KDMOL_LOG
            auto comp_list = comp.cast<py::list>();
            comp_list.attr("sort")();
            subgraph_nodes.append(comp_list);
        #endif
        subgraphs.push_back(gr.attr("subgraph")(comp));
    }

    checkpoint("depgraphSubs");
    #ifdef KDMOL_LOG
        log_routine("depgraphSubs");
        log_in("central_atom", central_idx);
        log_check("items", vangle_items);
        log_check("subgraphs", subgraph_nodes);
        log_finish();
    #endif


    for (int i = 0; i < vangle_items.size(); ++i) {
        auto [ vA, vB ] = vangle_items[i];
        // It does "gr[item[0]][item[1]]['length'] = value"
        gr.attr("__getitem__")(vA).attr("__getitem__")(vB).attr("__setitem__")("length", vangle_values[i]);
        checkpoint("depgraphLength");
        #ifdef KDMOL_LOG
            log_routine("depgraphLength");
            log_in("central_atom", central_idx);
            log_in("vA", vA);
            log_in("vB", vB);
            log_check("length", as_list(gr.attr("__getitem__")(vA).attr("__getitem__")(vB).attr("__getitem__")("length")));
            log_finish();
        #endif
    }

    coordinates_t new_coords;
    // new_coords.set_atom(this->central_idx, {0.0, 0.0, 0.0});
    for (const auto& subgr : subgraphs) {
        auto n_nodes = subgr.attr("number_of_nodes")().cast<int>();
        auto n_edges = subgr.attr("number_of_edges")().cast<int>();

        if ((n_nodes == 3) && (n_edges == 3))
            resolve_triangle(subgr, new_coords);
        else if (n_nodes - n_edges == 1)
            resolve_branched(subgr, new_coords);
        else
            throw std::runtime_error(fmt::format("Constraint group for atoms {} is too complicated (n_nodes={}, n_edges={})", 
                        py::repr(py::handle(py::cast(vangle_items))).cast<std::string>(), n_nodes, n_edges));
        
        auto nodes = py::list(py::handle(subgr).attr("nodes")).cast<std::vector<int>>();
        find_best_fit(new_coords, nodes);
    }

    for (const auto& subgr : subgraphs) {
        auto nodes = py::list(py::handle(subgr).attr("nodes")).cast<std::vector<int>>();
        for (const auto& atom : nodes)
            local_coords[atom] = new_coords[atom];
    }

    checkpoint("changevangFinish");
    #ifdef KDMOL_LOG
        {
        py::list coord_list;
        for (auto [atom, coords] : local_coords) {
            if (atom != central_idx) {
                py::list new_item;
                new_item.append(py::cast(atom));
                new_item.append(py::cast(coords(0)));
                new_item.append(py::cast(coords(1)));
                new_item.append(py::cast(coords(2)));
                coord_list.append(new_item);
            }
        }
        log_routine("changevangFinish");
        log_in("central_atom", central_idx);
        log_in("items", vangle_items);
        log_in("values", vangle_values);
        log_check("coord_SET1", coord_list);
        log_finish();
        }
    #endif

    #ifdef VALIDATION
        // TODO Some double checking would be nice
    #endif
}

// TODO Random changes in directions for given subs. Use RMSD to fit directions (changed vs. original) 

template <class T>
template <class G, class C>
inline void Polyhedron<T>::resolve_triangle(G& subgr, C& new_coords) {
    const auto nodes = py::handle(subgr).attr("nodes").cast<py::list>();
    const auto a_len = py::handle(subgr).attr("__getitem__")(nodes[0].cast<int>()).attr("__getitem__")(nodes[1].cast<int>()).attr("__getitem__")("length").cast<double>();
    const auto b_len = py::handle(subgr).attr("__getitem__")(nodes[1].cast<int>()).attr("__getitem__")(nodes[2].cast<int>()).attr("__getitem__")("length").cast<double>();
    const auto c_len = py::handle(subgr).attr("__getitem__")(nodes[2].cast<int>()).attr("__getitem__")(nodes[0].cast<int>()).attr("__getitem__")("length").cast<double>();
    const double gamma = acos((cos(c_len) - cos(a_len)*cos(b_len)) / sin(a_len) / sin(b_len));

    checkpoint("polyTriangle");
    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();

        py::list set_lengths;
        if (nodes[0] < nodes[1])
            set_lengths.append(py::make_tuple(nodes[0], nodes[1], a_len));
        else
            set_lengths.append(py::make_tuple(nodes[1], nodes[0], a_len));
        if (nodes[1] < nodes[2])
            set_lengths.append(py::make_tuple(nodes[1], nodes[2], b_len));
        else
            set_lengths.append(py::make_tuple(nodes[2], nodes[1], b_len));
        if (nodes[2] < nodes[0])
            set_lengths.append(py::make_tuple(nodes[2], nodes[0], c_len));
        else
            set_lengths.append(py::make_tuple(nodes[0], nodes[2], c_len));

        log_routine("polyTriangle");
        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("lengths_SET2", set_lengths);
        log_check("gamma", as_list(gamma));
        log_finish();
        }
    #endif

    new_coords.set_atom(nodes[0].cast<int>(), {1.0, 0.0, 0.0});
    so_matrix_t main_frame, dmatA, dmatB, vmat; // OPT: less matrices/product calcs should be possible
    get_dmat_sph(a_len, dmatA);
    get_dmat_sph(b_len, dmatB);
    get_vmat_sph(gamma, vmat);
    noalias(main_frame) = prod(dmatA, vmat);
    new_coords.set_atom(nodes[1].cast<int>(), {main_frame(0, 0), main_frame(1, 0), main_frame(2, 0)});
    main_frame = prod(main_frame, dmatB);
    new_coords.set_atom(nodes[2].cast<int>(), {main_frame(0, 0), main_frame(1, 0), main_frame(2, 0)});

    checkpoint("polyTriangleFin");
    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();

        py::list new_coord_list;
        for (auto node : nodes) {
            py::list new_item;
            new_item.append(node);
            new_item.append(py::cast(new_coords[node.cast<int>()](0)));
            new_item.append(py::cast(new_coords[node.cast<int>()](1)));
            new_item.append(py::cast(new_coords[node.cast<int>()](2)));
            new_coord_list.append(new_item);
        }
        log_routine("polyTriangleFin");
        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("DmatA", dmatA);
        log_check("DmatB", dmatB);
        log_check("Vmat", vmat);
        log_check("coordA", new_coords[nodes[1].cast<int>()]);
        log_check("coordB", new_coords[nodes[2].cast<int>()]);
        log_check("coords_SET1", new_coord_list);
        log_finish();
        }
    #endif
}

template <class T>
template <class G, class C>
inline void Polyhedron<T>::resolve_branched(G& subgr, C& new_coords) {
    const auto nodes = py::list(py::handle(subgr).attr("nodes")).cast<std::vector<int>>();
    const auto headnode = nodes[0];
    const auto py_headnode = py::cast(headnode);

    std::unordered_map<int, int> idx_to_bgl;
    std::vector<int> endnodes;
    for (const auto node : nodes) {
        py::int_ py_node = py::cast(node);
        idx_to_bgl[node] = idx_to_bgl.size();
        if ((headnode != node) && (py::list(subgr.attr("neighbors")(py_node)).size() == 1))
            endnodes.push_back(node);
    }

    dgr_t dgr(nodes.size());
    for (const auto enode : endnodes) {
        auto solvepath = py::list(Utils::PyModules::get_nx()->attr("shortest_path")(subgr, py_headnode, py::cast(enode)));
        for (int i = 0; i < solvepath.size() - 1; ++i) {
            auto vA = vertex(idx_to_bgl[solvepath[i].cast<int>()], dgr),
                 vB = vertex(idx_to_bgl[solvepath[i + 1].cast<int>()], dgr);
            if (!edge(vA, vB, dgr).second)
                add_edge(vA, vB, dgr);
        }
    }

    checkpoint("polyBranched");
    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();
        py::list endnodes_sorted;
        for (auto& item : endnodes)
            endnodes_sorted.append(item);
        endnodes_sorted.attr("sort")();

        py::list dgr_nodes, dgr_edges;
        boost::graph_traits <dgr_t>::vertex_iterator i, end;
        boost::property_map <dgr_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dgr);

        for (boost::tie(i, end) = vertices(dgr); i != end; ++i) {
            dgr_nodes.append(nodes[get(index_map, *i)]);
        }
        
        boost::graph_traits<dgr_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(dgr); ei != ei_end; ++ei)
            dgr_edges.append(py::make_tuple(nodes[get(index_map, source(*ei, dgr))],
                                            nodes[get(index_map, target(*ei, dgr))]));

        log_routine("polyBranched");
        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("endnodes", endnodes_sorted);
        log_check("dgr_nodes_SET1", dgr_nodes);
        log_check("dgr_edges_SET2", dgr_edges);
        log_finish();
        }
    #endif


    std::vector<std::pair<int, int>> todolist;
    std::unordered_map<int, so_matrix_t> frames;
    todolist.reserve(nodes.size() - 1);
    new_coords[headnode] = local_coords[headnode];
    boost::graph_traits<dgr_t>::out_edge_iterator out, out_end;
    boost::property_map <dgr_t, boost::vertex_index_t>::type dgr_index_map = get(boost::vertex_index, dgr);
    for (boost::tie(out, out_end) = out_edges(vertex(idx_to_bgl[headnode], dgr), dgr); out != out_end; ++out) {
        const auto nb = nodes[get(dgr_index_map, target(*out, dgr))];

        so_matrix_t dmat;
        const auto shift_length = py::handle(subgr).attr("__getitem__")(headnode).attr("__getitem__")(nb).attr("__getitem__")("length").cast<double>();
        get_dmat_sph(shift_length, dmat);
        auto startframe = get_frame_3d(local_coords, {headnode, nb});
        auto& frame = frames[nb];
        noalias(frame) = prod(startframe, dmat);
        new_coords.set_atom(nb, {frame(0, 0), frame(1, 0), frame(2, 0)});
        todolist.push_back(std::make_pair(headnode, nb));
    
        checkpoint("polyTodoCheck");
        #ifdef KDMOL_LOG
            {
            py::list nodes_sorted;
            for (auto& item : nodes)
                nodes_sorted.append(item);
            nodes_sorted.attr("sort")();

            log_routine("polyTodoCheck");
            log_in("central_atom", central_idx);
            log_in("nodes", nodes_sorted);
            log_in("nb", nb);
            log_in("items", this->vangle_items);
            log_in("values", this->vangle_values);
            log_check("todolist_SET2", todolist);
            log_check("start_frame", startframe);
            log_check("dmat", dmat);
            log_check("frame", frame);
            log_check("coords", new_coords[nb]);
            log_finish();
            }
        #endif
    }

    while(todolist.size() > 0) {
        auto [prevatom, curatom] = todolist.back();
        todolist.pop_back();
        for (boost::tie(out, out_end) = out_edges(vertex(idx_to_bgl[curatom], dgr), dgr); out != out_end; ++out) {
            const auto nb = nodes[get(dgr_index_map, target(*out, dgr))];
            const auto ang = get_sph_angle(prevatom, curatom, nb);
            const auto length = py::handle(subgr).attr("__getitem__")(curatom).attr("__getitem__")(nb).attr("__getitem__")("length").cast<double>();

            so_matrix_t dmat, vmat;
            get_dmat_sph(length, dmat);
            get_vmat_sph(ang, vmat);

            auto& new_frame = frames[nb];
            noalias(new_frame) = prod(frames[curatom], vmat);
            new_frame = prod(new_frame, dmat);
            new_coords.set_atom(nb, {new_frame(0, 0), new_frame(1, 0), new_frame(2, 0)});
            todolist.push_back(std::make_pair(curatom, nb));

            checkpoint("polyBranchedWhile");
            #ifdef KDMOL_LOG
                {
                py::list nodes_sorted;
                for (auto& item : nodes)
                    nodes_sorted.append(item);
                nodes_sorted.attr("sort")();

                log_routine("polyBranchedWhile");
                log_in("central_atom", central_idx);
                log_in("nodes", nodes_sorted);
                log_in("nb", nb);
                log_in("items", this->vangle_items);
                log_in("values", this->vangle_values);
                log_check("ang", as_list(ang));
                log_check("length", as_list(length));
                log_check("vmat", vmat);
                log_check("dmat", dmat);
                log_check("frame", new_frame);
                log_check("coords", new_coords[nb]);
                log_finish();
                }
            #endif
        }
    }
    checkpoint("polyBranchedFinal");
    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();

        py::list coords;
        for (auto node : nodes) {
            py::list new_item;
            new_item.append(py::cast(node));
            new_item.append(py::cast(new_coords[node](0)));
            new_item.append(py::cast(new_coords[node](1)));
            new_item.append(py::cast(new_coords[node](2)));
            coords.append(new_item);
        }

        log_routine("polyBranchedFinal");
        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("coords_SET1", coords);
        log_finish();
        }
    #endif
}

template <class T>
template <class C, template<class ...> class NodeContainer, class ... NodeArgs>
inline void Polyhedron<T>::find_best_fit(C& new_coords, const NodeContainer<int, NodeArgs...>& nodes) {
    umatrix_t new_coord_m(nodes.size(), 3);
    umatrix_t ref_coord_m(nodes.size(), 3);
    for (int i = 0; i < nodes.size(); ++i) {
        for (int j = 0; j < 3; ++j)
            new_coord_m(i, j) = new_coords.at(nodes[i])(j);
        for (int j = 0; j < 3; ++j)
            ref_coord_m(i, j) = coords_backup.at(nodes[i])(j);
    }

    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();

        py::list new_coord_list, ref_coord_list;
        for (int i = 0; i < nodes.size(); ++i) {
            py::list new_item, ref_item;
            new_item.append(py::cast(nodes[i]));
            new_item.append(py::cast(new_coord_m(i, 0)));
            new_item.append(py::cast(new_coord_m(i, 1)));
            new_item.append(py::cast(new_coord_m(i, 2)));
            new_coord_list.append(new_item);
            
            ref_item.append(py::cast(nodes[i]));
            ref_item.append(py::cast(ref_coord_m(i, 0)));
            ref_item.append(py::cast(ref_coord_m(i, 1)));
            ref_item.append(py::cast(ref_coord_m(i, 2)));
            ref_coord_list.append(ref_item);
        }

        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("new_coord_SET1", new_coord_list);
        log_check("ref_coord_SET1", ref_coord_list);
        log_finish();
        }
    #endif
    
    Utils::umatrix3d_t C_m, V, W, U;
    noalias(C_m) = prod(trans(new_coord_m), ref_coord_m);
    
    Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> C_ei(&C_m(0, 0));
    Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(C_ei);
    
    auto U_ei = svd.matrixU();
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            V(i, j) = U_ei(i, j);
    auto V_ei = svd.matrixV();
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            W(i, j) = V_ei(j, i);
    noalias(U) = prod(V, W);
    new_coord_m = prod(new_coord_m, U); // OPT: use noalias
    
    for (int i = 0; i < nodes.size(); ++i)
        new_coords.set_atom(nodes[i], {new_coord_m(i, 0), new_coord_m(i, 1), new_coord_m(i, 2)});
    
    checkpoint("polyRmsdFinal");
    #ifdef KDMOL_LOG
        {
        py::list nodes_sorted;
        for (auto& item : nodes)
            nodes_sorted.append(item);
        nodes_sorted.attr("sort")();

        py::list coord_list;
        for (auto node : nodes) {
            py::list new_item;
            new_item.append(py::cast(node));
            new_item.append(py::cast(new_coords[node](0)));
            new_item.append(py::cast(new_coords[node](1)));
            new_item.append(py::cast(new_coords[node](2)));
            coord_list.append(new_item);
        }

        log_routine("polyRmsdFinal");
        log_in("central_atom", central_idx);
        log_in("nodes", nodes_sorted);
        log_in("items", this->vangle_items);
        log_in("values", this->vangle_values);
        log_check("new_coord_list_SET1", coord_list);
        log_finish();
        }
    #endif
}

#endif // POLY_H
