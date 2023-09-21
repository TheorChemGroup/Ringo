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

#endif // POLY_H
