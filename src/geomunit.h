#ifndef GEOMUNIT_H
#define GEOMUNIT_H

#include "utils.h"

template <class MD, class DF, class DD, class AD>
class GeomUnit {
    public:
        using mol_data_t = MD;
        using dofs_t = DF;
        using discrete_dofs_t = DD;
        using assembly_data_t = AD;
        using coord_container_t = Utils::MapXyzContainer<std::map<int, Utils::uvector3d_t>>;
        using inner_container_t = typename assembly_data_t::inner_container_t;
        using problem_idxs_t = std::vector<int>;
        using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;
        using protected_xyzs_t = std::vector<Utils::uvector3d_t>;

        #ifdef OVERLAP_DETECTION
            using overlap_detector_t = Utils::OverlapDetector;
        #endif

        GeomUnit () : 
        #ifdef KDMOL_LOG
            idx(-1), ord(-1), 
        #endif
        headnode(-1) { }
        
        GeomUnit (py::object& subG, py::object& subLG, const int given_headnode, mol_data_t* mol_data);
        template <class C> aps_return_t apply_ps(C& coords, mol_data_t* mol_data, dofs_t* dofs);

        bool contains_problem(const int prob_idx) const noexcept;
        void index_owned_dofs(const active_idxs_t& own_dof_idxs, const active_idxs_t& own_discrete_idxs);

        template <class C> void protect_atoms(const C& conflicting_atoms);
        void finalize_atoms_protection();
        template <class C> inline void restore_protected(C& coords);

        #ifdef VALIDATION
            using validator_t = Utils::LocalGeometryValidator;
            validator_t conf_validator, ext_validator;
            void build_validator(mol_data_t* mol_data, dofs_t* dofs);
        #endif

        #ifdef OVERLAP_DETECTION
            overlap_detector_t overlap_detector;
        #endif

        #ifdef KDMOL_LOG
            int idx, ord;
        #endif
        int headnode, connect_node;
        problem_idxs_t problem_idxs;
        inner_container_t inner_container;
        active_idxs_t own_dof_idxs, own_discrete_idxs;
        
        std::vector<int> protected_atoms;
        protected_xyzs_t protected_xyzs;
        bool atom_protection_needed;
};


template <class MD, class DF, class DD, class AD>
GeomUnit<MD, DF, DD, AD>::GeomUnit (py::object& subG, py::object& subLG, const int given_headnode, mol_data_t* mol_data) :
    #ifdef KDMOL_LOG
        idx(-1), ord(-1), 
    #endif
    headnode(given_headnode), connect_node(-1), atom_protection_needed(false)
{
    for (const auto& py_node : subG.attr("nodes")()) {
        const auto node = py::handle(py_node).cast<py::int_>();
        auto nb_list = py::list(subG.attr("neighbors")(node));
        nb_list.attr("sort")();
        inner_container[node] = std::make_pair(nb_list[0].cast<int>(),
                                               nb_list[1].cast<int>());
    }

    for (const auto& py_idx : subLG.attr("nodes")()) {
        auto idx = py::handle(py_idx).cast<int>();
        problem_idxs.push_back(idx);
    }
    problem_idxs.shrink_to_fit();

    checkpoint("gu_innerframe");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms, log_frames;
        for (const auto& py_idx : subG.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();

        for (const auto& [ center, frame ] : inner_container) {
            py::list newitem;
            newitem.append(center);
            newitem.append(frame.first);
            newitem.append(frame.second);
            log_frames.append(newitem);
        }
        
        log_routine("gu_innerframe");
        log_in("atoms", log_atoms);
        log_check("innerframes_SET1", log_frames);
        log_finish();
        }
    #endif

    #ifdef OVERLAP_DETECTION
        {
        std::set<int> ext_atoms;
        for (const auto& [ atom, _ ] : inner_container) {
            ext_atoms.insert(atom);
            for (const auto& [ nbatom, __ ] : mol_data->polys[atom])
                ext_atoms.insert(nbatom);
        }
        assert(ext_atoms.size() != 0);
        overlap_detector = overlap_detector_t(mol_data->molgr, ext_atoms, mol_data->atom_symbols, "geomunit");
        }
    #endif
}

template <class MD, class DF, class DD, class AD>
bool GeomUnit<MD, DF, DD, AD>::contains_problem(const int prob_idx) const noexcept {
    auto found_iter = std::find(problem_idxs.cbegin(), problem_idxs.cend(), prob_idx);
    return found_iter != problem_idxs.cend();
}

#ifdef VALIDATION
template <class MD, class DF, class DD, class AD>
void GeomUnit<MD, DF, DD, AD>::build_validator(mol_data_t* mol_data, dofs_t* dofs) {
    std::set<int> atoms;
    for (const auto& [ atom, _ ] : inner_container)
        atoms.insert(atom);
    conf_validator = validator_t(*mol_data, *dofs, atoms);

    std::set<int> ext_atoms;
    for (const auto& [ atom, _ ] : inner_container) {
        ext_atoms.insert(atom);
        for (const auto& [ nbatom, __ ] : mol_data->polys[atom])
            ext_atoms.insert(nbatom);
    }
    ext_validator = validator_t(*mol_data, *dofs, ext_atoms);
}
#endif

template <class MD, class DF, class DD, class AD>
void GeomUnit<MD, DF, DD, AD>::index_owned_dofs(const active_idxs_t& own_dof_idxs, const active_idxs_t& own_discrete_idxs) {
    this->own_dof_idxs = own_dof_idxs;
    this->own_discrete_idxs = own_discrete_idxs;
}

template <class MD, class DF, class DD, class AD>
template <class C>
void GeomUnit<MD, DF, DD, AD>::protect_atoms(const C& given_atoms) {
    for (const auto atom : given_atoms) {
        auto it = inner_container.find(atom);
        if (it == inner_container.cend())
            protected_atoms.push_back(atom);
    }
}

template <class MD, class DF, class DD, class AD>
void GeomUnit<MD, DF, DD, AD>::finalize_atoms_protection() {
    protected_atoms.shrink_to_fit();
    protected_xyzs = protected_xyzs_t(protected_atoms.size());
    atom_protection_needed = (protected_atoms.size() > 0);
}

template <class MD, class DF, class DD, class AD>
template <class C>
inline void GeomUnit<MD, DF, DD, AD>::restore_protected(C& coords) {
    if (atom_protection_needed)
        for (int i = 0; i < protected_atoms.size(); ++i)
            coords[protected_atoms[i]] = protected_xyzs[i];
}

template <class MD, class DF, class DD, class AD>
template <class C>
aps_return_t GeomUnit<MD, DF, DD, AD>::apply_ps(C& coords, mol_data_t* mol_data, dofs_t* dofs) {
    if (atom_protection_needed)
        for (int i = 0; i < protected_atoms.size(); ++i)
            protected_xyzs[i] = coords[protected_atoms[i]];

    for (const auto [ atom, nbs ] : inner_container)
        mol_data->polys[atom].restore_neighbors(coords, nbs, mol_data);

    #ifdef VALIDATION
        if (!conf_validator.validate(coords, *mol_data, *dofs)) {
            mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
            restore_protected(coords);
            return aps_return_t::validationfail;
            // throw Utils::GeomFailure("Validation of geomunit conformation failed");
        }
        if (!ext_validator.validate(coords, *mol_data, *dofs)) {
            mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
            restore_protected(coords);
            return aps_return_t::validationfail;
            // throw Utils::GeomFailure("Validation of geomunit's extended coords failed");
        }
    #endif

    #ifdef OVERLAP_DETECTION
        if (!overlap_detector(coords)) {
            mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
            restore_protected(coords);
            return aps_return_t::overlap;
            // throw Utils::OverlapDetected(OVERLAP_GEOMUNIT_MESSAGE);
        }
    #endif
    restore_protected(coords);
    return aps_return_t::success;
}

#endif // GEOMUNIT_H
