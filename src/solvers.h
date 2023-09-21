#ifndef SOLVERS_H
#define SOLVERS_H

#include "pytlc/ringext.h"


template <class MD, class CD, class C>
class AbstractSolver {
    using mol_data_t = MD;
    using cyclic_dofs_t = CD;
    using extcoord_t = C;
    using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;

    public:
        using coord_t = Utils::MapXyzContainer<std::map<int, Utils::uvector3d_t>>;
        using dofs_container_t = std::vector<Utils::fulldof_t>;
        using marked_dofs_t = std::unordered_map<std::string, std::vector<Utils::fulldof_t>>;
        using requested_free_t = std::vector<Utils::dof_t>;
        
        extcoord_t* extcoord;

        virtual dofs_container_t mydofs(py::object& gr, py::list& py_consbonds, const requested_free_t& requested_free) = 0;
        virtual void finalize_solver(const int firstdof_idx, mol_data_t* mol_data,  cyclic_dofs_t* cyclic_dofs) = 0;
        virtual aps_return_t apply_ps(mol_data_t* mol_data) = 0;
        virtual extcoord_t get_coords() = 0;
        virtual const std::set<int> get_atoms_idxs() = 0;

        // Owned DOFs for refinement
        active_idxs_t own_dof_idxs, own_discrete_idxs;
        void index_owned_dofs(const active_idxs_t& own_dof_idxs, const active_idxs_t& own_discrete_idxs) {
            this->own_dof_idxs = own_dof_idxs;
            this->own_discrete_idxs = own_discrete_idxs;
        }
};


template <class MD, class CD, class C>
class IdentitySolver : public AbstractSolver<MD, CD, C> {
        using parent_t = AbstractSolver<MD, CD, C>;
    public:
        using coord_t = typename parent_t::coord_t;
    private:
        using mol_data_t = MD;
        using cyclic_dofs_t = CD;
        using extcoord_t = C;
        using fulldof_t = Utils::fulldof_t;
        using dofs_container_t = typename parent_t::dofs_container_t;
        using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;
        using requested_free_t = typename parent_t::requested_free_t;

        coord_t conformer;
        unsigned int natoms;
        dofs_container_t check_dihedrals; // This is needed only for testing
        double* dofvalue_pointer;
        int firstdof_idx;
        #ifdef KDMOL_LOG
            py::list py_atoms;
        #endif
    
    public:
        IdentitySolver() { }
        
        dofs_container_t mydofs(py::object& gr, py::list& py_consbonds, const requested_free_t& requested_free) override;
        void finalize_solver(const int firstdof_idx, mol_data_t* mol_data,  cyclic_dofs_t* cyclic_dofs) override;
        aps_return_t apply_ps(mol_data_t* mol_data) override;
        const std::set<int> get_atoms_idxs() override;
        coord_t get_coords() override;
};


class TLCSolverModel {
    public:
        TLCSolverModel() { }
        template <class B, class C> static std::tuple<bool, std::string, std::vector<Utils::dof_t>> can_apply(py::object& gr, B& consbonds, const C& requested_free);
    
    protected:
        template <class B, class M> static std::vector<unsigned int> get_constraints(py::object& gr, B& consbonds, const M& idxmap);
        template <class B, class M> static std::vector<unsigned int> format_requested_free(B& requested_free, const M& idxmap);
        static std::pair<std::vector<int>, std::unordered_map<int, int>> walk_cycle(py::object& gr);
};


template <class MD, class CD, class C>
class TLCSolver : public AbstractSolver<MD, CD, C>, public TLCSolverModel {
        using parent_t = AbstractSolver<MD, CD, C>;
    public:
        using coord_t = typename parent_t::coord_t;
    private:
        using mol_data_t = MD;
        using cyclic_dofs_t = CD;
        using extcoord_t = C;
        using fulldof_t = Utils::fulldof_t;
        using dofs_container_t = typename parent_t::dofs_container_t;
        using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;
        using requested_free_t = typename parent_t::requested_free_t;

        TLCExtSolver::axial_atoms_t axial_atoms; // OPT: Used only for postponed initialization
        TLCExtSolver::fulldof_data_t fulldof_data; // OPT: Used only for postponed initialization
        TLCExtSolver::fulldof_vector_t fulldof_list; // OPT: Used only for postponed initialization

        unsigned int nfixed, nfree, natoms;
        TLCExtSolver driver;
        std::vector<int> idx2atom;
        std::unordered_map<int, int> atom2idx;
        double* dofvalue_pointer;
        int* discrete_dofvalue_pointer;
        int firstdof_idx;
        TLCExtSolver::fulldof_data_t fulldofs_final;

    public:
        TLCSolver() { }
        coord_t get_coords() override;
        const auto& get_atoms() const { return idx2atom; }

        dofs_container_t mydofs(py::object& gr, py::list& py_consbonds, const requested_free_t& requested_free) override;
        void finalize_solver(const int firstdof_idx, mol_data_t* mol_data,  cyclic_dofs_t* cyclic_dofs) override;
        aps_return_t apply_ps(mol_data_t* mol_data) override;
        const std::set<int> get_atoms_idxs() override;
        void set_ddof_pointer(const int prob_idx, cyclic_dofs_t* cyclic_dofs);
};

template <class B, class M>
std::vector<unsigned int> TLCSolverModel::format_requested_free(B& requested_free, const M& idxmap) {
    std::set<unsigned int> idxset;
    for (const auto& bond : requested_free) {
        const int atA = idxmap.at(bond.first);
        const int atB = idxmap.at(bond.second);
        
        int bondidx;
        if (abs(atA - atB) == 1) {
            if (atA < atB)
                bondidx = atA;
            else
                bondidx = atB;
        } else { // Someone is equal to 0
            if (atA < atB)
                bondidx = atB;
            else
                bondidx = atA;
        }
        idxset.insert(bondidx);
    }
    std::vector<unsigned int> bondidxs(idxset.cbegin(), idxset.cend());
    return bondidxs;
}

template <class B, class M>
std::vector<unsigned int> TLCSolverModel::get_constraints(py::object& gr, B& consbonds, const M& idxmap) {
    std::set<unsigned int> idxset;
    for (const auto& item : consbonds) {
        std::vector<int> bond = py::handle(item).cast<std::vector<int>>();
        const int atA = idxmap.at(bond[0]);
        const int atB = idxmap.at(bond[1]);
        
        int bondidx;
        if (abs(atA - atB) == 1) {
            if (atA < atB)
                bondidx = atA;
            else
                bondidx = atB;
        } else { // Someone is equal to 0
            if (atA < atB)
                bondidx = atB;
            else
                bondidx = atA;
        }
        idxset.insert(bondidx);
    }

    for (const auto& py_edge : gr.attr("edges")()) {
        auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
        if (!gr.attr("__getitem__")(edge.first).attr("__getitem__")(edge.second).attr("__getitem__")("type").cast<bool>()) {
            int bondidx;
            const int atA = idxmap.at(edge.first);
            const int atB = idxmap.at(edge.second);
            if (abs(atA - atB) == 1) {
                if (atA < atB)
                    bondidx = atA;
                else
                    bondidx = atB;
            } else { // Someone is equal to 0
                if (atA < atB)
                    bondidx = atB;
                else
                    bondidx = atA;
            }
            // bondidxs.push_back(bondidx);
            idxset.insert(bondidx);
        }
    }
    
    std::vector<unsigned int> bondidxs(idxset.cbegin(), idxset.cend());
    return bondidxs;
}


template <class B, class C>
std::tuple<bool, std::string, std::vector<Utils::dof_t>> TLCSolverModel::can_apply(py::object& gr, B& consbonds, const C& requested_free) {
    const int nnodes = gr.attr("number_of_nodes")().cast<int>();
    const int nedges = gr.attr("number_of_edges")().cast<int>();
    const int ncycles = nedges - nnodes + 1;
    std::vector<Utils::dof_t> uncomplied_requests;
    if (ncycles <= 0)
        std::runtime_error(fmt::format("The provided cyclic part appears to have no rings (nedges - nnodes = {})", nedges - nnodes));
    else if (ncycles > 1) {
        const auto node_list = gr.attr("nodes").cast<py::list>().cast<std::vector<int>>();
        const auto warning_message = prepare_message(Utils::warning_codes["IK_NOT_APPLIED"], fmt::format("IK is unapplicable to rigid polycyclic (sub)fragment of {} atoms. This fragment will be kept as it is in the starting conformation.", nnodes), node_list);
        return std::make_tuple(false, warning_message, uncomplied_requests); // Return non-important message
    }

    const auto [ idx2atom, atom2idx ] = walk_cycle(gr);
    const auto bondidxs = get_constraints(gr, consbonds, atom2idx);

    // Very simple check if Nrotatable >= 6
    if (static_cast<int>(nedges) - static_cast<int>(bondidxs.size()) < 6) {
        const auto warning_message = prepare_message(Utils::warning_codes["IK_NOT_APPLIED"], fmt::format("IK is unapplicable to rigid {}-membered ring. This cycle will be kept as it is in the starting conformation.", nnodes), idx2atom);
        return std::make_tuple(false, warning_message, uncomplied_requests);
    }

    const auto requested_free_formatted = format_requested_free(requested_free, atom2idx);
    auto [ res, uncomplied_requests_raw ] = TLCExtSolver::check_applicable(nnodes, bondidxs, requested_free_formatted);

    if (res) {
        if (uncomplied_requests_raw.size() > 0) {
            for (const auto& bond : uncomplied_requests_raw)
                uncomplied_requests.push_back(std::make_pair(idx2atom[bond.first],
                                                             idx2atom[bond.second]));
        }
        return std::make_tuple(true, "", uncomplied_requests);
    } else {
        const auto warning_message = prepare_message(Utils::warning_codes["IK_NOT_APPLIED"]+"[important]", fmt::format("TLC can not be applied to flexible {}-membered cycle due to coterminal pair restriction. This cycle will be kept as it is in the starting conformation.", nnodes), idx2atom);
        return std::make_tuple(false, warning_message, uncomplied_requests);
    }
}

template <class MD, class CD, class C>
typename TLCSolver<MD, CD, C>::dofs_container_t TLCSolver<MD, CD, C>::mydofs(py::object& gr, py::list& py_consbonds, const requested_free_t& requested_free) {
    this->natoms = gr.attr("number_of_nodes")().cast<int>();
    std::tie(idx2atom, atom2idx) = this->walk_cycle(gr);
    const auto bondidxs = this->get_constraints(gr, py_consbonds, atom2idx);
    const auto requested_free_formatted = this->format_requested_free(requested_free, atom2idx);

    std::tie(fulldof_data, axial_atoms) = TLCExtSolver::propose_dofs(natoms, bondidxs, idx2atom, requested_free_formatted);

    checkpoint("tlc_init");
    #ifdef KDMOL_LOG
        py::list py_atoms_idx, py_axatoms;
        for (const auto i : axial_atoms)
            py_axatoms.append(py::cast(idx2atom.at(i)));
        py_axatoms.attr("sort")();

        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            py_atoms_idx.append(idx);
        }
        py_atoms_idx.attr("sort")();

        log_routine("tlc_init");
        log_in("atoms", py_atoms_idx);
        log_check("axatoms", py_axatoms);
        log_check("consbonds_SET2", py_consbonds);
        log_check("walk", idx2atom);
        log_finish();
    #endif

    
    fulldof_list.reserve(fulldof_data.at("variable").size() + fulldof_data.at("fixed").size());
    for (const auto& item : fulldof_data.at("variable")) {
        fulldof_t newitem = {idx2atom[std::get<0>(item)], idx2atom[std::get<1>(item)], idx2atom[std::get<2>(item)], idx2atom[std::get<3>(item)]};
        fulldof_list.push_back(newitem);
    }
    for (const auto& item : fulldof_data.at("fixed")) {
        fulldof_t newitem = {idx2atom[std::get<0>(item)], idx2atom[std::get<1>(item)], idx2atom[std::get<2>(item)], idx2atom[std::get<3>(item)]};
        fulldof_list.push_back(newitem);
    }
    return fulldof_list; // List of all dihedrals that should be provided (fixed + resolved). THE ORDER IS PRELIMINARY!
}

template <class MD, class CD, class C>
void TLCSolver<MD, CD, C>::finalize_solver (const int prob_idx,
                                            mol_data_t* mol_data,
                                            cyclic_dofs_t* cyclic_dofs)
{
    firstdof_idx = cyclic_dofs->fragment_starts[prob_idx];
    dofvalue_pointer = &cyclic_dofs->dofvalue_container[firstdof_idx];

    fulldofs_final = {
                        {"variable", dofs_container_t()},
                        {"solved", dofs_container_t()},
                        {"fixed", dofs_container_t()},
                    };
    for (const auto& dof : fulldof_data["solved"])
        fulldofs_final["solved"].push_back(dof);
    for (const auto& dof : fulldof_data["variable"])
        fulldofs_final["variable"].push_back(dof);
    for (const auto& dof : fulldof_data["fixed"]) {
        fulldof_t dof_gidx = {idx2atom[std::get<0>(dof)], idx2atom[std::get<1>(dof)],
                              idx2atom[std::get<2>(dof)], idx2atom[std::get<3>(dof)]};
        const auto dof_idx = cyclic_dofs->dofindex_in_frag(dof_gidx, prob_idx);
        auto doftype = cyclic_dofs->doftype_container[dof_idx];
        if (doftype == cyclic_dofs_t::dof_type_t::dep) {
            fulldofs_final["variable"].push_back(dof);
        } else if (doftype == cyclic_dofs_t::dof_type_t::fixed) {
            fulldofs_final["fixed"].push_back(dof);
        }
    }
    
    auto dof_data = TLCExtSolver::condense_dofs(fulldofs_final); // TODO CHECK THE ORDER OF VARIABLE AND FIXED DOFS!!!!!!!!!

    nfree = dof_data.at("variable").size();
    nfixed = dof_data.at("fixed").size();

    std::vector<double> length_data, vangle_data, fixed_dihedral_data;
    for (int i = 1; i < natoms; ++i)
        length_data.push_back(mol_data->bondlengths.at(mol_data->bond_data({idx2atom[i - 1], idx2atom[i]}).index));
    length_data.push_back(mol_data->bondlengths.at(mol_data->bond_data({idx2atom[natoms - 1], idx2atom[0]}).index));

    vangle_data.push_back(mol_data->polys.at(idx2atom[0]).get_vangle_rad(idx2atom[natoms - 1], idx2atom[1]));
    for (int i = 2; i < natoms; ++i)
        vangle_data.push_back(mol_data->polys.at(idx2atom[i - 1]).get_vangle_rad(idx2atom[i - 2], idx2atom[i]));
    vangle_data.push_back(mol_data->polys.at(idx2atom[natoms - 1]).get_vangle_rad(idx2atom[natoms - 2], idx2atom[0]));
    
    for (int i = 0; i < nfixed; ++i)
        fixed_dihedral_data.push_back(*(dofvalue_pointer + nfree + i));

    driver = TLCExtSolver(natoms, dof_data, this->axial_atoms,
                          length_data, vangle_data, fixed_dihedral_data);
    
    checkpoint("tlc_dofcheck");
    #ifdef KDMOL_LOG
        py::list py_vardata, py_fixeddata;
        for (int i = 0; i < fulldofs_final["variable"].size(); ++i) {
            py_vardata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<1>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<2>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<3>(fulldofs_final["variable"].at(i))]));
        }
        py_vardata.attr("sort");
        for (int i = 0; i < fulldofs_final["fixed"].size(); ++i) {
            py_fixeddata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["fixed"].at(i))],
                                             idx2atom[std::get<1>(fulldofs_final["fixed"].at(i))],
                                             idx2atom[std::get<2>(fulldofs_final["fixed"].at(i))],
                                             idx2atom[std::get<3>(fulldofs_final["fixed"].at(i))]));
        }
        py_fixeddata.attr("sort");
        
        py::list py_atoms = py::cast(idx2atom);
        py_atoms.attr("sort")();

        log_routine("tlc_dofcheck");
        log_in("atoms", py_atoms);
        log_check("variable_SET4", py_vardata);
        log_check("fixed_SET4", py_fixeddata);
        log_finish();
    #endif
}

template <class MD, class CD, class C>
void TLCSolver<MD, CD, C>::set_ddof_pointer(const int prob_idx, cyclic_dofs_t* cyclic_dofs) {
    discrete_dofvalue_pointer = &cyclic_dofs->molvalues_discrete_start[cyclic_dofs->cp2mol_discrete.at(prob_idx)];
    checkpoint("tlc_discrete_dof");
    #ifdef KDMOL_LOG
        {
        py::list py_atoms = py::cast(idx2atom);
        py_atoms.attr("sort")();
        
        log_routine("tlc_discrete_dof");
        log_in("atoms", py_atoms);
        log_in("atoms", py_atoms);
        log_finish();
        }
    #endif
}

template <class MD, class CD, class C>
const std::set<int> TLCSolver<MD, CD, C>::get_atoms_idxs() {
    return std::set<int>(idx2atom.cbegin(), idx2atom.cend());
}

template <class MD, class CD, class C>
aps_return_t TLCSolver<MD, CD, C>::apply_ps(mol_data_t* mol_data) {
    std::vector<double> free_dofs(nfree);
    std::copy(dofvalue_pointer, dofvalue_pointer + nfree, free_dofs.begin());
    auto& discrete_dof = *discrete_dofvalue_pointer;

    checkpoint("tlc_enter");
    #ifdef KDMOL_LOG
        {
        py::list py_dofdata;
        for (int i = 0; i < fulldofs_final["variable"].size(); ++i) {
            py_dofdata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<1>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<2>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<3>(fulldofs_final["variable"].at(i))],
                                             free_dofs[i]));
        }
        
        py::list py_atoms = py::cast(idx2atom);
        py_atoms.attr("sort")();


        log_routine("tlc_enter");
        log_in("atoms", py_atoms);
        log_in("call_idx", Logger::aps_count);
        log_check("param_data_SET4", py_dofdata);
        log_check("ddof_value", as_list(discrete_dof));
        log_finish();
        }
    #endif

    auto nsol = driver.solve(free_dofs);
    checkpoint(fmt::format("nsol = {}", nsol));
    if ((discrete_dof == -1) && (nsol != 0))
        discrete_dof = mol_data->refiners.randomize_ddof(nsol, *mol_data);
    else
        assertm(discrete_dof < nsol, "Not enough solutions");
    
    if (nsol == 0) {
        if (driver.are_some_removed()) {
            #ifdef KDMOL_LOG
                {
                py::list py_dofdata, py_coord;
                for (int i = 0; i < fulldofs_final["variable"].size(); ++i) {
                    py_dofdata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<1>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<2>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<3>(fulldofs_final["variable"].at(i))],
                                                    free_dofs[i]));
                }
                
                py::list py_atoms = py::cast(idx2atom);
                py_atoms.attr("sort")();

                // log_routine("tlc_checksol");
                // log_in("atoms", py_atoms);
                // log_in("call_idx", Logger::aps_count);
                // log_check("param_data_SET4", py_dofdata);
                // log_check("outcome", as_list("fail"));
                // log_finish();
                }
            #endif
            mol_data->refiners.set_active(this->own_dof_idxs, this->own_discrete_idxs);
            return aps_return_t::tlcfail;
            // throw Utils::TLCFailure("TLC failed");
        } else {
            #ifdef KDMOL_LOG
                {
                py::list py_dofdata, py_coord;
                for (int i = 0; i < fulldofs_final["variable"].size(); ++i) {
                    py_dofdata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<1>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<2>(fulldofs_final["variable"].at(i))],
                                                    idx2atom[std::get<3>(fulldofs_final["variable"].at(i))],
                                                    free_dofs[i]));
                }
                
                py::list py_atoms = py::cast(idx2atom);
                py_atoms.attr("sort")();

                // log_routine("tlc_checksol");
                // log_in("atoms", py_atoms);
                // log_in("call_idx", Logger::aps_count);
                // log_check("param_data_SET4", py_dofdata);
                // log_check("outcome", as_list("zero"));
                // log_finish();
                }
            #endif
            mol_data->refiners.set_active(this->own_dof_idxs, this->own_discrete_idxs);
            return aps_return_t::zero;
            // throw Utils::TLCZeroSolutions("TLC has no solutions");
        }
    }

    checkpoint("tlc_checksol");
    #ifdef KDMOL_LOG
        {
        py::list py_dofdata, py_coord;
        for (int i = 0; i < fulldofs_final["variable"].size(); ++i) {
            py_dofdata.append(py::make_tuple(idx2atom[std::get<0>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<1>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<2>(fulldofs_final["variable"].at(i))],
                                             idx2atom[std::get<3>(fulldofs_final["variable"].at(i))],
                                             free_dofs[i]));
        }

        const int sol_idx = discrete_dof;
        boost::numeric::ublas::matrix<double> coords_raw(natoms, 3);
        driver.solution_to_buffer(&coords_raw(0, 0), sol_idx);

        for (int i = 0; i < natoms; ++i) {
            py::list newitem;
            newitem.append(idx2atom[i]);
            newitem.append(coords_raw(i, 0));
            newitem.append(coords_raw(i, 1));
            newitem.append(coords_raw(i, 2));
            py_coord.append(newitem);
        }
        
        py::list py_atoms = py::cast(idx2atom);
        py_atoms.attr("sort")();

        // nsol does not match with reference implementation for cycle_2

        log_routine("tlc_checksol");
        log_in("atoms", py_atoms);
        log_in("call_idx", Logger::aps_count);
        log_check("param_data_SET4", py_dofdata);
        log_check("ddof_value", as_list(discrete_dof));
        // log_check("nsol", as_list(nsol)); // Extensive filtering is not implemented in ref python
        log_check("coords_SET1", py_coord);
        log_finish();
        }
    #endif
    return aps_return_t::success;
}

template <class MD, class CD, class C>
typename TLCSolver<MD, CD, C>::coord_t TLCSolver<MD, CD, C>::get_coords() {
    const int sol_idx = *discrete_dofvalue_pointer;
    // fmt::print("Getting soln {}\n", sol_idx);
    boost::numeric::ublas::matrix<double> coords_raw(natoms, 3);
    driver.solution_to_buffer(&coords_raw(0, 0), sol_idx);
    Utils::MapXyzContainer<std::map<int, Utils::uvector3d_t>> coords;
    for (int i = 0; i < natoms; ++i)
        coords.set_atom(idx2atom[i], {coords_raw(i, 0), coords_raw(i, 1), coords_raw(i, 2)});
    return coords;
}




template <class MD, class CD, class C>
typename IdentitySolver<MD, CD, C>::dofs_container_t IdentitySolver<MD, CD, C>::mydofs(py::object& gr, py::list& py_consbonds, const requested_free_t& requested_free) {
    py::list log_atoms;
    for (const auto& py_idx : gr.attr("nodes")())
        log_atoms.append(py_idx);
    assertm(requested_free.size() == 0, fmt::format("IdentitySolver got requests for free DOFs. They had to be filtered out way before this function is called. atoms={}", repr(log_atoms)));
    
    const int nnodes = gr.attr("number_of_nodes")().cast<int>();
    for (const auto& py_idx : gr.attr("nodes")()) {
        auto idx = py::handle(py_idx).cast<int>();
        conformer.set_atom(idx, {0.0, 0.0, 0.0});
    }

    for (const auto& bond : py_consbonds) {
        auto bond_list = py::handle(bond).cast<py::list>();
        auto bond_cpp = bond_list.cast<std::vector<unsigned int>>();
        #ifdef KDMOL_LOG
            if(bond_list.size() != 2)
                std::runtime_error("Unexpected size");
        #endif
        std::vector<unsigned int> sides;
        for (const auto& atom : bond_list) {
            auto nb_atoms_py = py::handle(gr.attr("neighbors")(atom)).cast<py::list>();
            std::vector<unsigned int> good_atoms;
            for (const auto& nbatom : nb_atoms_py)
                if (!bond_list.attr("__contains__")(nbatom).cast<bool>())
                    good_atoms.push_back(nbatom.cast<unsigned int>());
            sides.push_back(*std::min_element(good_atoms.cbegin(), good_atoms.cend()));
        }
        check_dihedrals.push_back({sides[0], bond_cpp[0], bond_cpp[1], sides[1]});
    }
    check_dihedrals.shrink_to_fit();
    checkpoint("idsolver_init");
    #ifdef KDMOL_LOG
        py::list py_atoms_idx;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            py_atoms_idx.append(idx);
        }
        py_atoms_idx.attr("sort")();

        this->py_atoms = py_atoms_idx;

        log_routine("idsolver_init");
        log_in("atoms", py_atoms_idx);
        log_check("consbonds_SET2", py_consbonds);
        log_check("params_SET4", check_dihedrals);
        log_finish();
    #endif
    return check_dihedrals;
}

template <class MD, class CD, class C>
void IdentitySolver<MD, CD, C>::finalize_solver (const int prob_idx,
                                                 mol_data_t* mol_data,
                                                 cyclic_dofs_t* cyclic_dofs)
{
    firstdof_idx = cyclic_dofs->fragment_starts[prob_idx];
    dofvalue_pointer = &cyclic_dofs->dofvalue_container[firstdof_idx];

    for (auto& [atom_idx, value] : conformer)
        value = mol_data->xyz[atom_idx];
    
    checkpoint("idsolver_finalize");
    #ifdef KDMOL_LOG
        for (int i = 0; i < check_dihedrals.size(); ++i) {
            auto myfulldof = check_dihedrals[i];
            Utils::dof_t mydof = {std::get<1>(myfulldof), std::get<2>(myfulldof)};
            auto extfulldof = cyclic_dofs->fulldof_container[firstdof_idx + i];
            auto extdof = cyclic_dofs->dof_container[firstdof_idx + i];
            if ((myfulldof != extfulldof) || (mydof != extdof)) {
                std::runtime_error("Error free");
            }
        }

        py::list py_fixeddih, py_coords;
        for (int i = 0; i < check_dihedrals.size(); ++i) {
            auto myfulldof = check_dihedrals[i];
            Utils::dof_t mydof = {std::get<1>(myfulldof), std::get<2>(myfulldof)};

            py::list newitem;
            newitem.append(py::cast(mydof.first));
            newitem.append(py::cast(mydof.second));
            newitem.append(py::cast(cyclic_dofs->dofvalue_container[firstdof_idx + i]));
            py_fixeddih.append(newitem);
        }

        for (auto& [atom_idx, value] : conformer) {
            py::list newitem;
            newitem.append(atom_idx);
            newitem.append(Utils::vector_to_list(value));
            py_coords.append(newitem);
        }
        
        log_routine("idsolver_finalize");
        log_in("atoms", py_atoms);
        log_check("fixeddih_SET2", py_fixeddih);
        log_check("coords_SET1", py_coords);
        log_finish();
    #endif
}

template <class MD, class CD, class C>
const std::set<int> IdentitySolver<MD, CD, C>::get_atoms_idxs() {
    std::set<int> res;
    for (const auto& [ key, _ ] : conformer)
        res.insert(key);
    return res;
}

template <class MD, class CD, class C>
aps_return_t IdentitySolver<MD, CD, C>::apply_ps(mol_data_t* mol_data) {
    checkpoint("idsolver_aps");
    #ifdef KDMOL_LOG
        py::list py_coords;
        for (auto& [atom_idx, value] : conformer) {
            py::list newitem;
            newitem.append(atom_idx);
            newitem.append(Utils::vector_to_list(value));
            py_coords.append(newitem);
        }
        
        log_routine("idsolver_aps");
        log_in("atoms", py_atoms);
        log_check("coords_SET1", py_coords);
        log_finish();
    #endif
    return aps_return_t::success;
}

template <class MD, class CD, class C>
typename IdentitySolver<MD, CD, C>::coord_t IdentitySolver<MD, CD, C>::get_coords() {
    return conformer;
}



#endif // SOLVERS_H
