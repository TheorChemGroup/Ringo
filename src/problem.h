#ifndef PROBLEM_H
#define PROBLEM_H

#include "solvers.h"
#include "utils.h"


class ProblemModel {
    public:
        using nxgraph_t = py::object;
        using bond_container_t = py::list;
        enum class solvertype_t { TLC = 1, Identity = 0, Undefined = -1 };
        solvertype_t solvertype;
        using requested_free_t = std::vector<Utils::dof_t>;
        using requested_free_pyt = py::list;
        using unfulfilled_free_t = std::vector<Utils::dof_t>;

    protected:
        nxgraph_t gr;
        bond_container_t consbonds;
        bool final_flag;
        std::string applicability_warning;
        requested_free_t requested_free;
        unfulfilled_free_t unfulfilled_free;

        inline void add_constraints() {
            for (const auto& py_edge : gr.attr("edges")()) {
                auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
                gr.attr("__getitem__")(edge.first).attr("__getitem__")(edge.second).attr("__setitem__")("type", py::cast(false));
                gr.attr("__getitem__")(edge.first).attr("__getitem__")(edge.second).attr("__setitem__")("unrestrained", py::cast(false));
            }
        }

        inline solvertype_t choose_solver() {
            const auto [ tlc_applicable, warning_message, uncomplied_requests ] = TLCSolverModel::can_apply(this->gr, this->consbonds, this->requested_free);

            // The message that appeared first is the most important and must be recorded
            if ((applicability_warning.size() == 0) && (warning_message.size() > 0) && (final_flag))
                applicability_warning = warning_message;


            if (tlc_applicable) {
                unfulfilled_free = unfulfilled_free_t(uncomplied_requests);
                return solvertype_t::TLC;
            } else {
                this->add_constraints();
                unfulfilled_free = unfulfilled_free_t(this->requested_free);
                return solvertype_t::Identity;
            }
        }


    public:
        ProblemModel () : solvertype(solvertype_t::Undefined), applicability_warning("") { }

        ProblemModel(nxgraph_t gr, bond_container_t consbonds, requested_free_pyt requested_free_py, const bool final_flag) :
            solvertype(solvertype_t::Undefined),
            gr(gr),
            consbonds(consbonds),
            final_flag(final_flag),
            applicability_warning(""),
            requested_free(requested_free_py.cast<requested_free_t>())
        {
            solvertype = choose_solver(); // Enable warnings that TLC cannot be applied 
        }

        ProblemModel(nxgraph_t gr, bond_container_t consbonds, const requested_free_t& requested_free, const bool final_flag) :
            solvertype(solvertype_t::Undefined),
            gr(gr),
            consbonds(consbonds),
            final_flag(final_flag),
            applicability_warning(""),
            requested_free(requested_free)
        {
            solvertype = choose_solver(); // Enable warnings that TLC cannot be applied 
        }

        int method() const { return static_cast<const int>(solvertype); }

        py::bool_ recheck_method() {
            const auto oldmethod = this->solvertype;
            this->solvertype = choose_solver();
            if ((oldmethod == solvertype_t::TLC) && (oldmethod != this->solvertype))
                add_constraints();
            return oldmethod == this->solvertype;
        }

        void record_warning() {
            if (applicability_warning.size() > 0) {
                assertm(final_flag, "record_warning call on non-final ProblemModel");
                Utils::add_message_to_feed(applicability_warning);
            }
        }

        requested_free_t& access_requested_freedofs() { return requested_free; }

        py::list get_unfulfilled_requests_py() const { return py::cast(unfulfilled_free); }
        
        const unfulfilled_free_t& get_unfulfilled_requests() const { return unfulfilled_free; }

        static void add_message_to_feed(const std::string json_contents) { Utils::add_message_to_feed(json_contents); }
        static py::dict get_warning_codes() { return py::cast(Utils::warning_codes); }
};

template <class MD, class DF, class DD, class CD, class C>
class Problem : public ProblemModel {
    using mol_data_t = MD;
    using dofs_t = DF;
    using discrete_dofs_t = DD;
    using cyclic_dofs_t = CD;
    using extcoord_t = C;

    using frame_t = Utils::frame_t;
    using dof_t = Utils::dof_t;
    using fulldof_t = Utils::fulldof_t;

    using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;

    using requested_free_t = std::vector<Utils::dof_t>;
    using requested_free_pyt = py::list;
    using unfulfilled_free_t = std::vector<Utils::dof_t>;

    public:
        using abstractsolver_t = AbstractSolver<MD, CD, C>;
        using tlcsolver_t = TLCSolver<MD, CD, C>;
        using idsolver_t = IdentitySolver<MD, CD, C>;
        using ProblemModel::solvertype_t;
        using marked_dofs_t = std::unordered_map<std::string, std::vector<Utils::fulldof_t>>;
        using dofs_container_t = typename abstractsolver_t::dofs_container_t;
        using owner_atoms_t = typename discrete_dofs_t::owner_atoms_t;

        abstractsolver_t* solver;

    private:
        using ProblemModel::ProblemModel;
        #ifdef VALIDATION
            using validator_t = Utils::LocalGeometryValidator;
            validator_t validator;
        #endif
        
        mol_data_t* mol_data;
        const dofs_t* dofs;
        cyclic_dofs_t* cyclic_dofs;

        active_idxs_t own_dof_idxs, own_discrete_idxs;

    public:
        Problem(ProblemModel::nxgraph_t gr, ProblemModel::bond_container_t consbonds, abstractsolver_t* solver,
                mol_data_t* mol_data, dofs_t* dofs, cyclic_dofs_t* cyclic_dofs, const requested_free_t& requested_free) :
                ProblemModel(gr, consbonds, requested_free, false), // false because all warnings are already recorded in final 'process_seq' call
                mol_data(mol_data),
                dofs(dofs),
                cyclic_dofs(cyclic_dofs),
                solver(solver)
        {
            checkpoint("problem_constr");
            assertm(this->unfulfilled_free.size() == 0, "unfulfilled DOF requirements are not filtered out on Problem obj initialization");
            
            #ifdef KDMOL_LOG
                py::list log_atoms, log_consbonds;
                for (const auto& py_idx : gr.attr("nodes")()) {
                    auto idx = py::handle(py_idx).cast<py::int_>();
                    log_atoms.append(idx);
                }
                log_atoms.attr("sort")();

                for (const auto& item : consbonds)
                    log_consbonds.append(item);
                log_consbonds.attr("sort")();

                if (this->solvertype == solvertype_t::Undefined)
                    throw std::runtime_error(fmt::format("The solver is undefined ({})", static_cast<int>(this->solvertype)));
                const auto solver_int = static_cast<int>(this->solvertype);

                log_routine("problem_constr");
                log_in("atoms", log_atoms);
                log_check("consbonds_SET2", log_consbonds);
                log_check("solver", as_list(solver_int));
                log_finish();
            #endif
            
            if ((this->solvertype == solvertype_t::TLC) && (!dynamic_cast<tlcsolver_t*>(solver))) {
                throw std::runtime_error("Solvertype and provided solver instance do not match");
            } else if ((this->solvertype == solvertype_t::Identity) && (!dynamic_cast<idsolver_t*>(solver))) {
                throw std::runtime_error("Solvertype and provided solver instance do not match");
            }
        }

        std::pair<dofs_container_t, marked_dofs_t> get_dofs(ProblemModel::nxgraph_t gr, ProblemModel::bond_container_t py_consbonds) {
            auto dofs_temp = solver->mydofs(gr, py_consbonds, this->requested_free);

            marked_dofs_t res = {
                                    {"free", dofs_container_t()},
                                    {"dep", dofs_container_t()},
                                    {"fixed", dofs_container_t()},
                                };
            auto dep_bonds = py_consbonds.cast<std::vector<std::pair<unsigned int, unsigned int>>>();
            for (const auto& dof : dofs_temp) {
                auto left = static_cast<unsigned int>(std::get<1>(dof));
                auto right = static_cast<unsigned int>(std::get<2>(dof));
                
                if (!gr.attr("__getitem__")(py::cast(left)).attr("__getitem__")(py::cast(right)).attr("__getitem__")("type").cast<bool>()) {
                    res["fixed"].push_back(dof);
                    continue;
                }

                auto itA = std::find(dep_bonds.begin(), dep_bonds.end(), std::make_pair(left, right));
                auto itB = std::find(dep_bonds.begin(), dep_bonds.end(), std::make_pair(right, left));
                if ((itA != dep_bonds.end()) || (itB != dep_bonds.end())) {
                    res["dep"].push_back(dof);
                    continue;
                }

                res["free"].push_back(dof);
            }

            dofs_container_t dofs;
            dofs.reserve(dofs_temp.size());
            for (const auto& dof : res["free"])
                dofs.push_back(dof);
            for (const auto& dof : res["dep"])
                dofs.push_back(dof);
            for (const auto& dof : res["fixed"])
                dofs.push_back(dof); // Fixed must be at the end!
            
            checkpoint("problem_dofscheck");
            #ifdef KDMOL_LOG
                {
                py::list log_atoms;
                for (const auto& py_idx : gr.attr("nodes")()) {
                    auto idx = py::handle(py_idx).cast<py::int_>();
                    log_atoms.append(idx);
                }
                log_atoms.attr("sort")();
                
                py::list log_free = py::cast(res.at("free")),
                         log_dep = py::cast(res.at("dep")),
                         log_fixed = py::cast(res.at("fixed"));

                log_routine("problem_dofscheck");
                log_in("atoms", log_atoms);
                log_check("dofs_SET4", dofs);
                log_check("free_SET4", log_free);
                log_check("dep_SET4", log_dep);
                log_check("fixed_SET4", log_fixed);
                log_finish();
                }
            #endif

            return std::make_pair(dofs, res);
        }

        const owner_atoms_t get_discrete_dof() const {
            owner_atoms_t res;
            if (this->solvertype == solvertype_t::TLC) {
                const auto atoms = dynamic_cast<tlcsolver_t*>(solver)->get_atoms();
                for (const auto atom : atoms)
                    res.insert(atom);
            }
            return res;
        }

        void finalize_solver (const int frag_idx) {
            solver->finalize_solver(frag_idx, mol_data, cyclic_dofs);
        }
        
        void set_dofs_pointer(cyclic_dofs_t* cyclic_dofs) { this->cyclic_dofs = cyclic_dofs; }

        void set_ddof_pointer (const int prob_idx) {
            assert(this->solvertype != solvertype_t::Identity);
            assert(this->solvertype != solvertype_t::Undefined);
            dynamic_cast<tlcsolver_t*>(solver)->set_ddof_pointer(prob_idx, cyclic_dofs);
        }

        const frame_t get_valid_frame(const int atom) const {
            auto nbs = py::handle(this->gr.attr("neighbors")(atom)).cast<py::list>().cast<std::vector<int>>();
            std::sort(nbs.begin(), nbs.end());
            return {atom, nbs[0], nbs[1]};
        }
        
        const fulldof_t get_valid_dihedral(const dof_t bond) const {
            dof_t sides;

            auto nbs = py::handle(this->gr.attr("neighbors")(bond.first)).cast<py::list>().cast<std::vector<int>>();
            std::sort(nbs.begin(), nbs.end());
            if (nbs[0] == bond.second)
                sides.first = nbs[1];
            else
                sides.first = nbs[0];

            nbs = py::handle(this->gr.attr("neighbors")(bond.second)).cast<py::list>().cast<std::vector<int>>();
            std::sort(nbs.begin(), nbs.end());
            if (nbs[0] == bond.first)
                sides.second = nbs[1];
            else
                sides.second = nbs[0];

            return {sides.first, bond.first, bond.second, sides.second};
        }

        #ifdef VALIDATION
            void build_validator() {
                validator = validator_t(*mol_data, *dofs, solver->get_atoms_idxs());
            }
        #endif

        void index_owned_dofs(const int prob_idx) {
            own_dof_idxs = active_idxs_t();
            auto [ start, end ] = cyclic_dofs->get_fragment_dofs(prob_idx);
            for (int i = start; i < end; ++i)
                if (cyclic_dofs->doftype_container[i] == cyclic_dofs_t::dof_type_t::free)
                    own_dof_idxs.push_back(dofs->free_indexer.at(cyclic_dofs->cp2mol_free[i]));
            own_dof_idxs.shrink_to_fit();

            checkpoint("tlc_refinedofs");
            #ifdef KDMOL_LOG
                if (own_dof_idxs.size() > 0) {
                    assert(solvertype != solvertype_t::Identity);
                    auto atoms = dynamic_cast<tlcsolver_t*>(solver)->get_atoms();

                    py::list py_atoms = py::cast(atoms), py_dofdata;
                    py_atoms.attr("sort")();

                    for (int i = start; i < end; ++i)
                        if (cyclic_dofs->doftype_container[i] == cyclic_dofs_t::dof_type_t::free) {
                            auto mol_idx = cyclic_dofs->cp2mol_free[i];
                            py_dofdata.append(py::make_tuple(
                                        dofs->side_container.at(mol_idx).first,
                                        dofs->bond_container.at(mol_idx).first,
                                        dofs->bond_container.at(mol_idx).second,
                                        dofs->side_container.at(mol_idx).second
                            ));
                        }
                    log_routine("tlc_refinedofs");
                    log_in("atoms", py_atoms);
                    log_check("refineps_SET4", py_dofdata);
                    log_finish();
                }
            #endif

            own_discrete_idxs = active_idxs_t();
            if (solvertype != solvertype_t::Identity) {
                own_discrete_idxs.push_back(cyclic_dofs->cp2mol_discrete.at(prob_idx));
                own_discrete_idxs.shrink_to_fit();
            } else {
                assert(own_dof_idxs.size() == 0);
            }
            solver->index_owned_dofs(own_dof_idxs, own_discrete_idxs);
        }

        const auto get_owned_dofs() {
            return std::make_pair(own_dof_idxs, own_discrete_idxs);
        }

        aps_return_t apply_ps() {
            auto result = solver->apply_ps(mol_data);
            if (result != aps_return_t::success)
                return result;
                
            #ifdef VALIDATION
                auto coords = solver->get_coords();
                if (!validator.validate(coords, *mol_data, *dofs)) {
                    mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
                    return aps_return_t::validationfail;
                    // throw Utils::GeomFailure("Validation of problem conformation failed");
                }
            #endif
            return aps_return_t::success;
        }

        auto get_coords() {
            return solver->get_coords();
        }
};

#endif // PROBLEM_H
