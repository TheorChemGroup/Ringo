#ifndef REFINERS_H
#define REFINERS_H


template <class DF, class AI>
class DihedralRefiner {
    private:
        using dofs_t = DF;
        using dihedral_distr_t = std::uniform_real_distribution<double>;
        using active_idxs_t = AI;
        dofs_t* dofs;
        dihedral_distr_t dihedral_distr;
        active_idxs_t active_idxs;

        using custom_sampling_limits_t = typename dofs_t::custom_sampling_limits_t;
        using custom_limits_data_t = typename dofs_t::custom_limits_data_t;
        bool all_default_limits;
        std::set<int> dofs_default_limits;
        using custom_distrs_t = std::unordered_map<int, dihedral_distr_t>;
        custom_distrs_t custom_distrs;

        #ifdef KDMOL_LOG
            py::list input_iterations;
            int iteration_idx;
        #endif

    public:
        DihedralRefiner () { }

        DihedralRefiner (dofs_t* dofs) :
            dofs(dofs),
            dihedral_distr(dihedral_distr_t(-M_PI, M_PI)),
            all_default_limits(true)
            #ifdef KDMOL_LOG
                , iteration_idx(-1)
            #endif
        { }

        #ifdef KDMOL_LOG
            void assign_unittests(const py::list& input_iterations) {
                py::set input_setcheck;
                for (const auto [ i, _ ] : dofs->free_indexer)
                    input_setcheck.add(py::cast(std::make_tuple(
                        dofs->side_container[i].first,
                        dofs->bond_container[i].first,
                        dofs->bond_container[i].second,
                        dofs->side_container[i].second
                    )));
                // Check that free dofs match the provided testing data
                py::set ref_setcheck;
                for (const auto& item : input_iterations[0]["param_data"])
                    ref_setcheck.add(py::handle(item).cast<py::tuple>());
                if (!input_setcheck.attr("__eq__")(ref_setcheck).cast<bool>())
                    throw std::runtime_error(fmt::format("{} != {}", repr(input_setcheck), repr(ref_setcheck)));
                
                this->input_iterations = input_iterations;
            }
        #endif

        void rebuild_active() { // This method is only called by user
            #ifdef KDMOL_LOG
                iteration_idx += 1; // -1 -> 0 on first call or just move to the next step
                auto iteration = py::handle(input_iterations[iteration_idx]).cast<py::dict>();
                auto param_data = py::handle(iteration["param_data"]).cast<py::list>();
                auto param_values = py::handle(iteration["param_values"]).cast<py::list>();
                assert(param_data.size() == param_values.size());

                active_idxs = active_idxs_t(param_data.size());
                for (int i = 0; i < param_data.size(); ++i) {
                    auto param_tuple = py::handle(param_data[i]).cast<py::tuple>();
                    bool dof_found = false;
                    for(const auto [dof_idx, free_idx] : dofs->free_indexer) {
                        auto dof_tuple = py::cast(std::make_tuple(
                            dofs->side_container[dof_idx].first,
                            dofs->bond_container[dof_idx].first,
                            dofs->bond_container[dof_idx].second,
                            dofs->side_container[dof_idx].second
                        ));
                        
                        if (py::handle(param_tuple.attr("__eq__")(dof_tuple)).cast<bool>()) {
                            active_idxs[i] = free_idx;
                            dof_found = true;
                        }
                    }
                    if (!dof_found)
                        fmt::print("Can't find DOF {}\n", repr(param_tuple));
                    assert(dof_found);
                }
            #else
                active_idxs = active_idxs_t(dofs->free_indexer.size());
                for (int j = 0; j < dofs->free_indexer.size(); ++j)
                    active_idxs[j] = j;
            #endif
        }

        inline void set_active(const std::vector<int>& dih_dofs) { active_idxs = dih_dofs; }
        inline int get_num_active() const { return active_idxs.size(); }

        template <class mol_data_t>
        void set_dofs(mol_data_t& mol_data) {
            #ifdef KDMOL_LOG
                auto iteration = py::handle(input_iterations[iteration_idx]).cast<py::dict>();
                auto param_data = py::handle(iteration["param_data"]).cast<py::list>();
                auto param_values = py::handle(iteration["param_values"]).cast<py::list>();
                assert(param_data.size() == param_values.size());
                assert(param_data.size() == active_idxs.size());

                for (int i = 0; i < param_data.size(); ++i)
                    dofs->dof_rawpointer[active_idxs[i]] = param_values[i].cast<double>();
            #else
                if (all_default_limits) {
                    for (const auto dof_idx : active_idxs)
                        dofs->dof_rawpointer[dof_idx] = dihedral_distr(mol_data.pcg_rng);
                } else {
                    for (const auto dof_idx : active_idxs) {
                        auto it = dofs_default_limits.find(dof_idx);
                        if (it != dofs_default_limits.end()) {
                            dofs->dof_rawpointer[dof_idx] = dihedral_distr(mol_data.pcg_rng);
                        } else {
                            dofs->dof_rawpointer[dof_idx] = Utils::standardize_dihedral(custom_distrs.at(dof_idx)(mol_data.pcg_rng));
                        }
                    }
                }
            #endif
        }

        void customize_sampling_limits(const custom_sampling_limits_t& custom_sampling_limits) {
            all_default_limits = false;
            dofs_default_limits.clear();
            for (const auto [ dof_idx, free_idx ] : dofs->free_indexer) { // The left one is full_dof
                auto it = custom_sampling_limits.find(free_idx);
                if (it == custom_sampling_limits.end()) {
                    dofs_default_limits.insert(free_idx);
                } else {
                    const auto dihedral = std::make_tuple(
                            dofs->side_container[dof_idx].first,
                            dofs->bond_container[dof_idx].first,
                            dofs->bond_container[dof_idx].second,
                            dofs->side_container[dof_idx].second
                        );
                }
            }

            for (const auto [ free_idx, limits ] : custom_sampling_limits) {
                custom_distrs[free_idx] = dihedral_distr_t(limits.first, limits.second);
            }
        }

        void customize_sampling_limits(const custom_limits_data_t& custom_limits_data) {
            all_default_limits = false;
            dofs_default_limits.clear();

            for (const auto [ dihedral, limits ] : custom_limits_data) {
                int found_free_idx = -1;
                for (const auto [ dof_idx, free_idx ] : dofs->free_indexer) {
                    const auto current_dihedral = std::make_tuple(
                        dofs->side_container[dof_idx].first,
                        dofs->bond_container[dof_idx].first,
                        dofs->bond_container[dof_idx].second,
                        dofs->side_container[dof_idx].second
                    );

                    if (current_dihedral == dihedral) {
                        found_free_idx = free_idx;
                        break;
                    } else if ((std::get<1>(current_dihedral) == std::get<1>(dihedral)) && (std::get<2>(current_dihedral) == std::get<2>(dihedral))) {
                        fmt::print("Mismatch of two PSs: {} vs. {}", repr(dihedral), repr(current_dihedral)); std::cout << std::endl;
                        throw std::runtime_error(fmt::format("Mismatch of two PSs: {} vs. {}", repr(dihedral), repr(current_dihedral)));
                    }
                }
                if (found_free_idx == -1) {
                    fmt::print("Mismatch of two PSs: couldn't find the dihedral {}", repr(dihedral)); std::cout << std::endl;
                }
                assertm(found_free_idx != -1, fmt::format("Mismatch of two PSs: couldn't find the dihedral {}", repr(dihedral)));

                custom_distrs[found_free_idx] = dihedral_distr_t(limits.first, limits.second);
            }

            for (const auto [ _, free_idx ] : dofs->free_indexer) {
                auto it = custom_distrs.find(free_idx);
                if (it == custom_distrs.end())
                    dofs_default_limits.insert(free_idx);
            }
        }

        custom_limits_data_t get_custom_limits() const {
            custom_limits_data_t custom_limits_data;
            for (const auto& [ free_idx, distr ] : custom_distrs) {
                int dof_idx = -1;
                for (const auto [ cdof, cfree ] : dofs->free_indexer)
                    if (free_idx == cfree) {
                        dof_idx = cdof;
                        break;
                    }
                const auto dihedral = std::make_tuple(
                        dofs->side_container[dof_idx].first,
                        dofs->bond_container[dof_idx].first,
                        dofs->bond_container[dof_idx].second,
                        dofs->side_container[dof_idx].second
                    );
                custom_limits_data[dihedral] = std::make_pair(distr.min(), distr.max());
            }
            return custom_limits_data;
        }

        const active_idxs_t& get_active() const { return active_idxs; }
        inline void set_successful() noexcept { active_idxs.clear(); }
};


template <class DD, class AI>
class DiscreteRefiner {
    private:
        using discrete_dofs_t = DD;
        using discrete_distr_t = std::uniform_int_distribution<int>;
        using active_idxs_t = AI;
        discrete_dofs_t* discrete_dofs;
        discrete_distr_t discrete_distr;
        active_idxs_t active_idxs;

        #ifdef KDMOL_LOG
            py::list input_iterations;
            int iteration_idx;
        #endif

    public:
        
        DiscreteRefiner () { }

        DiscreteRefiner (discrete_dofs_t* discrete_dofs) :
            discrete_dofs(discrete_dofs)
            #ifdef KDMOL_LOG
                , iteration_idx(-1)
            #endif
        { }

        #ifdef KDMOL_LOG
            void assign_unittests(const py::list& input_iterations) {
                // Check that discrete DOFs match the provided testing data
                py::set input_discrete_set;
                for (const auto& item : discrete_dofs->owner_container) {
                    auto list_item = py::handle(py::cast(item)).cast<py::list>();
                    list_item.attr("sort")();
                    input_discrete_set.add(list_item.cast<py::tuple>());
                }
                py::set ref_discrete_set;
                for (const auto& item : input_iterations[0]["ddof_data"]) {
                    auto list_item = py::handle(item).cast<py::list>();
                    list_item.attr("sort")();
                    ref_discrete_set.add(list_item.cast<py::tuple>());
                }
                if (!input_discrete_set.attr("__eq__")(ref_discrete_set).cast<bool>())
                    throw std::runtime_error(fmt::format("{} != {}", repr(input_discrete_set), repr(ref_discrete_set)));
                
                this->input_iterations = input_iterations;
            }
        #endif

        void rebuild_active() { // This method is only called by user
            #ifdef KDMOL_LOG
                iteration_idx += 1; // -1 -> 0 on first call or just move to the next step
                auto iteration = py::handle(input_iterations[iteration_idx]).cast<py::dict>();
                auto ddof_data = py::handle(iteration["ddof_data"]).cast<py::list>();
                auto ddof_values = py::handle(iteration["ddof_values"]).cast<py::list>();
                assert(ddof_data.size() == ddof_values.size());

                active_idxs = active_idxs_t(ddof_data.size());
                for (int i = 0; i < ddof_data.size(); ++i) {
                    auto ref_ddof = py::handle(ddof_data[i]).cast<py::set>();
                    bool ddof_found = false;
                    for (int j = 0; j < discrete_dofs->size(); ++j) {
                        auto cur_ddof = py::cast(discrete_dofs->owner_container[j]);
                        if (py::handle(ref_ddof.attr("__eq__")(cur_ddof)).cast<bool>()) {
                            active_idxs[i] = j;
                            ddof_found = true;
                        }
                    }
                    assert(ddof_found);
                }
            #else
                active_idxs = active_idxs_t(discrete_dofs->size());
                for (int j = 0; j < discrete_dofs->size(); ++j)
                    active_idxs[j] = j;
            #endif
        }

        inline void set_active(const std::vector<int>& discr_dofs) { active_idxs = discr_dofs; }
        inline int get_num_active() const { return active_idxs.size(); }

        template <class mol_data_t>
        void set_dofs (mol_data_t& mol_data) {
            #ifdef KDMOL_LOG
                auto iteration = py::handle(input_iterations[iteration_idx]).cast<py::dict>();
                auto ddof_data = py::handle(iteration["ddof_data"]).cast<py::list>();
                auto ddof_values = py::handle(iteration["ddof_values"]).cast<py::list>();
                assert(ddof_data.size() == ddof_values.size());
                assert(ddof_data.size() == active_idxs.size());

                for (int i = 0; i < ddof_data.size(); ++i)
                    discrete_dofs->dof_rawpointer[active_idxs[i]] = ddof_values[i].cast<int>();
            #else
                for (const auto ddof_idx : active_idxs)
                    discrete_dofs->dof_rawpointer[ddof_idx] = -1;
            #endif
        }

        template <class mol_data_t>
        inline int randomize_ddof (const int max_ddof, mol_data_t& mol_data) {
            #ifdef KDMOL_LOG
                throw std::runtime_error("Shouldn't go here");
            #endif
            discrete_distr = discrete_distr_t(0, max_ddof - 1);
            return discrete_distr(mol_data.pcg_rng);
        }

        const active_idxs_t& get_active() const { return active_idxs; }
        inline void set_successful() noexcept { active_idxs.clear(); }
};


template <class DF, class DD>
class Refiners {
    public:
        using active_idxs_t = std::vector<int>;

    private:
        using dofs_t = DF;
        using discrete_dofs_t = DD;

        using dihedral_refiner_t = DihedralRefiner<dofs_t, active_idxs_t>;
        using discrete_refiner_t = DiscreteRefiner<discrete_dofs_t, active_idxs_t>;

        dihedral_refiner_t dihedral_refiner;
        discrete_refiner_t discrete_refiner;

        using custom_sampling_limits_t = typename dofs_t::custom_sampling_limits_t;
        using custom_limits_data_t = typename dofs_t::custom_limits_data_t;
 
    public:
        Refiners () { }
 
        Refiners (dofs_t* dofs, discrete_dofs_t* discrete_dofs) :
            dihedral_refiner(dihedral_refiner_t(dofs)),
            discrete_refiner(discrete_refiner_t(discrete_dofs))
        { }

        void rebuild_active() {
            dihedral_refiner.rebuild_active();
            discrete_refiner.rebuild_active();
        }

        inline void set_active(const std::vector<int>& dih_dofs, const std::vector<int>& discr_dofs) {
            dihedral_refiner.set_active(dih_dofs);
            discrete_refiner.set_active(discr_dofs);
        }

        template <class mol_data_t>
        void set_dofs(mol_data_t& mol_data) {
            dihedral_refiner.set_dofs(mol_data);
            discrete_refiner.set_dofs(mol_data);
        }

        #ifdef KDMOL_LOG
            void assign_unittests(const py::list& input_iterations) {
                dihedral_refiner.assign_unittests(input_iterations);
                discrete_refiner.assign_unittests(input_iterations);
            }
        #endif

        inline const std::pair<const active_idxs_t&, const active_idxs_t&> get_active () const {
            return {dihedral_refiner.get_active(), discrete_refiner.get_active()};
        }

        inline void set_successful() noexcept {
            dihedral_refiner.set_successful();
            discrete_refiner.set_successful();
        }

        inline int get_num_active() const { return dihedral_refiner.get_num_active() + discrete_refiner.get_num_active(); }

        template <class mol_data_t>
        inline int randomize_ddof (const int max_ddof, mol_data_t& mol_data) {
            return discrete_refiner.randomize_ddof(max_ddof, mol_data);
        }

        inline void customize_sampling_limits(const custom_sampling_limits_t& custom_sampling_limits) {
            dihedral_refiner.customize_sampling_limits(custom_sampling_limits);
        }

        inline void customize_sampling_limits(const custom_limits_data_t& custom_sampling_limits) {
            dihedral_refiner.customize_sampling_limits(custom_sampling_limits);
        }

        inline custom_limits_data_t get_custom_limits() const {
            return dihedral_refiner.get_custom_limits();
        }

        // TODO const int size() { ... }
};


#endif // REFINERS_H
