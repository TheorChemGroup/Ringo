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

        template <class mol_data_t>
        inline int randomize_ddof (const int max_ddof, mol_data_t& mol_data) {
            discrete_distr = discrete_distr_t(0, max_ddof - 1);
            return discrete_distr(mol_data.pcg_rng);
        }
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

        template <class mol_data_t>
        inline int randomize_ddof (const int max_ddof, mol_data_t& mol_data) {
            return discrete_refiner.randomize_ddof(max_ddof, mol_data);
        }
};


#endif // REFINERS_H
