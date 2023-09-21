#ifndef RESOLVER_H
#define RESOLVER_H

#include "utils.h"

template <class MD, class AD, class CD>
class FusionResolver {
    public:
        using mol_data_t = MD;
        using assembly_data_t = AD;
        using cyclic_dofs_t = CD;
        using fulldof_t = Utils::fulldof_t;
        using dofpointer_t = typename cyclic_dofs_t::dofvalue_container_t::pointer;
        using dofconstpointer_t = typename cyclic_dofs_t::dofvalue_container_t::const_pointer;
        using coord_container_t = typename assembly_data_t::coord_container_t;
        using poly_t = typename mol_data_t::poly_t;

        FusionResolver () { }

        FusionResolver (const fulldof_t ind_dih, dofconstpointer_t ind_val,
                        const int depdof_idx, dofpointer_t dep_val,
                         mol_data_t* mol_data,
                         assembly_data_t& assembly_data,
                        cyclic_dofs_t& cyclic_dofs);
    
        void configure(const fulldof_t dep_dih_raw);
        void resolve_dependence();

    private:
        // Pointers
         poly_t *left_poly, *right_poly;
        dofpointer_t dep_val;
        dofconstpointer_t ind_val;

        fulldof_t ind_dof;
        bool calc_dih;
        double tor_sum;
    public:
        coord_container_t *coords;
        #ifdef KDMOL_DEBUG
            fulldof_t dep_dof;
        #endif
};

template <class MD, class AD, class CD>
FusionResolver<MD,AD,CD>::FusionResolver (const fulldof_t ind_dih, dofconstpointer_t ind_val,
                                          const int depdof_idx, dofpointer_t dep_val,
                                           mol_data_t* mol_data,
                                          assembly_data_t& assembly_data,
                                          cyclic_dofs_t& cyclic_dofs) :
    left_poly(&mol_data->polys[std::get<1>(ind_dih)]),
    right_poly(&mol_data->polys[std::get<2>(ind_dih)]),
    coords(&assembly_data.extg_coords),
    dep_val(dep_val), ind_val(ind_val), 
    ind_dof(ind_dih), calc_dih(ind_val == nullptr)
    #ifdef KDMOL_DEBUG
        , dep_dof(cyclic_dofs.fulldof_container[depdof_idx])
    #endif
{
    #ifdef KDMOL_DEBUG
        const auto dep_bond = cyclic_dofs.dof_container[depdof_idx];
        const auto ind_bond = std::make_pair(std::get<1>(ind_dih), std::get<2>(ind_dih));
        const auto ind_bond_rev = std::make_pair(std::get<2>(ind_dih), std::get<1>(ind_dih));
        if ((dep_bond != ind_bond) && (dep_bond != ind_bond_rev))
            throw std::runtime_error(fmt::format("Mismatch on FusionResolver construction. ind_dih={}, dep_dih={}", 
                                                 py::handle(py::repr(py::cast(ind_dih))).cast<std::string>(),
                                                 py::handle(py::repr(py::cast(cyclic_dofs.fulldof_container[depdof_idx]))).cast<std::string>()));
    #endif

    configure(cyclic_dofs.fulldof_container[depdof_idx]);
}

template <class MD, class AD, class CD>
void FusionResolver<MD,AD,CD>::configure(const fulldof_t dep_dih_raw) {
    fulldof_t dep_dih;
    auto dep_bond = std::make_pair(std::get<1>(dep_dih_raw), std::get<2>(dep_dih_raw));
    auto ind_bond = std::make_pair(std::get<1>(ind_dof), std::get<2>(ind_dof));
    if (dep_bond != ind_bond)
        dep_dih = std::make_tuple(std::get<3>(dep_dih_raw), std::get<2>(dep_dih_raw), std::get<1>(dep_dih_raw), std::get<0>(dep_dih_raw));
    else
        dep_dih = dep_dih_raw;
    
    auto left_impr = left_poly->get_sph_angle(std::get<0>(dep_dih), std::get<2>(ind_dof), std::get<0>(ind_dof));
    auto right_impr = right_poly->get_sph_angle(std::get<3>(dep_dih), std::get<1>(ind_dof), std::get<3>(ind_dof));
    tor_sum = left_impr + right_impr;

    checkpoint("resolver_config");
    #ifdef KDMOL_LOG
        int dof_present;
        if (calc_dih)
            dof_present = 0;
        else
            dof_present = 1;
            
        log_routine("resolver_config");
        log_in("dep_atoms", dep_dih_raw);
        log_in("ind_atoms", ind_dof);
        log_check("left_impr", as_list(left_impr));
        log_check("right_impr", as_list(right_impr));
        log_check("tor_sum", as_list(tor_sum));
        log_check("dof_present", as_list(dof_present));
        log_finish();
    #endif
}

template <class MD, class AD, class CD>
void FusionResolver<MD,AD,CD>::resolve_dependence() {
    double ind_value;
    if (calc_dih) {
        ind_value = coords->get_dihedral(std::get<0>(ind_dof), std::get<1>(ind_dof), std::get<2>(ind_dof), std::get<3>(ind_dof));
    } else {
        ind_value = *ind_val;
    }

    *dep_val = tor_sum + ind_value;

    while (*dep_val < -M_PI)
        *dep_val += 2*M_PI;
    while (*dep_val > M_PI)
        *dep_val -= 2*M_PI;
}

#endif // RESOLVER_H
