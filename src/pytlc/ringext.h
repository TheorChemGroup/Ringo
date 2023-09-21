#ifndef TLCEXTSOLVER_H
#define TLCEXTSOLVER_H

#include "tlcmain.h"
#include "rigidfrag.h"
#include "utils.h"


class TLCExtSolver {
    public:
        using axial_atoms_t = std::array<unsigned int, 3>;
        using value_container_t = std::vector<double>;

        using fulldof_atoms_t = Utils::fulldof_t;
        using fulldof_vector_t = std::vector<fulldof_atoms_t>;
        using fulldof_data_t = std::unordered_map<std::string, fulldof_vector_t>;
        using propose_dofs_t = std::pair<fulldof_data_t, axial_atoms_t>;

        using dof_atoms_t = Utils::dof_t;
        using dof_vector_t = std::vector<dof_atoms_t>;
        using dof_data_t = std::unordered_map<std::string, dof_vector_t>;

        TLCExtSolver() noexcept { }

        TLCExtSolver(const int n_atoms, const dof_data_t& dof_data, const axial_atoms_t& axial_atoms,
                     const value_container_t& length_data, const value_container_t& vangle_data,
                     const value_container_t& fixed_dihedral_data);

        int solve(const std::vector<double>& tangles);
        py::int_ py_solve(const py::list& tangles);
        void set_solution(py::array_t<double> solution, const int& sol_idx);

        // --- STATIC STUFF ---
        static std::pair<bool, std::vector<Utils::dof_t>> check_applicable(const int n_atoms, const std::vector<unsigned int>& constr_bonds, const std::vector<unsigned int>& requested_free_idxs);
        static propose_dofs_t propose_dofs(const int n_atoms, const std::vector<unsigned int> constr_bonds, std::vector<int> idx2atom, const std::vector<unsigned int>& requested_free_idxs); // OPT idx2atom is needed for debugging
        // Factory of TLCExtSolver objects for Python interface (Remove??)
        static py::object initialize(const int n_atoms, const py::dict& py_dof_data, const py::list& py_axial_atoms, const py::kwargs& kwargs);
        static dof_data_t condense_dofs(const fulldof_data_t& fulldofs);
        // --------------------

        #ifdef EXT_LOG
            py::str get_log() const
            { return py::cast(Utils::Logger::get_log()); }
        #endif

        void solution_to_buffer(double* buffer, const int& sol_idx);

        bool are_some_removed() const noexcept { return tlcsolver.some_removed; }

    private:
        #ifdef EXT_LOG
            static void log_metric_call(const TlcUtils::AtomGroupSimple& catoms, const TlcUtils::AxialAtomContainer& atompos, const double& curmetric) {
                Logger::log_rstart("calcmetric");
                Logger::log(" IN:: n_atoms = {}", catoms.size());
                Logger::log(" IN:: a = {}", atompos[0].get_num());
                Logger::log(" IN:: b = {}", atompos[1].get_num());
                Logger::log(" IN:: c = {}", atompos[2].get_num());
                Logger::log(" CHECK:: metric = {}", curmetric);
                Logger::log_rend("calcmetric");
            }
            std::string tangles_backup;
        #endif

        TLCMain tlcsolver;
        RigidFragContainer rfs;
        TlcUtils::AtomGroup catoms;
        TlcUtils::BondGroup cbonds;
};


#endif // TLCEXTSOLVER_H