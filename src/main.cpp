#include "molecule.h"
#include "problem.h"
#include "pytlc/ringext.h"
#include "pytlc/tlcmain.h"

#ifdef BUILD_PYXYZ
    #include "pyxyz/confpool.h"
    #include "pyxyz/molproxy.h"
#endif


PYBIND11_MODULE(ringo_base, m) {
    py::class_<Molecule>(m, "Molecule")
        .def(py::init<>())
        .def("sdf_constructor", &Molecule::constructor<const std::string, false>)
        .def("graph_constructor", &Molecule::constructor<py::dict, true>)
        .def("set_requested_dofs", &Molecule::set_requested_dofs_py)
        .def("apply_ps", &Molecule::apply_ps_py)
        .def("get_xyz", &Molecule::get_xyz_py)
        .def("get_symbols", &Molecule::get_symbols_py)
        .def("reconfigure", &Molecule::reconfigure)
        .def("get_ps", &Molecule::get_ps_py)
        .def("get_discrete_ps", &Molecule::get_discrete_ps_py)
        .def("require_best_sequence", &Molecule::require_best_sequence)
        .def("customize_sampling", &Molecule::customize_sampling_py)
        .def("get_biggest_ringfrag_atoms", &Molecule::get_biggest_ringfrag_atoms)
        .def("molgraph_access", &Molecule::molgraph_access_py)
        .def("get_num_flexible_rings", &Molecule::get_num_flexible_rings)
        .def("get_num_rigid_rings", &Molecule::get_num_rigid_rings)
        #ifdef KDMOL_LOG
            .def("get_log", &Molecule::get_log)
        #endif
        ;
        
    py::class_<TLCExtSolver>(m, "TLCExtSolver")
        .def(py::init<>())
        // .def_static("propose_dofs", &TLCExtSolver::propose_dofs_py)
        .def_static("initialize", &TLCExtSolver::initialize)
        .def("solve", &TLCExtSolver::py_solve)
        .def("set_solution", &TLCExtSolver::set_solution)
        #ifdef EXT_LOG
            .def("get_log", &TLCExtSolver::get_log)
        #endif
        ;
    
    py::class_<TLCMain>(m, "TLCMain")
        .def(py::init<>())
        .def("solve", &TLCMain::solve)
        .def("set_solution", &TLCMain::set_solution)
        #ifdef TLC_LOG
            .def("get_log", &TLCMain::get_log)
        #endif
        ;
    
    py::class_<ProblemModel>(m, "Problem")
        .def(py::init<py::object, py::list, py::list, const bool>())
        .def("recheck_method", &ProblemModel::recheck_method)
        .def("record_warning", &ProblemModel::record_warning)
        .def("get_unfulfilled_requests", &ProblemModel::get_unfulfilled_requests_py)
        .def_static("add_message_to_feed", &ProblemModel::add_message_to_feed)
        .def_static("get_warning_codes", &ProblemModel::get_warning_codes)
        .def_property_readonly("method", &ProblemModel::method)
        ;

    #ifdef KDMOL_LOG
        m.def("log", &Utils::Logger::log_py);
        m.def("log_rstart", &Utils::Logger::log_rstart_py);
        m.def("log_rend", &Utils::Logger::log_rend_py);
    #endif

    #ifdef KDMOL_LOG
        m.def("run_unit_tests", &run_unit_tests);
    #endif

    #ifdef BUILD_PYXYZ
        // Classes related to pyxyz - Confpool and MolProxy
        py::class_<Confpool>(m, "Confpool")
            .def(py::init<const py::object&>(), py::arg("copy") = py::none())
            .def("__getitem__", &Confpool::__getitem__)
            .def("__setitem__", &Confpool::__setitem__)
            .def("__getattr__", &Confpool::__getattr__)
            .def("__setattr__", &Confpool::__setattr__)
            .def("__delitem__", &Confpool::__delitem__)
            .def("__len__", &Confpool::__len__)
            .def("__iter__", &Confpool::__iter__)
            .def("__next__", &Confpool::__next__)
            .def("include_from_file", &Confpool::include_from_file)
            .def("include_from_xyz", &Confpool::include_from_xyz_py)
            .def("include_subset", &Confpool::include_subset)
            .def("float_from_descr", &Confpool::float_from_descr)
            .def("filter", &Confpool::filter, py::arg(), py::arg("inplace") = true)
            .def("count", &Confpool::count)
            .def("upper_cutoff", &Confpool::upper_cutoff, py::arg(), py::arg(), py::arg("inplace") = true)
            .def("lower_cutoff", &Confpool::lower_cutoff, py::arg(), py::arg(), py::arg("inplace") = true)
            .def("sort", &Confpool::sort, py::arg(), py::arg("ascending") = true, py::arg("inplace") = true)
            .def("generate_connectivity", &Confpool::generate_connectivity)
            .def("generate_isomorphisms", &Confpool::generate_isomorphisms)
            .def("get_num_isomorphisms", &Confpool::get_num_isomorphisms)
            .def("rmsd_filter", &Confpool::rmsd_filter)
            .def("get_rmsd_matrix", &Confpool::get_rmsd_matrix)
            .def("get_connectivity", &Confpool::get_connectivity)
            .def("save", &Confpool::save)
            .def("as_table", &Confpool::as_table)
            ;
    
        py::class_<MolProxy>(m, "MolProxy")
            .def("__getitem__", &MolProxy::__getitem__)
            .def("__setitem__", &MolProxy::__setitem__)
            .def("__getattr__", &MolProxy::__getattr__)
            .def("__setattr__", &MolProxy::__setattr__)
            .def("l", &MolProxy::l)
            .def("v", &MolProxy::v)
            .def("z", &MolProxy::z)
            .def("rmsd", &MolProxy::rmsd)
            .def("get_connectivity", &MolProxy::get_connectivity)
        ;
    #endif

    // Record build options
    #ifdef BUILD_PYXYZ
        m.attr("use_pyxyz") = true;
    #else
        m.attr("use_pyxyz") = false;
    #endif

    #ifdef BUILD_OVERLAP_DETECTION
        m.attr("use_overlap_detection") = true;
        m.def("get_vdw_radii", &Utils::get_vdw_radii);
        m.def("set_vdw_radii", &Utils::set_vdw_radii);
        m.def("set_radius_multiplier", &Utils::set_radius_multiplier);
    #else
        m.attr("use_overlap_detection") = false;
    #endif

    m.def("get_status_feed", &Utils::get_status_feed);
    m.def("clear_status_feed", &Utils::clear_status_feed);
    m.attr("warning_codes") = Utils::warning_codes;
    m.attr("build_flags") = RINGO_BUILDFLAGS;
}
