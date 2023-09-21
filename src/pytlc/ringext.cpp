#include "ringext.h"

using TlcUtils::AtomContainerDynamic;
using TlcUtils::BondContainerDynamic;
using TlcUtils::AtomItem;
using TlcUtils::AxialAtomContainer;
using TlcUtils::AtomGroupSimple;
using TlcUtils::BondGroupSimple;


std::pair<bool, std::vector<Utils::dof_t>> TLCExtSolver::check_applicable(const int n_atoms, const std::vector<unsigned int>& constr_bonds,
                                    const std::vector<unsigned int>& requested_free_idxs)
{
    AtomGroupSimple catoms(n_atoms);
    BondGroupSimple cbonds(n_atoms);
    TlcUtils::connect(catoms, cbonds);

    auto constrained_bonds = BondContainerDynamic(constr_bonds, cbonds);
    for(auto bond : constrained_bonds)
        bond->free = false;
    
    std::vector<Utils::dof_t> uncomplied_dofs;

    if (requested_free_idxs.size() > 0) {
        auto requested_free = BondContainerDynamic(requested_free_idxs, cbonds);

        AtomContainerDynamic proper_atoms; // Both bonds are easily dependent
        AtomContainerDynamic onereqested_atoms; // One bond was requested to be free
        AtomContainerDynamic tworequested_atoms; // Both bonds were requested to be free

        // fmt::print("Enforcing {}", repr(requested_free_idxs)); std::cout << std::endl;
        for(int i = 1; i < n_atoms + 1; ++i) {
            auto& atom = catoms[i];
            const auto& left_bond = atom.get_left_bond();
            const auto& right_bond = atom.get_right_bond();
            if ((!left_bond.free) || (!right_bond.free))
                continue; // At least one neighb bond is not free => not a potental axial atom
            
            const bool left_requested = requested_free.contains(left_bond);
            const bool right_requested = requested_free.contains(right_bond);
            if (left_requested && right_requested)
                tworequested_atoms += atom;
            else if (left_requested || right_requested)
                onereqested_atoms += atom;
            else
                proper_atoms += atom;
        }

        AxialAtomContainer atompos_main, atompos;
        atompos_main.set(&catoms[0], 0);
        atompos_main.set(&catoms[0], 1);
        atompos_main.set(&catoms[0], 2);
        AtomContainerDynamic pa = catoms.collect_proper_atoms();
        auto pa_size = pa.size();
        double distance_metric = 0;
        int request_metric = std::numeric_limits<int>::min();
        for (int i = 0; i < pa_size; ++i) {
            atompos.set(&pa[i], 0);
            for (int j = i + 1; j < pa_size; ++j) { // QUESTIONABLE AF!!!
                atompos.set(&pa[j], 1);
                for (int k = j + 1; k < pa_size; ++k) {
                    atompos.set(&pa[k], 2);
                    if (atompos.is_sensible(catoms)) {
                        int current_request_metric = 0;
                        for (auto& atom : atompos) {
                            if (onereqested_atoms.contains(atom))
                                current_request_metric -= 1;
                            else if (tworequested_atoms.contains(atom))
                                current_request_metric -= 2;
                        }

                        const auto current_distance_metric = atompos.get_metric(catoms);
                        const bool requests_improved = (request_metric < current_request_metric);
                        const bool requests_equal = (request_metric == current_request_metric);
                        const bool distance_improved = (distance_metric < current_distance_metric);
                        if (requests_improved || (requests_equal && distance_improved)) {
                            atompos_main = atompos;
                            distance_metric = current_distance_metric;
                            request_metric = current_request_metric;
                        }
                    }
                }
            }
        }

        bool applicable = !((atompos_main[0] == catoms[0]) && (atompos_main[1] == catoms[0]));
        if (applicable) {
            for (auto* atom : atompos_main) {
                if (proper_atoms.contains(atom))
                    continue;
                else if (onereqested_atoms.contains(atom)) {
                    auto& left_bond = atom->get_left_bond();
                    auto& right_bond = atom->get_right_bond();
                    const bool left_requested = requested_free.contains(left_bond);
                    const bool right_requested = requested_free.contains(right_bond);
                    if (left_requested)
                        uncomplied_dofs.push_back(left_bond.atom_idxs());
                    else
                        uncomplied_dofs.push_back(right_bond.atom_idxs());
                } else { // if (tworequested_atoms.contains(atom))
                    auto& left_bond = atom->get_left_bond();
                    auto& right_bond = atom->get_right_bond();
                    uncomplied_dofs.push_back(left_bond.atom_idxs());
                    uncomplied_dofs.push_back(right_bond.atom_idxs());
                }
            }
        }

        return std::make_pair(applicable, uncomplied_dofs);

    } else {
        AxialAtomContainer atompos_main;
        bool first_added = false;
        for(int i = 1; i < n_atoms + 1; ++i) {
            AtomItem& atom = catoms[i];
            if(atom.can_be_axial()) {
                if (i == 1)
                    first_added = true;
                if (i == n_atoms && first_added)
                    break;
                atompos_main.add_atom(atom);
                if(atompos_main.is_full())
                    break;
                i += 1; // Skip the next atom
            }
        }

        bool applicable = atompos_main.is_full() && atompos_main.is_sensible(catoms);
        return std::make_pair(applicable, uncomplied_dofs);
    }
}


TLCExtSolver::propose_dofs_t TLCExtSolver::propose_dofs(const int n_atoms, const std::vector<unsigned int> constr_bonds,
                                std::vector<int> idx2atom, const std::vector<unsigned int>& requested_free_idxs)
{
    AtomGroupSimple catoms(n_atoms);
    BondGroupSimple cbonds(n_atoms);
    TlcUtils::connect(catoms, cbonds);
    #ifdef EXT_LOG
        TlcUtils::check_topology(catoms, cbonds);
    #endif
    
    auto constrained_bonds = BondContainerDynamic(constr_bonds, cbonds);
    auto requested_free = BondContainerDynamic(requested_free_idxs, cbonds);

    // auto idx2atom_py = py::cast(idx2atom);
    for(auto bond : constrained_bonds) {
        bond->free = false;
        // if (idx2atom_py.attr("__contains__")(py::cast(3)).cast<bool>() && idx2atom_py.attr("__contains__")(py::cast(13)).cast<bool>()) {
        //     fmt::print("Rofl = {} - {}\n", py::repr(py::cast(idx2atom[bond->get_left_atom().get_num()])).cast<std::string>(),
        //                                    py::repr(py::cast(idx2atom[bond->get_right_atom().get_num()])).cast<std::string>());
        // }
    }
    // if (idx2atom_py.attr("__contains__")(py::cast(3)).cast<bool>() && idx2atom_py.attr("__contains__")(py::cast(4)).cast<bool>()) {
    //     throw std::runtime_error(fmt::format("Walk = {}", py::repr(idx2atom_py).cast<std::string>()));
    // }

    // Ensure that requested_free does not overlap with constrained_bonds, then temporarily set them as fixed
    for(auto bond : requested_free) {
        assertm(bond->free, "Requested DOF is already set as fixed");
        bond->free = false;
    }

    AxialAtomContainer atompos_main;
    bool first_added = false;
    for(int i = 1; i < n_atoms + 1; ++i) {
        AtomItem& atom = catoms[i];
        if(atom.can_be_axial()) {
            if (i == 1)
                first_added = true;
            if (i == n_atoms && first_added)
                break;
            atompos_main.add_atom(atom);
            if(atompos_main.is_full())
                break;
            i += 1; // Skip the next atom
        }
    }
    // if(!(atompos_main.is_full()) || !(atompos_main.is_sensible(catoms)))
    //     throw std::runtime_error("TLC can not be applied to the current ring and arrangment of fixed bonds");

    double metric, curmetric = atompos_main.get_metric(catoms);
    #ifdef EXT_LOG
        log_metric_call(catoms, atompos_main, curmetric);
    #endif

    AtomContainerDynamic pa = catoms.collect_proper_atoms();
    // std::vector<int> test(pa.size());
    // for(int i = 0; i < pa.size(); ++i)
    //     test[i] = pa[i].get_num();
    // fmt::print("PA = {}\n", py::repr(py::cast(test)).cast<std::string>());

    auto pa_size = pa.size();
    AxialAtomContainer atompos;
    axial_atoms_t testA;
    for(int f = 0; f < 3; ++f)
        testA[f] = atompos_main[f].get_num();
    py::list myatoms_py = py::cast(testA);
    // fmt::print("Start = {} ", py::repr(myatoms_py).cast<std::string>());
    for (int i = 0; i < pa_size; ++i) {
        atompos.set(&pa[i], 0);
        for (int j = i + 2; j < pa_size; ++j) {
            atompos.set(&pa[j], 1);
            for (int k = j + 2; k < pa_size; ++k) {
                atompos.set(&pa[k], 2);

                axial_atoms_t testB;
                for(int f = 0; f < 3; ++f)
                    testB[f] = atompos[f].get_num();
                py::list curatoms_py = py::cast(testB);
                // fmt::print("Test = {} ", py::repr(py::cast(testB)).cast<std::string>());

                if (atompos.is_sensible(catoms)) {
                    metric = atompos.get_metric(catoms);
                    if ((metric + 0.00001 < curmetric) || ((abs(curmetric - metric) < 0.00001) && py::handle(curatoms_py.attr("__lt__")(myatoms_py)).cast<bool>())) {
                        // if (metric + 0.00001 < curmetric)
                        //     fmt::print("Accept metric {} < {}\n", metric, curmetric);
                        // else
                        //     fmt::print("Accept list {} < {}\n", py::repr(curatoms_py).cast<std::string>(), py::repr(myatoms_py).cast<std::string>());
                        myatoms_py = py::cast(testB);
                        atompos_main = atompos;
                        curmetric = metric;
                    }
                }
            }
        }
    }

    // Switch requested free bonds back to be free
    for(auto bond : requested_free)
        bond->free = true;

    #ifdef EXT_LOG
        TlcUtils::check_topology(catoms, cbonds);
    #endif

    // Deduce solved bonds
    for (auto &bond : cbonds)
        bond.solved = false;
    for(auto& atom : atompos_main)
        atom->make_bonds_solved();
    
    fulldof_data_t dof_data = {
        { "variable", fulldof_vector_t() },
        { "solved", fulldof_vector_t() },
        { "fixed", fulldof_vector_t() }
    };

    for(auto& bond : cbonds) {
        std::pair<unsigned int, unsigned int> my_atoms = bond.atom_idxs();
        auto new_item = std::make_tuple((bond.get_left_atom() - 1).get_num(), my_atoms.first, my_atoms.second, (bond.get_right_atom() + 1).get_num());
        if (bond.free && !bond.solved)
            dof_data["variable"].push_back(fulldof_atoms_t(new_item));
        else if(bond.free && bond.solved)
            dof_data["solved"].push_back(fulldof_atoms_t(new_item));
        else
            dof_data["fixed"].push_back(fulldof_atoms_t(new_item));
    }

    axial_atoms_t axial_atoms;
    for(int i = 0; i < 3; ++i)
        axial_atoms[i] = atompos_main[i].get_num();
    
    // auto idx2atom_py = py::cast(idx2atom);
    // if (idx2atom_py.attr("__contains__")(py::cast(3)).cast<bool>() && idx2atom_py.attr("__contains__")(py::cast(13)).cast<bool>()) {
    //     fmt::print("Rofl = {}\n", py::repr(py::cast(axial_atoms)).cast<std::string>());
    //     throw std::runtime_error(fmt::format("Walk = {}", py::repr(idx2atom_py).cast<std::string>()));
    // }

    return std::make_pair(dof_data, axial_atoms);
}

TLCExtSolver::dof_data_t TLCExtSolver::condense_dofs(const fulldof_data_t& fulldofs) {
    TLCExtSolver::dof_data_t dof_data;
    for (const auto& [key, fulldof_list] : fulldofs) {
        auto& dof_list = dof_data[key];
        for (const auto fulldof : fulldof_list) {
            dof_list.push_back(std::make_pair(std::get<1>(fulldof), std::get<2>(fulldof)));
        }
    }
    return dof_data;
}

TLCExtSolver::TLCExtSolver(const int n_atoms, const dof_data_t& dof_data, const axial_atoms_t& axial_atoms,
                           const value_container_t& length_data, const value_container_t& vangle_data,
                           const value_container_t& fixed_dihedral_data)
{
    catoms = TlcUtils::AtomGroup(n_atoms);
    cbonds = TlcUtils::BondGroup(n_atoms);
    TlcUtils::connect(catoms, cbonds);
    #ifdef EXT_LOG
        TlcUtils::check_topology(catoms, cbonds);
    #endif
    
    // Default: free=true and solved=false
    for(auto& idxpair : dof_data.at("solved"))
        catoms[idxpair.first].get_right_bond().solved = true;

    for(auto& idxpair : dof_data.at("fixed"))
        catoms[idxpair.first].get_right_bond().free = false;
    
    cbonds.markdown_dihedrals(dof_data, catoms);

    TlcUtils::AtomContainerStatic<3> axatoms(axial_atoms, catoms);
    catoms.init_tlc_map(axial_atoms);

    const std::array<std::pair<unsigned int, unsigned int>, 3> rf_idsx = {std::make_pair(0, 1),
                                                                          std::make_pair(1, 2),
                                                                          std::make_pair(2, 0)};
    for (int i = 0; i < 3; ++i)
        rfs[i].initialize(axatoms[rf_idsx[i].first], axatoms[rf_idsx[i].second], catoms, cbonds);
    
    cbonds.set_length_data(length_data);
    catoms.set_vangle_data(vangle_data);
    cbonds.set_fixed_dihedral_data(fixed_dihedral_data);

    // Partial construction of cbonds.tlc_lengths and catoms.tlc_vangles
    for(auto& rf : rfs)
        rf.write_tlc_data_init(catoms, cbonds);
}

py::object TLCExtSolver::initialize(const int n_atoms, const py::dict& py_dof_data, const py::list& py_axial_atoms, const py::kwargs& kwargs) {
    auto dof_data = py_dof_data.cast<dof_data_t>();
    if(py_axial_atoms.size() != 3)
        throw std::runtime_error("Expected a list of THREE axial atoms");
    auto axial_atoms = py_axial_atoms.cast<axial_atoms_t>();
    
    if (kwargs.attr("__contains__")("params").cast<bool>()) {
        auto param_list = kwargs["params"].cast<py::list>();
        if (param_list.size() != 3)
            throw std::runtime_error(fmt::format("Expected 3 arguments in param list. Got {}", param_list.size()));

        auto length_data = param_list[0].cast<value_container_t>();
        auto vangle_data = param_list[1].cast<value_container_t>();
        for(auto& item : vangle_data)
            item *= DEG2RAD;
        auto fixed_dihedral_data = param_list[2].cast<value_container_t>();
        for(auto& item : fixed_dihedral_data)
            item *= DEG2RAD;
        return py::cast(TLCExtSolver(n_atoms, dof_data, axial_atoms, length_data, vangle_data, fixed_dihedral_data));
        
    } else if (kwargs.attr("__contains__")("xyz").cast<bool>())
        throw std::runtime_error("Not implemented");
    else
        throw std::runtime_error("Either 'xyz' or 'params' keyword argument must be provided to this function");
}

int TLCExtSolver::solve(const std::vector<double>& free_tangles) {
    cbonds.set_variable_dihedral_data(free_tangles);

    // #ifdef EXT_LOG
    //     log_routine("solveA");
    //     log_in("n_atoms", catoms.size());
    //     log_check("vangles_SET1", catoms.list_vangles());
    //     log_finish();
    //     log_in("free_tangles", free_tangles);
    //     log_check("lengths", cbonds.list_lengths()); // Not correct!!!!
    //     log_check("dihedrals", cbonds.list_dihedrals_full());
    // #endif

    for(auto& rf : rfs)
        rf.write_tlc_data_work(catoms, cbonds);

    // #ifdef EXT_LOG
    //     log_routine("solveB");
    //     log_in("n_atoms", catoms.size());
    //     log_in("free_tangles", free_tangles);
    //     log_check("bonds", cbonds.repr_tlc_lengths());
    //     log_check("vangles", catoms.repr_tlc_vangles());
    //     log_check("dihedrals", cbonds.repr_tlc_dihedrals());
    //     log_finish();
    // #endif

    int nsol = tlcsolver.solve_internal(cbonds.tlc_lengths, catoms.tlc_vangles, cbonds.tlc_dihedrals);

    // #ifdef EXT_LOG
    //     log_routine("solveC");
    //     log_in("n_atoms", catoms.size());
    //     log_in("free_tangles", free_tangles);
    //     log_check("ret_value", nsol);
    //     log_finish();
    //     tangles_backup = py::repr(py::cast(free_tangles)).cast<std::string>();
    // #endif

    return nsol;
}

py::int_ TLCExtSolver::py_solve(const py::list& py_tangles) {
    auto free_tangles = py_tangles.cast<std::vector<double>>();
    for(auto& item : free_tangles)
        item *= DEG2RAD;
    return solve(free_tangles);
}

void TLCExtSolver::set_solution(py::array_t<double> solution, const int& sol_idx) {
    py::buffer_info buffer = solution.request();
    double *buffer_raw = static_cast<double*>(buffer.ptr);
    const auto natoms = catoms.size();
    if (buffer.ndim != 2)
        throw std::runtime_error("Input should be 2-D NumPy array");

    if ((solution.shape()[0] != natoms) || (solution.shape()[1] != 3))
        throw std::runtime_error(fmt::format("Dimensions of np.array are not right: ({}; {}) instead of ({}; 3)", solution.shape()[0], solution.shape()[1], natoms));

    solution_to_buffer(buffer_raw, sol_idx);
}

void TLCExtSolver::solution_to_buffer(double* buffer, const int& sol_idx) {
    #ifdef EXT_LOG
        using FMatrix = boost::numeric::ublas::fixed_matrix<double, 9, 3>;
        FMatrix conf;
        tlcsolver.solution_to_buffer(&conf(0, 0), sol_idx);

        // Logger::log_rstart("setsolConf");
        // Logger::log(fmt::format(" IN:: n_atoms = {}", catoms.size()));
        // Logger::log(fmt::format(" IN:: sol_idx = {}", sol_idx));
        // Logger::log(fmt::format(" IN:: free_tangles = {}", tangles_backup));
        // Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix(conf)));
        // Logger::log_rend("setsolConf");

        for(unsigned int i = 0; i < 3 * catoms.size(); ++i)
            buffer[i] = 0.0;
    #endif
    
    tlcsolver.solution_to_buffer(buffer, sol_idx, catoms.tlc_map);
    
    #ifdef EXT_LOG
        // Logger::log_rstart("setsolTLC");
        // Logger::log(fmt::format(" IN:: n_atoms = {}", catoms.size()));
        // Logger::log(fmt::format(" IN:: free_tangles = {}", tangles_backup));
        // Logger::log(fmt::format(" IN:: sol_idx = {}", sol_idx));
        // Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix_buffer(buffer, catoms.size(), 3)));
        // Logger::log_rend("setsolTLC");
    #endif

    for(auto& rf : rfs)
        rf.reconstruct(buffer, catoms, cbonds);
    
    #ifdef EXT_LOG
        // Logger::log_rstart("setsolMain");
        // Logger::log(fmt::format(" IN:: n_atoms = {}", catoms.size()));
        // Logger::log(fmt::format(" IN:: free_tangles = {}", tangles_backup));
        // Logger::log(fmt::format(" IN:: sol_idx = {}", sol_idx));
        // Logger::log(fmt::format(" CHECK:: xyzs = {}", TlcUtils::repr_matrix_buffer(buffer, catoms.size(), 3)));
        // Logger::log_rend("setsolMain");
    #endif
}
