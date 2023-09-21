#ifndef CYCLICPART_H
#define CYCLICPART_H

#include "geomunit.h"
#include "fusion_resolver.h"
#include "fused_assembler.h"
#include "generic_fragment.h"
#include "problem.h"
#include "utils.h"

// Data:
// A: pointer to molecule's vector of indepentent DOFs
// B: vector of all dofs grouped by CPs

// First of all, copy from A to B all free DOFs
// LBs know where to look in B OR what dihedral to compute, and put resolved DOFs in B
// Problem objects look for their segment in B

// Info about linking bonds: generate_ig -> attrs in LG edges -> CXX version of LG ->
// -> Problem Constructor -> Solver Constructor -> All kinds of DOFs -> Construct B -> Construct linking bonds

// Dihedral resolution requires only fusion via bonds
// Linking requires fusion via bonds AND atom (spiro)


template <class P, class C>
struct AssemblyData {
    template <class T> using container_t = std::vector<T>;
    using problem_t = P;
    using coord_container_t = C;

    using link_t = std::pair<int, int>;
    using idx_t = int;
    using idx2link_map_t = std::vector<link_t>;
    idx2link_map_t idx2link_map;
    using link2idx_map_t = std::map<link_t, idx_t>;
    link2idx_map_t link2idx_map;

    using atom_container_t = container_t<int>;
    using ring_atoms_t = container_t<atom_container_t>;
    ring_atoms_t ring_atoms;

    using solvertype_t = ProblemModel::solvertype_t;
    using abstractsolver_t = typename problem_t::abstractsolver_t;
    using tlcsolver_t = typename problem_t::tlcsolver_t;
    using idsolver_t = typename problem_t::idsolver_t;
    using tlc_solvers_t = container_t<tlcsolver_t>;
    using id_solvers_t = container_t<idsolver_t>;
    tlc_solvers_t tlc_solvers;
    id_solvers_t id_solvers;

    using inner_t = std::pair<int, int>;
    using inner_container_t = std::unordered_map<int, inner_t>;
    coord_container_t extg_coords;

    using idx2seq_t = std::vector<int>;
    idx2seq_t idx2seq;

    AssemblyData () { }

    void set_nsolvers(const int ntlc, const int nid) {
        tlc_solvers.reserve(ntlc);
        id_solvers.reserve(nid);
    }

    abstractsolver_t* create_solver(solvertype_t solvertype) {
        switch(solvertype) {
            case solvertype_t::TLC:
                tlc_solvers.push_back(tlcsolver_t());
                return static_cast<abstractsolver_t*>(&tlc_solvers.back());
                break;
            case solvertype_t::Identity:
                id_solvers.push_back(idsolver_t());
                return static_cast<abstractsolver_t*>(&id_solvers.back());
                break;
            default:
                throw std::runtime_error("!!!");
                break;
        }
    }

    void finalize() {
        ring_atoms.shrink_to_fit();
    }
};

struct CyclicDofs {
    using value_t = double;
    using index_t = int;
    using fouratoms_t = std::tuple<index_t, index_t, index_t, index_t>;
    template <class T> using container_t = std::vector<T>; // Container must be ordered
    template <class A, class B> using map_t = std::unordered_map<A, B>;

    using dof_t = Utils::dof_t;
    using fulldof_t = Utils::fulldof_t;
    using dof_container_t = container_t<dof_t>;
    using fulldof_container_t = container_t<fulldof_t>;
    dof_container_t dof_container;
    fulldof_container_t fulldof_container;

    using dofvalue_container_t = container_t<value_t>;
    dofvalue_container_t dofvalue_container;
    
    enum dof_type_t : int { free = 2, dep = 1, fixed = 0 };
    using doftype_container_t = container_t<dof_type_t>;
    doftype_container_t doftype_container;

    using fragment_starts_t = container_t<index_t>;
    fragment_starts_t fragment_starts;

    using exposability_container_t = container_t<bool>;
    exposability_container_t exposability_container;

    using cp2mol_map_t = map_t<int,int>;
    cp2mol_map_t cp2mol_full, cp2mol_free, cp2mol_discrete;
    const double* molvalues_start;
    int* molvalues_discrete_start;

    CyclicDofs() { }

    void new_dof(const fulldof_t& dofatoms, const dof_type_t dof_type, const value_t value, const bool exposable) {
        fulldof_container.push_back(dofatoms);
        dof_container.push_back({std::get<1>(dofatoms), std::get<2>(dofatoms)});
        doftype_container.push_back(dof_type);
        dofvalue_container.push_back(value);
        exposability_container.push_back(exposable);
    }

    #ifdef KDMOL_LOG
        fulldof_container_t get_dofdata(const dof_type_t curtype) {
            fulldof_container_t res;
            for (int i = 0; i < fulldof_container.size(); ++i) {
                if (doftype_container.at(i) == curtype)
                    res.push_back(fulldof_container.at(i));
            }
            return res;
        }
    #endif

    void finalize() {
        dof_container.shrink_to_fit();
        fulldof_container.shrink_to_fit();
        dofvalue_container.shrink_to_fit();
        doftype_container.shrink_to_fit();
        exposability_container.shrink_to_fit();

        nfixed = 0;
        nfree = 0;
        ndep = 0;
        for (const auto& item : doftype_container) {
            if (item == fixed)
                nfixed++;
            else if (item == dep)
                ndep++;
            else
                nfree++;
        }
        nexposed = 0;
        for (const auto& item : exposability_container)
            if (item) nexposed++;
    }

    std::pair<int, int> get_fragment_dofs(const int frag_idx) {
        auto start = fragment_starts.at(frag_idx);
        int end;
        if (fragment_starts.size() != frag_idx + 1)
            end = fragment_starts.at(frag_idx + 1);
        else
            end = dof_container.size();
        return std::make_pair(start, end);
    }

    inline int dofindex_in_frag(fulldof_t& search_dof) const {
        int count = std::count(fulldof_container.cbegin(), fulldof_container.cend(), search_dof);

        if (count != 1)
            throw std::runtime_error(fmt::format("The DOF {} is not unique in this cyclic part", repr(search_dof), repr(fulldof_container)));

        auto found_iter = std::find(fulldof_container.cbegin(), fulldof_container.cend(), search_dof);
        const int index = std::distance(fulldof_container.cbegin(), found_iter);
        return index;
    }

    int dofindex_in_frag(const fulldof_t& search_dof, const int frag_index) const {
        int end;
        if (fragment_starts.size() != frag_index + 1)
            end = fragment_starts[frag_index + 1];
        else
            end = dof_container.size();

        int founddof_idx = -1;
        for (int i = fragment_starts[frag_index]; i < end; ++i) {
            if (fulldof_container[i] == search_dof) {
                founddof_idx = i;
                break;
            }
        }

        if (founddof_idx == -1)
            throw std::runtime_error(fmt::format("Bug found (cp.dep_index 2) dof={}", repr(search_dof)));
        if (doftype_container[founddof_idx] == dof_type_t::free)
            throw std::runtime_error(fmt::format("Free DOF (cp.dep_index 2) dof={}", repr(search_dof)));
        return founddof_idx;
    }

    inline int get_nfree() const noexcept { return nfree; }
    
    inline int get_ndep() const noexcept { return ndep; }
    
    inline int get_nfixed() const noexcept { return nfixed; }
    
    inline int size() const noexcept { return fulldof_container.size(); }

    template <class M> using mol_dofs_t = std::vector<std::tuple<fulldof_t, M, double>>;
    
    template <class M>
    mol_dofs_t<M> get_dofs() const { // Return only fixed and free DOFs
        mol_dofs_t<M> res;
        res.reserve(nfixed + nfree);
        for (int i = 0; i < fulldof_container.size(); ++i)
            if (exposability_container[i]) {
                M moltype;
                if (doftype_container[i] == dof_type_t::fixed)
                    moltype = M::fixed;
                else
                    moltype = M::free;
                res.push_back({fulldof_container[i], moltype, dofvalue_container[i]});
            }
        return res;
    }

    private:
        int nfixed, ndep, nfree;

    public:
        int nexposed, ndiscrete;
};


struct LgVertexProps
{ };

struct LgEdgeProps
{
    int link_idx;
    int resolver_start;
    int n_resolvers;
};

struct DgVertexProps
{ };

struct DgEdgeProps
{
    int link_idx;
};

template<class G, class MD, class DF, class DD, class D>
class CyclicPart : public GenericFragment<G, MD, DF, DD, D> {
    using parent_t = GenericFragment<G, MD, DF, DD, D>;
    public:
        using mol_subgraph_t = G;
        using main_graph_t = G;
        using mol_data_t = MD;
        using dofs_t = DF;
        using discrete_dofs_t = DD;
        using dfg_t = D;
        using pygraph_t = G;
        using fulldof_t = Utils::fulldof_t;

        using cyclic_dofs_t = CyclicDofs;
        using dof_type_t = typename cyclic_dofs_t::dof_type_t;

        using coord_container_t = Utils::MapXyzContainer<std::map<int, Utils::uvector3d_t>>;

        using problem_t = Problem<mol_data_t, dofs_t, discrete_dofs_t, cyclic_dofs_t, coord_container_t>;
        using problem_container_t = std::vector<problem_t>;
        using solvertype_t = typename problem_t::solvertype_t;

        using assembly_data_t = AssemblyData<problem_t, coord_container_t>;

        using build_seq_t = std::vector<int>;

        using lgvertex_t = LgVertexProps;
        using lgedge_t = LgEdgeProps;
        using lg_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, lgvertex_t, lgedge_t>;

        using dgvertex_t = DgVertexProps;
        using dgedge_t = DgEdgeProps;
        using dg_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, dgvertex_t, dgedge_t>;

        using geomunit_t = GeomUnit<mol_data_t, dofs_t, discrete_dofs_t, assembly_data_t>;
        using geomunit_container_t = std::vector<geomunit_t>;

        using resolver_t = FusionResolver<mol_data_t, assembly_data_t, cyclic_dofs_t>;
        using resolver_container_t = std::vector<resolver_t>;
        using assembler_t = PolycycleAssembler<mol_data_t, assembly_data_t>;
        using assembler_container_t = std::vector<assembler_t>;

        using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;

        #ifdef VALIDATION
            using validator_t = Utils::LocalGeometryValidator;
            validator_t conf_validator, ext_validator, local_validator;
        #endif

        CyclicPart() { }
        
        CyclicPart(mol_data_t* mol_data, dofs_t* dofs, discrete_dofs_t* discrete_dofs, mol_subgraph_t& gr, const std::string& filename);

        void build_infrastructure(mol_subgraph_t& gr);

        void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg, const int given_index) override;

        void configure_assembly(const main_graph_t& molgr, const dfg_t& dfg) override
        { this->configure_assembly(molgr, dfg, this->index); }

        aps_return_t apply_ps();

        void update_configuration(const dfg_t& dfg);

        inline int get_some_atom() const
        { return py::list(gr.attr("nodes")())[0].cast<int>(); }

        inline auto get_dofs() const;

        void map_dofvalues(const int start, const int size);
        void map_discrete_dofvalues(const int start, const int size);

        using owner_container_t = typename discrete_dofs_t::owner_container_t;
        const owner_container_t get_discrete_dofs() const;

        void index_owned_dofs();

        int get_num_tlcsolvers() const noexcept;
        int get_num_idsolvers() const noexcept;

    private:
        template <class C>
        void add_dofset(const std::pair<typename problem_t::dofs_container_t,
                                        typename problem_t::marked_dofs_t>& newdofs,
                        const int idx, C& coord, const py::object& gr);

        mol_subgraph_t gr;
        // pygraph_t ig, lg;
        lg_t lg;
        dg_t dg;
        int at0_idx, at1_idx, at2_idx;
        bool recalc_frame;
        build_seq_t build_sequence;
        assembly_data_t assembly_data;
        cyclic_dofs_t cyclic_dofs;
        problem_container_t problem_container;
        geomunit_container_t geomunit_container;
        resolver_container_t resolver_container;
        assembler_container_t assembler_container;

        active_idxs_t own_dof_idxs, own_discrete_idxs;

        int index;
        std::string filename;
};


template<class G, class MD, class DF, class DD, class D>
CyclicPart<G, MD, DF, DD, D>::CyclicPart(mol_data_t* mol_data, dofs_t* dofs, discrete_dofs_t* discrete_dofs, mol_subgraph_t& gr, const std::string& filename) :
    parent_t(mol_data, dofs, discrete_dofs),
    gr(gr),
    at0_idx(-1),
    at1_idx(-1),
    at2_idx(-1),
    index(-1),
    filename(filename)
{

    checkpoint("cp_constr");
    #ifdef KDMOL_LOG
        py::list py_atoms_idx, py_atoms_xyz, py_bonds_num, py_atoms_sym, py_bonds_type;
        for (const auto& py_idx : gr.attr("nodes")()) {
            py::list new_xyz;
            auto idx = py::handle(py_idx).cast<int>();
            py_atoms_idx.append(idx);
            auto xyz = this->mol_data->xyz[idx];
            new_xyz.append(xyz[0]);
            new_xyz.append(xyz[1]);
            new_xyz.append(xyz[2]);
            py_atoms_xyz.append(new_xyz);
            py_atoms_sym.append(mol_data->atom_symbols[idx]);
        }
        py_atoms_idx.attr("sort")();
        
        for (const auto& py_edge : gr.attr("edges")()) {
            auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
            auto edge_idx = mol_data->bond_data(edge).index;
            py_bonds_num.append(py_edge);
            py_bonds_type.append(py::cast(mol_data->bondtypes[edge_idx]));
        }

        log_routine("cp_constr");
        log_in("atoms", py_atoms_idx);
        log_check("bonds_SET2", py_bonds_num);
        log_finish();
    #endif

    build_infrastructure(gr);
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::build_infrastructure(mol_subgraph_t& gr) {
    for (const auto& py_edge : gr.attr("edges")()) {
        auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
        auto edge_idx = this->mol_data->bond_data(edge).index;
        gr.attr("__getitem__")(edge.first).attr("__getitem__")(edge.second).attr("__setitem__")("type", py::cast(this->mol_data->bondtypes[edge_idx]));
    }
    auto ig = Utils::PyModules::get_pyutils()->attr("generate_ig")(gr);

    checkpoint("cp_checkIG");
    #ifdef KDMOL_LOG
        {
        py::list py_atoms_idx, py_ignodes, py_igedges, py_combonds, py_comatoms, py_ndof, py_atoms_part;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            py_atoms_idx.append(idx);
        }
        py_atoms_idx.attr("sort")();
        
        for (const auto& py_idx : ig.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            py_ignodes.append(idx);

            py_ndof.append(py::make_tuple(
                py_idx, ig.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("NDOF")
            ));

            py::list newitem;
            newitem.append(idx);
            py::list atom_idxs;
            for (const auto& py_atom_idx : ig.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("graph").attr("nodes")()) {
                atom_idxs.append(py_atom_idx);
            }
            atom_idxs.attr("sort")();
            newitem.append(atom_idxs);
            py_atoms_part.append(newitem);
        }
        py_ignodes.attr("sort")();
        
        for (const auto& py_edge : ig.attr("edges")()) {
            auto edgelist = py::handle(py_edge).cast<py::list>();
            py_igedges.append(edgelist);

            py_combonds.append(py::make_tuple(
                edgelist[0], edgelist[1],
                ig.attr("__getitem__")(edgelist[0]).attr("__getitem__")(edgelist[1]).attr("__getitem__")("common_bonds")
            ));
            
            py_comatoms.append(py::make_tuple(
                edgelist[0], edgelist[1],
                ig.attr("__getitem__")(edgelist[0]).attr("__getitem__")(edgelist[1]).attr("__getitem__")("link")
            ));
        }

        log_routine("cp_checkIG");
        log_in("atoms", py_atoms_idx);
        log_check("nodes", py_ignodes);
        log_check("edges_SET2", py_igedges);
        log_check("commonbonds_SET2", py_combonds);
        log_check("commonatoms_SET2", py_comatoms);
        log_check("ndofs_SET1", py_ndof);
        log_check("atoms_SET1", py_atoms_part);
        log_finish();
        }
    #endif

    std::vector<Utils::dof_t> requested_free;
    if (this->dofs->requested_free.size() > 0)
        for (const auto& bond : this->dofs->requested_free)
            if (py::handle(gr.attr("has_edge")(bond.first, bond.second)).cast<bool>())
                requested_free.push_back(bond);
    auto py_requested_free = py::cast(requested_free);

    auto py_build_sequence = Utils::PyModules::get_pyutils()->attr("get_strategy")(ig, py::cast(this->filename), py_requested_free, this->dofs->require_best_sequence);
    build_sequence = py::handle(py_build_sequence).cast<build_seq_t>();

    auto pylg = py::handle(Utils::PyModules::get_pyutils()->attr("process_seq")(py_build_sequence, ig, py_requested_free)).cast<pygraph_t>();

    // TODO Might move this to AssemblyData methods
    const auto lg_size = py::handle(pylg.attr("number_of_nodes")()).cast<int>();
    lg = lg_t(lg_size);
    auto& link2idx = assembly_data.link2idx_map;
    auto& idx2link = assembly_data.idx2link_map;
    idx2link.reserve(lg_size);

    // Remove all unfulfilled DOF requests out of the requested_free
    for (int idx = 0; idx < lg_size; ++idx) {
        auto py_idx = py::cast(idx);
        auto prob_model = py::handle(pylg.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("ikprob")).cast<ProblemModel>();
        auto& requested_free = prob_model.access_requested_freedofs();
        const auto& unfulfilled_free = prob_model.get_unfulfilled_requests();

        for (const auto& bond : unfulfilled_free) {
            auto it = std::find(requested_free.begin(), requested_free.end(), bond);
            if (it == requested_free.end()) {
                it = std::find(requested_free.begin(), requested_free.end(), std::make_pair(bond.second, bond.first));
                assertm(it != requested_free.end(), "Cannot find the DOF request");
            }
            requested_free.erase(it);
        }
        pylg.attr("nodes").attr("__getitem__")(py_idx).attr("__setitem__")("ikprob", py::cast(prob_model));
    }

    assembly_data.idx2seq = typename assembly_data_t::idx2seq_t(lg_size);
    for (int i = 0; i < lg_size; ++i)
        assembly_data.idx2seq[build_sequence[i]] = i;

    int edge_idx = 0;
    for (const auto& py_edge : pylg.attr("edges")()) {
        auto [ vA, vB ] = py::handle(py_edge).cast<std::pair<int, int>>();
        lg_t::edge_descriptor edge_descr = add_edge(vA, vB, lg).first;
        lg[edge_descr].link_idx = edge_idx;
        lg[edge_descr].n_resolvers = 0;
        idx2link.push_back({vA, vB});
        link2idx[{vA, vB}] = edge_idx;
        edge_idx++;
    }

    assembly_data.ring_atoms = typename assembly_data_t::ring_atoms_t(lg_size);
    problem_container = problem_container_t(lg_size);
    cyclic_dofs.fragment_starts = typename cyclic_dofs_t::fragment_starts_t(lg_size);
    auto ntlc = 0;
    auto nid = 0;
    for (int idx = 0; idx < lg_size; ++idx) {
        const auto py_idx = py::cast(idx);
        auto prob_model = py::handle(pylg.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("ikprob")).cast<ProblemModel>();
        switch(prob_model.solvertype) {
            case solvertype_t::TLC:
                ntlc++;
                break;
            case solvertype_t::Identity:
                nid++;
                break;
            default:
                throw std::runtime_error("!!!");
                break;
        }
    }
    // Vectors of solvers MUST be created before we start ditrubuting pointers
    assembly_data.set_nsolvers(ntlc, nid);

    // Create solvers and DOF set
    for (int idx = 0; idx < lg_size; ++idx) {
        const auto py_idx = py::cast(idx);
        auto g = ig.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("graph");
        assembly_data.ring_atoms[idx] = py::handle(g.attr("nodes")()).cast<py::list>().cast<typename assembly_data_t::atom_container_t>();

        auto linking_bonds = py::handle(pylg.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("linking_bonds")).cast<py::list>();
        
        auto prob_model = py::handle(pylg.attr("nodes").attr("__getitem__")(py_idx).attr("__getitem__")("ikprob")).cast<ProblemModel>();
        if (prob_model.solvertype == solvertype_t::Identity) {
            assertm(!this->mol_data->xyzless_init, fmt::format("XYZ-less conformer generation is impossible because the ring {} cannot be generated through inverse kinematics", repr(assembly_data.ring_atoms[idx])));
        }
        auto* solver_pointer = assembly_data.create_solver(prob_model.solvertype);
        auto& requested_free = prob_model.access_requested_freedofs();
        problem_container[idx] = problem_t(g, linking_bonds, solver_pointer, this->mol_data, this->dofs, &cyclic_dofs, requested_free);
        add_dofset(problem_container[idx].get_dofs(g, linking_bonds), idx, this->mol_data->xyz, g);
    }
    cyclic_dofs.finalize(); // Only now can distribute pointers
    
    for (int idx = 0; idx < lg_size; ++idx)
        problem_container[idx].finalize_solver(idx);
    
    cyclic_dofs.ndiscrete = assembly_data.tlc_solvers.size();
    #ifdef VALIDATION
        assert(cyclic_dofs.ndiscrete == get_discrete_dofs().size());
    #endif

    checkpoint("cycpart_dofscheck");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();
        
        py::list log_dofs = py::cast(cyclic_dofs.fulldof_container),
                 log_free = py::cast(cyclic_dofs.get_dofdata(dof_type_t::free)),
                 log_dep = py::cast(cyclic_dofs.get_dofdata(dof_type_t::dep)),
                 log_fixed = py::cast(cyclic_dofs.get_dofdata(dof_type_t::fixed));

        log_routine("cycpart_dofscheck");
        log_in("atoms", log_atoms);
        log_check("dofs_SET4", log_dofs);
        log_check("free_SET4", log_free);
        log_check("dep_SET4", log_dep);
        log_check("fixed_SET4", log_fixed);
        log_finish();
        }
    #endif

    // Create DG and geomunits
    const auto nx = Utils::PyModules::get_nx();
    auto pydg = nx->attr("DiGraph")();
    auto undirdg = nx->attr("Graph")();
    pydg.attr("add_nodes_from")(py::cast(build_sequence));
    undirdg.attr("add_nodes_from")(py::cast(build_sequence));
    geomunit_container = geomunit_container_t(lg_size);
    for (int i = 0; i < build_sequence.size(); ++i) {
        py::list build_subsequence;
        for (int j = 0; j < i + 1; ++j)
            build_subsequence.append(py::cast(build_sequence[j]));
        
        auto subIG_total = ig.attr("subgraph")(build_subsequence);
        auto headnode = build_sequence[i];
        auto headnode_py = py::cast(headnode);
        auto conn_comps = py::list(nx->attr("connected_components")(subIG_total));
        py::list subLG_nodes;
        for (const auto& item : conn_comps) {
            auto nodes = item.cast<py::list>();
            if (py::handle(nodes.attr("__contains__")(headnode_py)).cast<bool>()) {
                subLG_nodes = nodes;
                break;
            }
        }
        if (subLG_nodes.size() == 0) {
            throw std::runtime_error(fmt::format("Suction! headnode={} build_sequence={} build_subsequence={} conn_comps={}", headnode, 
                                        py::handle(py::repr(py::cast(build_sequence))).cast<std::string>(),
                                        py::repr(build_subsequence).cast<std::string>(),
                                        py::repr(conn_comps).cast<std::string>()));
        }
        auto subLG = pylg.attr("subgraph")(subLG_nodes);
        py::list subG_nodes;
        for (const auto& item : subLG.attr("nodes")()) {
            auto lg_node = py::handle(item).cast<int>();
            subG_nodes.attr("__iadd__")(py::cast(assembly_data.ring_atoms[lg_node]));
        }
        auto subG = gr.attr("subgraph")(subG_nodes);
        
        checkpoint("cp_checkDG_A");
        #ifdef KDMOL_LOG
            {
            py::list log_atoms, log_subg_nodes, log_sublg_nodes, log_subg_edges, log_sublg_edges;
            for (const auto& py_idx : gr.attr("nodes")()) {
                auto idx = py::handle(py_idx).cast<py::int_>();
                log_atoms.append(idx);
            }
            log_atoms.attr("sort")();

            for (const auto& py_idx : subG.attr("nodes")()) {
                auto idx = py::handle(py_idx).cast<py::int_>();
                log_subg_nodes.append(idx);
            }
            log_subg_nodes.attr("sort")();

            for (const auto& py_idx : subLG.attr("nodes")()) {
                auto idx = py::handle(py_idx).cast<py::int_>();
                log_sublg_nodes.append(idx);
            }
            log_sublg_nodes.attr("sort")();
            
            for (const auto& py_edge : subG.attr("edges")()) {
                auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
                log_subg_edges.append(py::cast(edge));
            }
            
            for (const auto& py_edge : subLG.attr("edges")()) {
                auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
                log_sublg_edges.append(py::cast(edge));
            }

            log_routine("cp_checkDG_A");
            log_in("atoms", log_atoms);
            log_in("headnode", headnode);
            log_check("subLG_nodes", log_sublg_nodes);
            log_check("subG_nodes", log_subg_nodes);
            log_check("subG_edges_SET2", log_subg_edges);
            log_check("subLG_edges_SET2", log_sublg_edges);
            log_finish();
            }
        #endif

        geomunit_container[i] = geomunit_t(subG, subLG, headnode, this->mol_data);

        if (subLG_nodes.size() > 1) {
            py::int_ incoming = -1;
            std::vector<py::int_> incoming_nodes;
            for (int j = i - 1; j >= 0; --j) {
                auto nb = py::cast(build_sequence[j]);
                if ((py::handle(subLG.attr("has_node")(nb)).cast<bool>()) &&
                    (!py::handle(nx->attr("has_path")(undirdg, headnode_py, nb)).cast<bool>()) &&
                    (py::list(pydg.attr("neighbors")(nb)).size() == 0))
                {
                    incoming_nodes.push_back(nb.cast<py::int_>());
                }
            }

            assertm(incoming_nodes.size() > 0, "Bug found");
            
            for (const auto& incoming_node : incoming_nodes) {
                pydg.attr("add_edge")(incoming_node, headnode_py);
                undirdg.attr("add_edge")(incoming_node, headnode_py);
            }

            for (const auto& incoming_node_py : incoming_nodes) {
                const auto incoming_node = incoming_node_py.cast<int>();
                auto& cur_geomunit = geomunit_container[assembly_data.idx2seq[incoming_node]];

                for (int j = i - 1; j >= 0; --j) {
                    auto nb = build_sequence[j];
                    auto nb_py = py::cast(nb);
                    if ((cur_geomunit.contains_problem(nb)) &&
                        (py::handle(subLG.attr("has_edge")(headnode_py, nb_py)).cast<bool>()))
                    {
                        cur_geomunit.connect_node = nb;
                    }
                }
            }
        }

        /*
        if subLG.number_of_nodes() > 1:
        incoming_nodes = []
        for nb in reversed(self.build_strategy[:i]):
            if subLG.has_node(nb) and not nx.has_path(undirDG, headnode, nb) and len(list(nx.neighbors(self.DG, nb))) == 0:
                incoming_nodes.append(nb)
        assert len(incoming_nodes) > 0
        for node in incoming_nodes:
            self.DG.add_edge(node, headnode)
            undirDG.add_edge(node, headnode)
        for incoming in incoming_nodes:
            for nb in reversed(self.build_strategy[:i]):
                if self.DG.nodes[incoming]['gu'].LG.has_node(nb) and subLG.has_edge(headnode, nb):
                    self.DG.nodes[incoming]['gu'].connect_node = nb
        */
    }

    dg = dg_t(lg_size);
    int link_idx = 0;
    for (const auto& edge : pydg.attr("edges")()) {
        auto [ vA, vB ] = py::handle(edge).cast<std::pair<int, int>>();
        dg[add_edge(vA, vB, dg).first].link_idx = link_idx++;
    }

    checkpoint("cp_checkDGidx");
    #ifdef KDMOL_LOG
        { // Assign 'ord' and 'idx' to each geomunit 
        std::vector<std::pair<int, int>> todo_nodes = {{build_sequence[build_sequence.size() - 1], 0}};
        std::vector<int> idxs = {0};
        while (todo_nodes.size() > 0) {
            auto [ mynode, myord ] = todo_nodes.back();
            todo_nodes.pop_back();
            geomunit_container[assembly_data.idx2seq[mynode]].ord = myord;
            geomunit_container[assembly_data.idx2seq[mynode]].idx = idxs.at(myord);
            idxs.at(myord) += 1;

            boost::property_map <dg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dg);
            boost::graph_traits<dg_t>::in_edge_iterator ei, ei_end;
            for (tie(ei, ei_end) = in_edges(vertex(mynode, dg), dg); ei != ei_end; ++ei) {
                auto othernode = get(index_map, source(*ei, dg));
                todo_nodes.push_back({othernode, myord + 1});
            }

            if (idxs.size() == myord + 1)
                idxs.push_back(0);
        }
        for (const auto& geomunit : geomunit_container) {
            if ((geomunit.ord == -1) || (geomunit.idx == -1))
                throw std::runtime_error(fmt::format("Geomunit(head={}) has ord={}, idx={}", geomunit.headnode, geomunit.ord, geomunit.idx));
        }

        py::list log_dg_idxs, log_atoms;
        for (const auto& geomunit : geomunit_container) {
            py::list newitem;
            newitem.append(py::cast(geomunit.headnode));
            newitem.append(py::cast(geomunit.ord));
            newitem.append(py::cast(geomunit.idx));
            log_dg_idxs.append(newitem);
        }

        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();

        log_routine("cp_checkDGidx");
        log_in("atoms", log_atoms);
        log_check("dg_data_SET1", log_dg_idxs);
        log_finish();
        }
    #endif
    
    

    checkpoint("cp_checkDG");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms, log_pydg_nodes, log_pydg_edges, log_dg_nodes, log_dg_edges;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();

        for (const auto& py_idx : pydg.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_pydg_nodes.append(idx);
        }
        log_pydg_nodes.attr("sort")();
        
        for (const auto& py_edge : pydg.attr("edges")()) {
            auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
            log_pydg_edges.append(py::cast(edge));
        }

        boost::graph_traits <dg_t>::vertex_iterator i, end;
        boost::property_map <dg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dg);

        for (boost::tie(i, end) = vertices(dg); i != end; ++i) {
            log_dg_nodes.append(get(index_map, *i));
        }
        log_dg_nodes.attr("sort")();
        
        boost::graph_traits<dg_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(dg); ei != ei_end; ++ei)
            log_dg_edges.append(py::make_tuple(get(index_map, source(*ei, dg)),
                                               get(index_map, target(*ei, dg))));
        
        log_routine("cp_checkDG");
        log_in("atoms", log_atoms);
        log_check("pynodes", log_pydg_nodes);
        log_check("pyedges_SET2", log_pydg_edges);
        log_check("nodes", log_dg_nodes);
        log_check("edges_SET2", log_dg_edges);
        log_finish();
        }
    #endif

    { // Construct resolvers
        boost::property_map <lg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, lg);
        boost::graph_traits<lg_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(lg); ei != ei_end; ++ei) {
            int dep_node = get(index_map, source(*ei, lg)), ind_node = get(index_map, target(*ei, lg));
            py::list bonds_py = ig.attr("__getitem__")(py::cast(dep_node)).attr("__getitem__")(py::cast(ind_node)).attr("__getitem__")("common_bonds");

            auto edge_descr = edge(dep_node, ind_node, lg).first;
            lg[edge_descr].resolver_start = resolver_container.size();

            int n_resolvers = 0;
            auto common_bonds = bonds_py.cast<std::vector<std::pair<unsigned int, unsigned int>>>();
            for (const auto bond : common_bonds) {
                // Attempt finding independent DOF
                fulldof_t ind_dih;
                std::vector<int> indep_params;
                indep_params.reserve(1);
                {
                    auto [frag_start, frag_end] = cyclic_dofs.get_fragment_dofs(ind_node);
                    for (int i = frag_start; i < frag_end; ++i) {
                        auto dof_bond = cyclic_dofs.dof_container[i];
                        auto dof_bond_rev = std::make_pair(dof_bond.second, dof_bond.first);
                        if (bond == dof_bond || bond == dof_bond_rev)
                            indep_params.push_back(i);
                    }
                }
                if (indep_params.size() > 1)
                    throw std::runtime_error(fmt::format("Bug: several dep DOFs found"));
                // Find the corresponding independent dihedral
                typename resolver_t::dofconstpointer_t ind_val;
                if (indep_params.size() == 1) {
                    const auto indepdof_idx = indep_params[0];
                    ind_dih = cyclic_dofs.fulldof_container[indepdof_idx];
                    ind_val = &cyclic_dofs.dofvalue_container[indepdof_idx];
                } else {
                    ind_dih = problem_container[ind_node].get_valid_dihedral(bond);
                    ind_val = nullptr;
                }
                
                // Find the dependent DOF
                std::vector<int> dep_params;
                dep_params.reserve(1);
                {    
                    auto [frag_start, frag_end] = cyclic_dofs.get_fragment_dofs(dep_node);
                    for (int i = frag_start; i < frag_end; ++i) {
                        auto dof_bond = cyclic_dofs.dof_container[i];
                        auto dof_bond_rev = std::make_pair(dof_bond.second, dof_bond.first);
                        if (bond == dof_bond || bond == dof_bond_rev)
                            dep_params.push_back(i);
                    }
                }
                if (dep_params.size() > 1)
                    throw std::runtime_error(fmt::format("Bug: several dep DOFs found"));
                else if (dep_params.size() == 0)
                    throw std::runtime_error(fmt::format("Bug: unable to find the dep DOF"));
                const auto depdof_idx = dep_params[0];
                typename resolver_t::dofpointer_t dep_val = &cyclic_dofs.dofvalue_container[depdof_idx];;

                auto resolver = resolver_t(ind_dih, ind_val, depdof_idx, dep_val, this->mol_data, assembly_data, cyclic_dofs);
                resolver_container.push_back(resolver);
                n_resolvers++;
            }
            lg[edge_descr].n_resolvers = n_resolvers;
        }
        resolver_container.shrink_to_fit();

        checkpoint("cp_checkResolvers");
        #ifdef KDMOL_LOG
            {
            py::list log_atoms, log_resolvers;
            for (const auto& py_idx : gr.attr("nodes")()) {
                auto idx = py::handle(py_idx).cast<py::int_>();
                log_atoms.append(idx);
            }
            log_atoms.attr("sort")();

            boost::property_map <lg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, lg);
            boost::graph_traits<lg_t>::out_edge_iterator ei, ei_end;
            for (const auto& geomunit : geomunit_container) {
                py::list depdofs;
                for (tie(ei, ei_end) = out_edges(vertex(geomunit.headnode, lg), lg); ei != ei_end; ++ei) {
                    auto edge = *ei;
                    for (int i = lg[edge].resolver_start; i < lg[edge].resolver_start + lg[edge].n_resolvers; ++i){
                        depdofs.append(py::cast(resolver_container[i].dep_dof));
                    }
                }
                py::list newitem;
                newitem.append(py::cast(geomunit.headnode));
                newitem.append(depdofs);
                log_resolvers.append(newitem);
            }

            log_routine("cp_checkResolvers");
            log_in("atoms", log_atoms);
            log_check("depdihs_SET1", log_resolvers);
            log_finish();
            }
        #endif
    }

    { // Check the assembly sequence for corner cases. Use of single XYZ container for extG sometimes causes data corruption
        std::unordered_map<int, std::set<int>> atom_echelons; // Group of atoms, that must be non-overlapping
        std::set<int> cp_atoms;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            cp_atoms.insert(idx);
        }
        
        boost::property_map <dg_t, boost::vertex_index_t>::type dg_index_map = get(boost::vertex_index, dg);
        boost::graph_traits<dg_t>::in_edge_iterator dgei, dgei_end;
        for (auto& geomunit : geomunit_container) {
            std::vector<int> connected_geomunits;
            
            for (tie(dgei, dgei_end) = in_edges(vertex(geomunit.headnode, dg), dg); dgei != dgei_end; ++dgei) {
                auto edge_idx = dg[*dgei].link_idx;
                auto ind_node = get(dg_index_map, source(*dgei, dg)), dep_node = get(dg_index_map, target(*dgei, dg));
                auto ind_headnode = geomunit_container[assembly_data.idx2seq[ind_node]].headnode;
                
                assertm(ind_headnode == ind_node, "ASDHAKLSF");
                connected_geomunits.push_back(ind_headnode);
            }

            for(int i = 0; i < connected_geomunits.size(); i++) {
                for(int j = i+1; j < connected_geomunits.size(); j++) {
                    std::set<int> conflicting_atoms;
                    const auto i_headnode = connected_geomunits[i];
                    const auto j_headnode = connected_geomunits[j];
                    std::set_intersection(atom_echelons.at(i_headnode).cbegin(), atom_echelons.at(i_headnode).cend(),
                                          atom_echelons.at(j_headnode).cbegin(), atom_echelons.at(j_headnode).cend(),
                                          std::inserter(conflicting_atoms, conflicting_atoms.begin()));
                    if (conflicting_atoms.size() > 0) {
                        geomunit_container[assembly_data.idx2seq[i_headnode]].protect_atoms(conflicting_atoms);
                        geomunit_container[assembly_data.idx2seq[j_headnode]].protect_atoms(conflicting_atoms);
                    }
                }
            }
            
            atom_echelons[geomunit.headnode];
            for (const auto [ catom, _ ] : geomunit.inner_container) {
                for (const auto& [ atom, _ ] : this->mol_data->polys.at(catom)) {
                    auto it = cp_atoms.find(atom);
                    if (it != cp_atoms.cend()) {
                        atom_echelons.at(geomunit.headnode).insert(atom);
                    }
                }
            }
        }
        for (auto& geomunit : geomunit_container)
            geomunit.finalize_atoms_protection();
    }

    { // Construct assemblers
        assembler_container = assembler_container_t(num_edges(dg));
        boost::property_map <dg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dg);
        boost::graph_traits<dg_t>::edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = edges(dg); ei != ei_end; ++ei) {
            auto edge_idx = dg[*ei].link_idx;
            auto ind_node = get(index_map, source(*ei, dg)), dep_node = get(index_map, target(*ei, dg));
            auto incoming_node = geomunit_container[assembly_data.idx2seq[ind_node]].connect_node;
            auto common_atoms = py::handle(ig.attr("__getitem__")(incoming_node).attr("__getitem__")(dep_node).attr("__getitem__")("link")).cast<std::vector<int>>();
            auto central_atom = *std::min_element(common_atoms.begin(), common_atoms.end());
            Utils::frame_t dep_frame = problem_container[dep_node].get_valid_frame(central_atom);
            Utils::frame_t ind_frame = problem_container[incoming_node].get_valid_frame(central_atom);
            assembler_container[edge_idx] = assembler_t(ind_frame, dep_frame, this->mol_data, assembly_data);
        }
    }
}

template<class G, class MD, class DF, class DD, class D>
template <class C>
void CyclicPart<G, MD, DF, DD, D>::add_dofset(const std::pair<typename problem_t::dofs_container_t,
                                                          typename problem_t::marked_dofs_t>& newdofs,
                                          const int idx, C& coord, const py::object& gr)
{
    auto& [ dofs, doftypes ] = newdofs;
    cyclic_dofs.fragment_starts[idx] = cyclic_dofs.fulldof_container.size();
    for (const auto& dof : dofs) {
        dof_type_t curtype;
        if (std::find(doftypes.at("free").cbegin(), doftypes.at("free").cend(), dof) != doftypes.at("free").cend())
            curtype = dof_type_t::free;
        else if (std::find(doftypes.at("dep").cbegin(), doftypes.at("dep").cend(), dof) != doftypes.at("dep").cend())
            curtype = dof_type_t::dep;
        else if (std::find(doftypes.at("fixed").cbegin(), doftypes.at("fixed").cend(), dof) != doftypes.at("fixed").cend())
            curtype = dof_type_t::fixed;
        else
            throw std::runtime_error("WFT");
        
        auto dihedral = 0.0;
        if (!this->mol_data->xyzless_init)
            dihedral = coord.get_dihedral(std::get<0>(dof), std::get<1>(dof), std::get<2>(dof), std::get<3>(dof));
        else
            dihedral = this->mol_data->get_provided_dihedral(std::get<0>(dof), std::get<1>(dof), std::get<2>(dof), std::get<3>(dof));
        const bool unrestrained = py::handle(gr.attr("__getitem__")(std::get<1>(dof)).attr("__getitem__")(std::get<2>(dof)).attr("__getitem__")("unrestrained")).cast<bool>();
        const bool exposable = (curtype != dof_type_t::dep) && (unrestrained);
        cyclic_dofs.new_dof(dof, curtype, dihedral, exposable);
    }
}

template<class G, class MD, class DF, class DD, class D>
inline auto CyclicPart<G, MD, DF, DD, D>::get_dofs() const {
    return cyclic_dofs.get_dofs<typename dofs_t::dof_type>();
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::map_dofvalues(const int start, const int size) {
    // This fixes the broken pointers
    for (const auto& geomunit : geomunit_container)
        problem_container[geomunit.headnode].solver->extcoord = &assembly_data.extg_coords;
    for (auto& resolver : resolver_container)
        resolver.coords = &assembly_data.extg_coords;
    for (auto& assembler : assembler_container)
        assembler.coords = &assembly_data.extg_coords;
    for (auto& problem : problem_container)
        problem.set_dofs_pointer(&cyclic_dofs);
    
    if (cyclic_dofs.nexposed > 0) {
        cyclic_dofs.molvalues_start = &this->dofs->dofvalue_container.at(0);
        auto& cp2mol_full = cyclic_dofs.cp2mol_full;
        auto& cp2mol_free = cyclic_dofs.cp2mol_free;
        for (int cp_idx = 0; cp_idx < cyclic_dofs.fulldof_container.size(); ++cp_idx) {
            if (!cyclic_dofs.exposability_container[cp_idx])
                continue;

            auto dof = cyclic_dofs.fulldof_container[cp_idx];
            int mol_idx = -1;
            for (int i = start; i < start + size; ++i) {
                fulldof_t cur_dof = {
                    this->dofs->side_container.at(i).first,
                    this->dofs->bond_container.at(i).first,
                    this->dofs->bond_container.at(i).second,
                    this->dofs->side_container.at(i).second
                };
                if (cur_dof == dof) {
                    if (mol_idx != -1)
                        throw std::runtime_error(fmt::format("Attempting to bind several MOL-level dofs to one dihedral {}", repr(dof)));
                    mol_idx = i;
                }
            }
            
            if (mol_idx == -1)
                throw std::runtime_error(fmt::format("Exposed CP-DOF {} was not assigned at MOL-level", repr(dof)));
            cp2mol_full[cp_idx] = mol_idx;
            if (this->dofs->doftype_container.at(mol_idx) == dofs_t::dof_type::free)
                cp2mol_free[cp_idx] = mol_idx;
        }
    } else {
        cyclic_dofs.molvalues_start = nullptr;
    }

    #ifdef VALIDATION
        { // Initialize validators when cyclic dofs are added in MOL-level container
        std::set<int> main_atoms;
        for (auto [ atom, _ ] : geomunit_container[geomunit_container.size() - 1].inner_container)
            main_atoms.insert(atom);
        conf_validator = validator_t(*(this->mol_data), *(this->dofs), main_atoms);

        std::set<int> ext_atoms;
        for (auto [ atom, _ ] : geomunit_container[geomunit_container.size() - 1].inner_container)
            for (auto& [ nbatom, __ ]: this->mol_data->polys[atom])
                ext_atoms.insert(nbatom);
        ext_validator = validator_t(*(this->mol_data), *(this->dofs), ext_atoms);

        for (auto& problem : problem_container)
            problem.build_validator();
        for (auto& geomunit : geomunit_container)
            geomunit.build_validator(this->mol_data, this->dofs);
        }
    #endif
}

template<class G, class MD, class DF, class DD, class D>
const typename CyclicPart<G, MD, DF, DD, D>::owner_container_t CyclicPart<G, MD, DF, DD, D>::get_discrete_dofs() const {
    owner_container_t res;
    for (const auto& problem : problem_container) {
        auto owner_atoms = problem.get_discrete_dof();
        if (owner_atoms.size() != 0)
            res.push_back(owner_atoms);
    }
    return res;
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::map_discrete_dofvalues(const int start, const int size) {
    checkpoint("cycpart_discrete_dofs");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms, log_ddof;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();

        for (int i = 0; i < problem_container.size(); ++i) {
            const auto& problem = problem_container[i];
            auto owner_atoms = problem.get_discrete_dof();
            if (owner_atoms.size() == 0)
                continue;

            auto py_atom_set = py::handle(py::cast(owner_atoms)).cast<py::list>();
            py_atom_set.attr("sort")();

            py::list newitem;
            newitem.append(i);
            newitem.append(py_atom_set);
            log_ddof.append(newitem);
        }

        log_routine("cycpart_discrete_dofs");
        log_in("atoms", log_atoms);
        log_check("discretedofs_SET1", log_ddof);
        log_finish();
        }
    #endif

    if (cyclic_dofs.ndiscrete > 0) {
        cyclic_dofs.molvalues_discrete_start = &this->discrete_dofs->dofvalue_container.at(0);
        auto& cp2mol = cyclic_dofs.cp2mol_discrete;
        for (int cp_idx = 0; cp_idx < problem_container.size(); ++cp_idx) {
            auto owner_atoms = problem_container[cp_idx].get_discrete_dof();
            if (owner_atoms.size() == 0)
                continue;

            int mol_idx = -1;
            for (int i = start; i < start + size; ++i) {
                auto& cur_atoms = this->discrete_dofs->owner_container.at(i);
                if (cur_atoms == owner_atoms) {
                    if (mol_idx != -1)
                        throw std::runtime_error(fmt::format("Attempting to bind several MOL-level discrete dofs to one set of atoms {}", repr(owner_atoms)));
                    mol_idx = i;
                }
            }
            
            if (mol_idx == -1)
                throw std::runtime_error(fmt::format("Discrete CP-DOF was not assigned at MOL-level. Atoms={}", repr(owner_atoms)));
            cp2mol[cp_idx] = mol_idx;
        }
    } else {
        cyclic_dofs.molvalues_discrete_start = nullptr;
    }

    for (int prob_idx = 0; prob_idx < problem_container.size(); ++prob_idx) {
        if (problem_container[prob_idx].solvertype != problem_t::solvertype_t::Identity) {
            assert(problem_container[prob_idx].solvertype != problem_t::solvertype_t::Undefined);
            problem_container[prob_idx].set_ddof_pointer(prob_idx);
        }
    }
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::configure_assembly(const main_graph_t& molgr, const dfg_t& dfg, const int given_index) {
    if (this->index == -1)
        this->index = given_index;
    else if (this->index != given_index)
        throw std::runtime_error(fmt::format("Mismatch of indices given by MOL: this->index={} but got {}", repr(this->index), repr(given_index)));

    // Initialize coordinates
    for (auto [ atom, _ ] : geomunit_container[geomunit_container.size() - 1].inner_container)
        for (auto& [ nbatom, __ ]: this->mol_data->polys[atom])
            assembly_data.extg_coords.set_atom(nbatom, {0.0, 0.0, 0.0});

    auto myvertex = vertex(index, dfg);
    int out_count = 0;
    typename boost::graph_traits<dfg_t>::out_edge_iterator out, out_end;
    typename dfg_t::edge_descriptor out_edge;
    for (boost::tie(out, out_end) = out_edges(myvertex, dfg); out != out_end; ++out) {
        out_count++;
        out_edge = *out;
    }

    if (out_count == 1) {
        recalc_frame = true;
        int param_idx = dfg[out_edge].param_idx;
        
        if (this->mol_data->belongs_to_fragment(this->dofs->bond_container[param_idx].second, this->index)) {
            at0_idx = this->dofs->bond_container[param_idx].first;
            at1_idx = this->dofs->bond_container[param_idx].second;
            at2_idx = this->dofs->side_container[param_idx].second;
        } else {
            at0_idx = this->dofs->bond_container[param_idx].second;
            at1_idx = this->dofs->bond_container[param_idx].first;
            at2_idx = this->dofs->side_container[param_idx].first;
        }
    } else if (out_count == 0) {
        recalc_frame = false;
    } else {
        throw std::runtime_error("Experiencing problems with the tree");
    }

    #ifdef KDMOL_LOG
        py::int_ recalc_int;
        if (recalc_frame)
            recalc_int = 1;
        else
            recalc_int = -1;

        // log_routine("cpconfigureB");
        // log_in("atoms", py_atoms);
        // log_check("recalc_frame", as_list(recalc_int.cast<int>()));
        // log_check("at0_idx", as_list(at0_idx));
        // log_check("at1_idx", as_list(at1_idx));
        // log_check("at2_idx", as_list(at2_idx));
        // log_finish();
    #endif

    { // Assign bonds that belong to this cyclic part
        for (const auto& py_edge : gr.attr("edges")()) { // Inner bonds of cyclic part
            auto edge = py::handle(py_edge).cast<std::pair<int, int>>();
            this->mol_data->assign_parent_fragment(edge, this->index);
        }

        // Bonds that do not connect to other fragments (like C-H)
        const auto& frag_coords = this->mol_data->fg_container[this->index];
        for (const auto& py_at_idx : gr.attr("nodes")()) {
            auto at_idx = py::handle(py_at_idx).cast<int>();
            for (const auto& [nb_idx, _] : this->mol_data->edge_indexer.atom_data(at_idx))
                if (frag_coords.find(nb_idx) != frag_coords.end())
                    this->mol_data->assign_parent_fragment({at_idx, nb_idx}, this->index);
        }

        // Outcoming Bonds in dfg that connect to other fragments
        typename boost::property_map <dfg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dfg);
        typename boost::graph_traits<dfg_t>::in_edge_iterator in, in_end;
        for (boost::tie(in, in_end) = in_edges(myvertex, dfg); in != in_end; ++in) {
            auto in_edge = *in;
            auto in_node_idx = get(index_map, source(in_edge, dfg));
            this->mol_data->assign_parent_fragment(this->dofs->bond_container[dfg[in_edge].param_idx], this->index);
        }
    }
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::index_owned_dofs() {
    // Cache idxs of free dihedrals for refinement
    own_dof_idxs = active_idxs_t();
    for (const auto [_, mol_idx] : cyclic_dofs.cp2mol_free)
        own_dof_idxs.push_back(this->dofs->free_indexer[mol_idx]);
    own_dof_idxs.shrink_to_fit();

    // Cache DDOF idxs for refinement
    own_discrete_idxs = active_idxs_t();
    for (const auto [_, mol_idx] : cyclic_dofs.cp2mol_discrete)
        own_discrete_idxs.push_back(mol_idx);
    own_discrete_idxs.shrink_to_fit();

    for (int prob_idx = 0; prob_idx < problem_container.size(); ++prob_idx)
        problem_container[prob_idx].index_owned_dofs(prob_idx);
    
     for (auto& geomunit : geomunit_container) {
        auto [ own_dof_idxs, own_discrete_idxs ] = problem_container[geomunit.headnode].get_owned_dofs();
        if (own_dof_idxs.size() + own_discrete_idxs.size() == 0) {
            own_dof_idxs = active_idxs_t();
            own_discrete_idxs = active_idxs_t();
            for (const auto problem_idx : geomunit.problem_idxs) {
                const auto [ new_dof_idxs, new_discrete_idxs ] = problem_container[problem_idx].get_owned_dofs();
                std::copy(new_dof_idxs.cbegin(), new_dof_idxs.cend(), std::back_inserter(own_dof_idxs));
                std::copy(new_discrete_idxs.cbegin(), new_discrete_idxs.cend(), std::back_inserter(own_discrete_idxs));
            }
        }
        geomunit.index_owned_dofs(own_dof_idxs, own_discrete_idxs);
     }
}

template<class G, class MD, class DF, class DD, class D>
aps_return_t CyclicPart<G, MD, DF, DD, D>::apply_ps() {
    // Copy values of free DOFs
    for (const auto [cp_idx, mol_idx] : cyclic_dofs.cp2mol_free) {
        cyclic_dofs.dofvalue_container.at(cp_idx) = cyclic_dofs.molvalues_start[mol_idx];

        #ifdef VALIDATION
            auto my_dof = cyclic_dofs.fulldof_container[cp_idx];
            auto mol_dof = std::make_tuple(
                this->dofs->side_container.at(mol_idx).first,
                this->dofs->bond_container.at(mol_idx).first,
                this->dofs->bond_container.at(mol_idx).second,
                this->dofs->side_container.at(mol_idx).second
            );
            if (my_dof != mol_dof)
                throw std::runtime_error(fmt::format("Mismatch of dofs my the map: cp={} mol={}", repr(my_dof), repr(mol_dof)));
            if (cyclic_dofs.doftype_container.at(cp_idx) != dof_type_t::free)
                throw std::runtime_error(fmt::format("For some reason non-free DOF {} was assigned by apply_ps()", repr(my_dof)));
            if (this->dofs->doftype_container.at(mol_idx) != dofs_t::dof_type::free) // This would be a really weird case
                throw std::runtime_error(fmt::format("For some reason non-free DOF {} was assigned by apply_ps()", repr(mol_dof)));
        #endif
    }

    checkpoint("cpaps_checkInputDOF");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();
        
        py::list input_data;
        for (int i = 0; i < cyclic_dofs.size(); ++i) {
            if (cyclic_dofs.doftype_container.at(i) == dof_type_t::free) {
                input_data.append(py::make_tuple(std::get<0>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<1>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<2>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<3>(cyclic_dofs.fulldof_container.at(i)),
                                                 cyclic_dofs.dofvalue_container[i]));
            }
        }

        log_routine("cpaps_checkInputDOF");
        log_in("call_idx", Logger::aps_count);
        log_in("atoms", log_atoms);
        log_check("param_data_SET4", input_data);
        log_finish();
        }
    #endif

    boost::property_map <lg_t, boost::vertex_index_t>::type lg_index_map = get(boost::vertex_index, lg);
    boost::graph_traits<lg_t>::out_edge_iterator lgei, lgei_end;
    boost::property_map <dg_t, boost::vertex_index_t>::type dg_index_map = get(boost::vertex_index, dg);
    boost::graph_traits<dg_t>::in_edge_iterator dgei, dgei_end;
    for (auto& geomunit : geomunit_container) {
        for (tie(lgei, lgei_end) = out_edges(vertex(geomunit.headnode, lg), lg); lgei != lgei_end; ++lgei) {
            auto edge = *lgei;
            for (int i = lg[edge].resolver_start; i < lg[edge].resolver_start + lg[edge].n_resolvers; ++i)
                resolver_container[i].resolve_dependence();
        }

        auto problem_result = problem_container[geomunit.headnode].apply_ps();
        if (problem_result != aps_return_t::success)
            return problem_result;
        
        auto dep_coords = problem_container[geomunit.headnode].get_coords();
        for (tie(dgei, dgei_end) = in_edges(vertex(geomunit.headnode, dg), dg); dgei != dgei_end; ++dgei) {
            auto edge_idx = dg[*dgei].link_idx;
            auto ind_node = get(dg_index_map, source(*dgei, dg)), dep_node = get(dg_index_map, target(*dgei, dg));
            auto ind_gu = geomunit_container[assembly_data.idx2seq[ind_node]];
            auto incoming_node = geomunit_container[assembly_data.idx2seq[ind_node]].connect_node;

            if (incoming_node == -1)
                throw std::runtime_error("Incoming node was not assigned");
            assembler_container[edge_idx].connect_cycles(dep_coords, ind_gu.inner_container
                #ifdef KDMOL_LOG
                    , geomunit.ord
                #endif
                );
        }

        for (auto& [atom, xyz] : dep_coords)
            assembly_data.extg_coords[atom] = xyz;

        auto gu_result = geomunit.apply_ps(assembly_data.extg_coords, this->mol_data, this->dofs);
        if (gu_result != aps_return_t::success)
            return gu_result;

        checkpoint("cpaps_interm");
        #ifdef KDMOL_LOG
            {
            py::list log_atoms, log_coords;
            for (const auto& py_idx : gr.attr("nodes")()) {
                auto idx = py::handle(py_idx).cast<py::int_>();
                log_atoms.append(idx);
            }
            log_atoms.attr("sort")();

            for (auto& [ atom, _ ] : geomunit.inner_container) {
                py::list newitem;
                newitem.append(atom);
                newitem.append(assembly_data.extg_coords[atom](0));
                newitem.append(assembly_data.extg_coords[atom](1));
                newitem.append(assembly_data.extg_coords[atom](2));
                log_coords.append(newitem);
            }

            log_routine("cpaps_interm");
            log_in("call_idx", Logger::aps_count);
            log_in("atoms", log_atoms);
            log_in("head", geomunit.headnode);
            log_check("coords_SET1", log_coords);
            log_finish();
            }
        #endif
    }

    checkpoint("cpaps_checkInner");
    #ifdef KDMOL_LOG
        {
        py::list log_atoms, input_data, log_coords;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<py::int_>();
            log_atoms.append(idx);
        }
        log_atoms.attr("sort")();

        for (int i = 0; i < cyclic_dofs.size(); ++i) {
            if (cyclic_dofs.doftype_container.at(i) == dof_type_t::free)
                input_data.append(py::make_tuple(std::get<0>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<1>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<2>(cyclic_dofs.fulldof_container.at(i)),
                                                 std::get<3>(cyclic_dofs.fulldof_container.at(i)),
                                                 cyclic_dofs.dofvalue_container[i]));
        }

        for (auto& [ atom, _ ] : geomunit_container[geomunit_container.size() - 1].inner_container) {
            py::list newitem;
            newitem.append(atom);
            newitem.append(assembly_data.extg_coords[atom](0));
            newitem.append(assembly_data.extg_coords[atom](1));
            newitem.append(assembly_data.extg_coords[atom](2));
            log_coords.append(newitem);
        }

        log_routine("cpaps_checkInner");
        log_in("call_idx", Logger::aps_count);
        log_in("atoms", log_atoms);
        log_check("param_data_SET4", input_data);
        log_check("coords_SET1", log_coords);
        log_finish();
        }
    #endif

    #ifdef VALIDATION
        if (!conf_validator.validate(assembly_data.extg_coords, *(this->mol_data), *(this->dofs))) {
            this->mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
            return aps_return_t::validationfail;
            // throw Utils::GeomFailure("Validation of cyclic part conformation failed");
        }
        if (!ext_validator.validate(assembly_data.extg_coords, *(this->mol_data), *(this->dofs))) {
            this->mol_data->refiners.set_active(own_dof_idxs, own_discrete_idxs);
            return aps_return_t::validationfail;
            // throw Utils::GeomFailure("Validation of extended cyclic part failed");
        }
    #endif
    return aps_return_t::success;
}

template<class G, class MD, class DF, class DD, class D>
void CyclicPart<G, MD, DF, DD, D>::update_configuration(const dfg_t& dfg) {
    se_matrix_t startframe, inv_startframe;
    if (this->recalc_frame) {
        startframe = this->mol_data->polys[at1_idx].get_frame(assembly_data.extg_coords,
                                            {at0_idx, at2_idx}, at1_idx);
        for (int i = 0; i < 3; ++i) {
            startframe(i, 0) = -startframe(i, 0);
            startframe(i, 2) = -startframe(i, 2);
        }
    } else {
        startframe.assign(boost::numeric::ublas::identity_matrix<typename se_matrix_t::value_type>(startframe.size1()));
    }
    Utils::invert_matrix(startframe, inv_startframe);

    #ifdef KDMOL_LOG
        py::list loc_coords;
    #endif
    for (const auto& atom_idx : this->mol_data->fg_container[this->index]) {
        Utils::uvector4d_t tempcoord;
        for(int i = 0; i < 3; ++i)
            tempcoord(i) = assembly_data.extg_coords.at(atom_idx)(i);
        tempcoord(3) = 1.0;
        noalias(this->mol_data->local_xyz.get_4d(atom_idx)) = prod(inv_startframe, tempcoord);

        #ifdef KDMOL_LOG
            py::list newitem;
            newitem.append(atom_idx);
            newitem.append(Utils::vector_to_list(this->mol_data->local_xyz.get_4d(atom_idx)));
            loc_coords.append(newitem);
        #endif
    }

    #ifdef VALIDATION
        if (!conf_validator.validate(this->mol_data->local_xyz, *(this->mol_data), *(this->dofs))) {
            // return aps_return_t::validationfail;
            throw Utils::GeomFailure("Validation of the local coordinates failed");
        }
    #endif

    auto myvertex = vertex(index, dfg);
    typename boost::property_map <dfg_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, dfg);
    typename boost::graph_traits<dfg_t>::in_edge_iterator in, in_end;
    #ifdef KDMOL_LOG
        py::list transitions;
    #endif
    for (boost::tie(in, in_end) = in_edges(myvertex, dfg); in != in_end; ++in) {
        auto in_edge = *in;
        auto in_node_idx = get(index_map, source(in_edge, dfg));
        int param_idx = dfg[in_edge].param_idx;
        se_matrix_t newframe;
        if (this->mol_data->belongs_to_fragment(this->dofs->bond_container[param_idx].second, this->index)) {
            newframe = this->mol_data->polys[this->dofs->bond_container[param_idx].second].get_frame(assembly_data.extg_coords,
                                            {this->dofs->bond_container[param_idx].first,
                                             this->dofs->side_container[param_idx].second},
                                             this->dofs->bond_container[param_idx].first);
        } else {
            newframe = this->mol_data->polys[this->dofs->bond_container[param_idx].first].get_frame(assembly_data.extg_coords,
                                            {this->dofs->bond_container[param_idx].second,
                                             this->dofs->side_container[param_idx].first},
                                             this->dofs->bond_container[param_idx].second);
        }
        noalias(this->mol_data->transition_map[param_idx]) = prod(inv_startframe, newframe);

        #ifdef KDMOL_LOG
            py::list newitem;
            newitem.append(in_node_idx);
            newitem.append(Utils::matrix_to_list(this->mol_data->transition_map[param_idx]));
            transitions.append(newitem);
        #endif
    }

    checkpoint("cpupdate");
    #ifdef KDMOL_LOG
        py::list py_atoms, py_inner_data, py_outer_data;
        for (const auto& py_idx : gr.attr("nodes")()) {
            auto idx = py::handle(py_idx).cast<int>();
            py_atoms.append(idx);
        }
        py_atoms.attr("sort")();

        log_routine("cpupdate");
        log_in("atoms", py_atoms);
        log_in("call_idx", Logger::aps_count);
        log_check("at0_idx", as_list(at0_idx));
        log_check("at1_idx", as_list(at1_idx));
        log_check("at2_idx", as_list(at2_idx));
        log_check("startframe", startframe);
        log_check("loccoords_SET1", loc_coords);
        log_check("transitions_SET1", transitions);
        log_finish();
    #endif
}

template<class G, class MD, class DF, class DD, class D>
int CyclicPart<G, MD, DF, DD, D>::get_num_tlcsolvers() const noexcept {
    return assembly_data.tlc_solvers.size();
}

template<class G, class MD, class DF, class DD, class D>
int CyclicPart<G, MD, DF, DD, D>::get_num_idsolvers() const noexcept {
    return assembly_data.id_solvers.size();
}

#endif // CYCLICPART_H
