#ifndef MOLECULE_H
#define MOLECULE_H

#include "polyhedron.h"
#include "cyclicpart.h"
#include "fragment.h"
#include "refiners.h"


#if defined(VALIDATION) || defined(WEAK_VALIDATION)
    #define MOL_VALIDATE
#endif
#if defined(OVERLAP_DETECTION) || defined(OVERLAP_DETECTION_FINAL)
    #define MOL_OVERLAP_DETECTION
#endif

template <class P, class G, class DF, class DD>
struct MolData {
    using coord_container_t = Utils::BoostXyzContainer<3>;
    coord_container_t xyz;
    using localcoord_container_t = Utils::Boost4DAdapter<Utils::BoostXyzContainer<4>>;
    localcoord_container_t local_xyz;

    using main_graph_t = G;
    main_graph_t molgr;

    template <class T> using base_container_t = std::vector<T>;

    using atom_symbol_t = std::string;
    using atom_symbols_t = base_container_t<atom_symbol_t>;
    atom_symbols_t atom_symbols;
    
    using poly_t = P;
    base_container_t<poly_t> polys;

    using bondtype_t = int;
    base_container_t<bondtype_t> bondtypes;
    using bondlength_t = double;
    base_container_t<bondlength_t> bondlengths;

    using fragatoms_t = std::set<int>;
    using fg_container_t = base_container_t<fragatoms_t>;
    fg_container_t fg_container;
    using frame_container_t = base_container_t<se_matrix_t>;
    frame_container_t frame_container;
    using frag_neighbors_t = std::vector<std::vector<int>>;
    frag_neighbors_t frag_neighbors;

    using indexer_t = Utils::EdgeIndexer;
    indexer_t edge_indexer;
    using atom_in_frag_t = base_container_t<int>;
    atom_in_frag_t atom_in_frag;

    using transition_map_t = std::map<int, se_matrix_t>;
    transition_map_t transition_map;

    using pcg_t = pcg32; // For random numbers
    pcg_t pcg_rng;

    using dofs_t = DF;
    using discrete_dofs_t = DD;
    using refiners_t = Refiners<dofs_t, discrete_dofs_t>;
    refiners_t refiners;

    bool xyzless_init;
    std::map<std::tuple<int, int, int, int>, double> provided_dihedrals;
    static constexpr double DEFAULT_DIHEDRAL = 0.0;
    
    MolData() { }
    
    MolData(const int& natoms, const int& nbonds, dofs_t& dofs, discrete_dofs_t& discrete_dofs) :
        xyz(natoms),
        local_xyz(natoms),
        atom_symbols(natoms),
        polys(natoms),
        bondtypes(nbonds),
        bondlengths(nbonds),
        refiners(refiners_t(&dofs, &discrete_dofs)) // Pass here a pointer to 'pcg_rng'??
    { }

    inline int size() const noexcept { return atom_symbols.size(); }

    inline int number_of_neighbors(const indexer_t::atom_t& atom_idx) const noexcept {
        return edge_indexer.atom_data(atom_idx).size();
    }

    inline int get_some_neighbor(const indexer_t::atom_t& atom_idx, const indexer_t::atom_t& banned_idx) const {
        #ifdef KDMOL_DEBUG
            indexer_t::atom_t min_idx = -1;
            for (const auto& [nb_idx, _] : edge_indexer.atom_data(atom_idx))
                if ((nb_idx != banned_idx) && ((min_idx == -1) || (nb_idx < min_idx)))
                    min_idx = nb_idx;
            if (min_idx == -1)
                throw std::runtime_error(fmt::format("Error while searching for appropriate neighbor of {} (banned = {})", atom_idx, banned_idx));
            return min_idx; // Return the smallest index for easier debugging against reference code
        #else
            for (const auto [nb_idx, _] : edge_indexer.atom_data(atom_idx))
                if (nb_idx != banned_idx)
                    return nb_idx;
            throw std::runtime_error(fmt::format("Error while searching for appropriate neighbor of {} (banned = {})", atom_idx, banned_idx));
        #endif
    }

    using bond_t = std::pair<indexer_t::atom_t, indexer_t::atom_t>;

    inline const auto& bond_data(const bond_t& bond) const noexcept {
        return edge_indexer.atom_data(bond.first).at(bond.second);
    }

    inline auto& set_bond_data(const bond_t& bond) noexcept {
        return edge_indexer.set_atom_data(bond.first).at(bond.second);
    }

    inline void assign_parent_fragment(const bond_t& edge, const int& frag_idx) {
        #ifdef KDMOL_DEBUG
            if ((bond_data(edge).parent_frag != -1) && (bond_data(edge).parent_frag != frag_idx))
                throw std::runtime_error(fmt::format("Broken index for bond {}-{}", edge.first, edge.second));
            if ((bond_data({edge.second, edge.first}).parent_frag != -1) && (bond_data({edge.second, edge.first}).parent_frag != frag_idx))
                throw std::runtime_error(fmt::format("Broken index for bond {}-{}", edge.second, edge.first));
        #endif
        set_bond_data(edge).parent_frag = frag_idx;
        set_bond_data({edge.second, edge.first}).parent_frag = frag_idx;
    }

    inline const bondtype_t& get_bondtype(const bond_t& bond) const noexcept {
        return bondtypes[bond_data(bond).index];
    }

    inline bool belongs_to_fragment(const int& atom_idx, const int& frag_idx) const noexcept {
        const auto& frag_atoms = fg_container[frag_idx];
        return frag_atoms.find(atom_idx) != frag_atoms.cend();
    }

    void set_provided_dihedrals (py::dict user_dihedrals) {
        for (auto item : user_dihedrals) {
            auto key = item.first.cast<std::tuple<int, int, int, int>>();
            auto value = item.second.cast<double>();
            provided_dihedrals[std::make_tuple(std::get<0>(key) - 1,
                                               std::get<1>(key) - 1,
                                               std::get<2>(key) - 1,
                                               std::get<3>(key) - 1)] = value;
        }
    }

    double get_provided_dihedral (const int a, const int b, const int c, const int d) const {
        // Create a tuple with the input values
        std::tuple<int, int, int, int> key{a, b, c, d};

        // Check if the tuple is a key in the provided dihedrals map
        auto it = provided_dihedrals.find(key);
        if (it != provided_dihedrals.end()) {
            // Return the corresponding value if the tuple exists in the map
            return it->second * DEG2RAD;
        } else {
            if (get_bondtype({b, c}) != 1)
                fmt::print("[WARNING] Dihedral {}-{}-{}-{} was not provided. Taking the default value = {}\n", a+1, b+1, c+1, d+1, DEFAULT_DIHEDRAL);
            return DEFAULT_DIHEDRAL * DEG2RAD;
        }
    }
};

struct Dofs { // Structure is as follow --NONCYCLIC--|--CP1--|--CP2--, etc.
    using value_t = double;
    using index_t = int;
    using bondatoms_t = std::pair<index_t, index_t>;
    using sideatoms_t = std::pair<index_t, index_t>;
    template <class T> using container_t = std::vector<T>; // Container must be ordered
    
    using bond_container_t = container_t<bondatoms_t>;
    bond_container_t bond_container;
    using side_container_t = container_t<sideatoms_t>;
    side_container_t side_container;
    using dofvalue_container_t = container_t<value_t>;
    dofvalue_container_t dofvalue_container;

    enum dof_type : bool { fixed = false, free = true };
    using doftype_container_t = container_t<dof_type>;
    doftype_container_t doftype_container;

    using free_indexer_t = std::map<int, int>; // Idx in dofvalue_container -> Idx in dof_rawpointer_t
    free_indexer_t free_indexer;
    
    using dof_rawpointer_t = value_t*;
    dof_rawpointer_t dof_rawpointer;
    using dof_nparray_t = py::array_t<value_t>;
    dof_nparray_t dof_nparray;
    using dof_vector_t = std::vector<value_t>;
    dof_vector_t dof_vector;

    using dof_indexer_t = std::map<bondatoms_t, int>;
    dof_indexer_t dof_indexer;

    int noncyclic_size;
    using cpints_container_t = container_t<int>;
    cpints_container_t cpdof_starts, cpdof_sizes;

    using requested_free_t = container_t<bondatoms_t>;
    requested_free_t requested_free;
    
    bool require_best_sequence;

    using custom_sampling_limits_t = std::unordered_map<int, std::pair<double, double>>;
    using custom_limits_data_t = std::map< Utils::fulldof_t, std::pair<double, double> >;

    Dofs() : dof_rawpointer(nullptr), require_best_sequence(false) { }

    void new_dof(const bondatoms_t& bondatoms, const sideatoms_t& sideatoms, const dof_type doftype, const value_t value) {
        bond_container.push_back(bondatoms);
        side_container.push_back(sideatoms);
        doftype_container.push_back(doftype);
        dofvalue_container.push_back(value);

        auto found_iter = dof_indexer.find(bondatoms);
        if (found_iter != dof_indexer.end())
            throw std::runtime_error(fmt::format("Double addition of DOF at MOL-level {}", repr(bondatoms)));
        dof_indexer[bondatoms] = bond_container.size() - 1;

        bondatoms_t bondatoms_inv = {bondatoms.second, bondatoms.first};
        found_iter = dof_indexer.find(bondatoms_inv);
        if (found_iter != dof_indexer.end())
            throw std::runtime_error(fmt::format("Double addition of DOF at MOL-level {}", repr(bondatoms_inv)));
        dof_indexer[bondatoms_inv] = bond_container.size() - 1;
    }

    void finalize() {
        bond_container.shrink_to_fit();
        side_container.shrink_to_fit();
        doftype_container.shrink_to_fit();
        dofvalue_container.shrink_to_fit();

        nfixed = 0;
        nfree = 0;
        for (const auto& item : doftype_container) {
            if (item == fixed)
                nfixed++;
            else
                nfree++;
        }
        // Generate numpy array in future
    }

    inline int get_nfree() const noexcept { return nfree; }
    
    inline int get_nfixed() const noexcept { return nfixed; }
    
    inline int size() const noexcept { return bond_container.size(); }

    inline int param_from_bond(const bondatoms_t& bond) const noexcept { return dof_indexer.at(bond); }

    bool contain(const Utils::fulldof_t dof) {
        bondatoms_t bond = {std::get<1>(dof), std::get<2>(dof)};
        auto found_iter = dof_indexer.find(bond);
        if (found_iter != dof_indexer.cend()) {
            const auto idx = dof_indexer[bond];
            sideatoms_t sides = {std::get<0>(dof), std::get<3>(dof)};
            if (sides != side_container[idx])
                throw std::runtime_error(fmt::format("Mismatch of sides my={} added={} for dof={}", repr(side_container[idx]), repr(sides), repr(dof)));
        }
        return found_iter != dof_indexer.cend();
    }
    
    template <class M>
    inline double get_tmat(const int& param_idx, M& res, const int& mult) const noexcept {
        res.assign(boost::numeric::ublas::identity_matrix<typename M::value_type>(res.size1()));
        double ang = -dofvalue_container[param_idx] * mult;
        const double cos_ang = cos(ang);
        const double sin_ang = sin(ang);
        res(1, 1) = cos_ang; res(1, 2) = -sin_ang;
        res(2, 1) = sin_ang; res(2, 2) = cos_ang;
        return -ang;
    }

    bool pointer_is_good() const noexcept { return ((size() == 0) || (dof_rawpointer != nullptr)); }

    private:
        int nfixed, nfree;
};

struct DiscreteDofs { // Structure is as follow --CP1--|--CP2--|--CP3--, etc. TODO: implement noncyclic ddofs
    using value_t = int;
    using index_t = int;
    template <class T> using container_t = std::vector<T>; // Container must be ordered
    
    using owner_atoms_t = std::set<int>;
    using owner_container_t = container_t<owner_atoms_t>;
    owner_container_t owner_container;

    using dofvalue_container_t = container_t<value_t>;
    dofvalue_container_t dofvalue_container;
    
    using dof_rawpointer_t = value_t*;
    dof_rawpointer_t dof_rawpointer;
    using dof_nparray_t = py::array_t<value_t>;
    dof_nparray_t dof_nparray;
    using dof_vector_t = std::vector<value_t>;
    dof_vector_t dof_vector;

    // TODO Implement container for maximal allowed values of each dof

    // using dof_indexer_t = std::map<bondatoms_t, int>;
    // struct set_compare {
    //     bool operator()(const std::set<int>& lhs, const std::set<int>& rhs) const {
    //         return lhs < rhs;
    //     }
    // };
    // std::map<std::set<int>, std::string, set_compare> my_map;

    using cpints_container_t = container_t<index_t>;
    cpints_container_t cpdof_starts, cpdof_sizes;

    DiscreteDofs() : dof_rawpointer(nullptr) { }

    void new_dof(const owner_atoms_t& owner_atoms) {
        owner_container.push_back(owner_atoms);
        dofvalue_container.push_back(-1);
    }

    inline int size() const noexcept { return dofvalue_container.size(); }

    void finalize() {
        owner_container.shrink_to_fit();
        dofvalue_container.shrink_to_fit();
    }

    bool pointer_is_good() const noexcept { return ((size() == 0) || (dof_rawpointer != nullptr)); }
};

struct DfgVertexProps
{ };

template <class paramidx_t>
struct DfgEdgeProps
{
    paramidx_t param_idx;
};

class Molecule {
    // OPT: Don't use Python-graphs
    using main_graph_t = py::object; // These are considered
    using subgraph_t = main_graph_t; // to be the same in GenericFragment and its childred

    using dofs_t = Dofs;
    using discrete_dofs_t = DiscreteDofs;

    using local_xyz_t = boost::numeric::ublas::fixed_vector<double, 3>;
    using poly_t = Polyhedron<local_xyz_t>;
    using mol_data_t = MolData<poly_t, main_graph_t, dofs_t, discrete_dofs_t>;

    // Directed Graph of Fragments = dfg
    using dfgvertex_t = DfgVertexProps;
    using dfgedge_t = DfgEdgeProps<dofs_t::index_t>;
    using dfg_t = boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, dfgvertex_t, dfgedge_t>;

    using cyclic_part_t = CyclicPart<main_graph_t, mol_data_t, dofs_t, discrete_dofs_t, dfg_t>;
    using cp_container_t = std::vector<cyclic_part_t>;

    using fragment_t = Fragment<main_graph_t, mol_data_t, dofs_t, discrete_dofs_t, dfg_t>;
    using fragment_container_t = std::vector<fragment_t>;
    using generic_fragment_t = GenericFragment<main_graph_t, mol_data_t, dofs_t, discrete_dofs_t, dfg_t>;
    using fragref_container_t = std::vector<generic_fragment_t*>;

    using todo_queue_t = std::vector<int>;


    #ifdef MOL_VALIDATE
        using validator_t = Utils::GeometryValidator;
    #endif

    #ifdef MOL_OVERLAP_DETECTION
        using overlap_detector_t = Utils::OverlapDetector;
    #endif

    
    public:
        using active_idxs_t = typename mol_data_t::refiners_t::active_idxs_t;
        using owner_container_t = discrete_dofs_t::owner_container_t;
        using requested_free_t = dofs_t::requested_free_t;
        using custom_sampling_limits_t = dofs_t::custom_sampling_limits_t;
        using custom_limits_data_t = dofs_t::custom_limits_data_t;

        Molecule() { }
        template <class C, bool xyzless> void constructor(C input_object, py::module_& nx, py::module_& pyutils);

        aps_return_t apply_ps();
        py::int_ apply_ps_py() { return py::cast(static_cast<int>(apply_ps())); }
        void reconfigure(const py::dict& upd_data);
        void set_requested_dofs_py(py::list requested_dofs) { dofs.requested_free = requested_dofs.cast<requested_free_t>(); }
        void set_requested_dofs(requested_free_t requested_dofs) { dofs.requested_free = requested_dofs; }
        requested_free_t get_requested_dofs() { return dofs.requested_free; }

        std::pair<std::vector<Utils::fulldof_t>, double*> get_ps();
        py::tuple get_ps_py();
        std::pair<discrete_dofs_t::owner_container_t, int*> get_discrete_ps();
        py::tuple get_discrete_ps_py();
        inline const typename mol_data_t::coord_container_t& get_xyz() const noexcept { return mol_data.xyz; }
        py::array_t<double> get_xyz_py() const;
        inline const typename mol_data_t::atom_symbols_t& get_symbols() const noexcept { return mol_data.atom_symbols; }
        py::list get_symbols_py() const { return py::cast(mol_data.atom_symbols); }
        int get_num_flexible_rings() const;
        int get_num_rigid_rings() const;
        inline int get_size() const noexcept { return mol_data.size(); }
        inline void init_rng(const int seed) { mol_data.pcg_rng = typename mol_data_t::pcg_t(seed); }

        inline const typename mol_data_t::coord_container_t& const_coord_access() const noexcept { return mol_data.xyz; }
        inline typename mol_data_t::coord_container_t& coord_access() noexcept { return mol_data.xyz; }
        inline const typename mol_data_t::atom_symbols_t& const_symbols_access() const noexcept { return mol_data.atom_symbols; }
        inline main_graph_t& molgraph_access() noexcept { return mol_data.molgr; }
        main_graph_t molgraph_access_py() { return mol_data.molgr; }
        
        void rebuild_active() { mol_data.refiners.rebuild_active(); }
        void refine_dofs() { mol_data.refiners.set_dofs(mol_data); }
        int num_active_dofs() const noexcept { return mol_data.refiners.get_num_active(); }
        inline const auto get_active () const { return mol_data.refiners.get_active(); }
        bool is_xyzless() const noexcept { return mol_data.xyzless_init; }

        inline const std::string get_filename() const noexcept { return this->filename; }

        void require_best_sequence() { dofs.require_best_sequence = true; }

        void customize_sampling_py(py::dict custom_dof_limits) {
            mol_data.refiners.customize_sampling_limits(custom_dof_limits.cast<custom_sampling_limits_t>());
        }

        
        inline void customize_sampling(const custom_limits_data_t& custom_dof_limits) {
            mol_data.refiners.customize_sampling_limits(custom_dof_limits);
        }

        inline custom_limits_data_t get_custom_limits() const { return mol_data.refiners.get_custom_limits(); }

        py::list get_biggest_ringfrag_atoms() const;

        #ifdef KDMOL_LOG
            void assign_unittests(const py::list& input_iterations)
            { mol_data.refiners.assign_unittests(input_iterations); }

            Utils::fulldof_t dof_as_tuple(const int given_free_idx) const {
                for (const auto [ dof_idx, free_idx ] : dofs.free_indexer)
                    if (free_idx == given_free_idx)
                        return {dofs.side_container[dof_idx].first,
                                dofs.bond_container[dof_idx].first,
                                dofs.bond_container[dof_idx].second,
                                dofs.side_container[dof_idx].second};
                throw std::runtime_error(fmt::format("Free dof #{} not found", given_free_idx));
            }

            py::list ddof_as_list(int idx) const {
                auto list_item = py::handle(py::cast(discrete_dofs.owner_container[idx])).cast<py::list>();
                list_item.attr("sort")();
                return list_item;
            }

            py::str get_log() const { return py::cast(Utils::Logger::get_log()); }
        #endif

    private:
        py::module_ nx_o, pyutils_o;
        py::module_ *nx, *pyutils;
        void read_sdf(const std::string& filename);
        void read_graph(py::object input_graph);
        void gen_cyclic_parts();
        void gen_dofs();
        void gen_fragment_graph();
        
        mol_data_t mol_data;
        cp_container_t cyclic_parts;
        dofs_t dofs;
        discrete_dofs_t discrete_dofs;
        dfg_t dfg;
        fragment_container_t frag_container;
        fragref_container_t frag_refs;

        int headnode;
        todo_queue_t todo_queue;

        // Indices of free non-cyclic DOF (own for this class).
        // These will be refined if error happens within Molecule class
        active_idxs_t noncyclic_dof_idxs;

        #ifdef MOL_VALIDATE
            validator_t main_validator;
        #endif
        
        #ifdef MOL_OVERLAP_DETECTION
            overlap_detector_t overlap_detector;
        #endif

        std::string filename;
};

#endif // MOLECULE_H
