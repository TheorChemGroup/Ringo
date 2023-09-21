#ifndef TLCEXTUTILS_H
#define TLCEXTUTILS_H

#include "../utils.h"

namespace TlcUtils {

    class BondItem;
    class AtomContainerDynamic;
    // template<std::size_t N> class AtomContainerStatic;

    class AtomItem {
        public:
            AtomItem() noexcept { }
            inline unsigned int get_num() const noexcept { return num; }
            inline void set_num(unsigned int newnum) noexcept { num = newnum; }

            inline void set_bonds(BondItem *left, BondItem *right) noexcept {
                left_bond = left;
                right_bond = right;
            }
            
            AtomItem operator++(int) noexcept;
            AtomItem operator--(int) noexcept;
            AtomItem operator+(const int& step) const noexcept;
            AtomItem operator-(const int& step) const noexcept;
            AtomItem& operator+=(const int& step) noexcept;
            AtomItem& operator-=(const int& step) noexcept;
            bool operator==(AtomItem const& b) const noexcept { return this->num == b.num; }
            bool operator!=(AtomItem const& b) const noexcept { return this->num != b.num; }

            #ifdef EXT_LOG
                std::string repr();
                py::list pyrepr();
            #endif

            inline void make_bonds_solved() noexcept;
            inline bool can_be_axial() noexcept;
            inline BondItem& get_left_bond() noexcept { return *left_bond; }
            inline BondItem& get_right_bond() noexcept { return *right_bond; }
            inline const BondItem& get_left_bond() const noexcept { return *left_bond; }
            inline const BondItem& get_right_bond() const noexcept { return *right_bond; }
            
        private:
            unsigned int num;
            BondItem *left_bond, *right_bond;
    };


    class BondItem {
        public:
            BondItem() noexcept : free(true), solved(false) { }
            inline unsigned int get_num() const noexcept { return num; }
            inline void set_num(unsigned int newnum) noexcept { num = newnum; }

            inline AtomItem& get_right_atom() noexcept { return *right_atom; }
            inline AtomItem& get_left_atom() noexcept { return *left_atom; }

            BondItem operator++(int) noexcept;
            BondItem operator--(int) noexcept;
            bool operator==(BondItem const& b) const noexcept { return this->num == b.num; }
            bool operator!=(BondItem const& b) const noexcept { return this->num != b.num; }

            inline void set_atoms(AtomItem *left, AtomItem *right) noexcept {
                left_atom = left;
                right_atom = right;
            }
            bool free, solved;

            #ifdef EXT_LOG
                std::string repr();
                py::list pyrepr();
            #endif

            // operator++ , operator-- , 
            // resursive : operator+ , operator-

            inline std::pair<unsigned int, unsigned int> atom_idxs() noexcept
                { return std::pair<unsigned int, unsigned int>(left_atom->get_num(), right_atom->get_num()); }

            inline std::pair<unsigned int, unsigned int> atom_idxs_reversed() noexcept 
                { return std::pair<unsigned int, unsigned int>(right_atom->get_num(), left_atom->get_num()); }

        private:
            unsigned int num;
            AtomItem *left_atom, *right_atom;
    };


    template<typename T>
    class CyclicGroup {
        public:
            CyclicGroup (const int& size) noexcept : 
                ord(size),
                items(std::vector<T>(size))
            { 
                unsigned int i = 0;
                for(auto item : items) {
                    items[i].set_num(i);
                    i += 1;
                }
            }

            // Impose cyclicity
            inline T& operator [](int idx) noexcept {
                while (idx < 0)
                    idx += ord;
                return items[idx % ord];
            } 
            
            inline const T& operator[](int idx) const noexcept {
                while (idx < 0)
                    idx += ord;
                return items[idx % ord];
            }

            inline typename std::vector<T>::iterator begin() noexcept { return items.begin(); }
            inline typename std::vector<T>::iterator end() noexcept { return items.end(); }
            inline typename std::vector<T>::const_iterator cbegin() noexcept { return items.cbegin(); }
            inline typename std::vector<T>::const_iterator cend() noexcept { return items.cend(); }

            inline unsigned int size() const noexcept { return ord; }

            unsigned int get_distance(const unsigned int& a, const unsigned int& b) const noexcept {
                unsigned int dist = abs(static_cast<int>(b - a)) - 1;
                if (b < a)
                    dist = ord - dist - 2;
                return dist;
            }
        
        protected:
            unsigned int ord;
            std::vector<T> items;
    };


    class AtomGroupSimple : public CyclicGroup<AtomItem> {
        public:
            AtomGroupSimple(const int& n_atoms) noexcept
                : CyclicGroup<AtomItem>(n_atoms)
            { }

            AtomContainerDynamic collect_proper_atoms() noexcept;
    };


    class BondGroupSimple : public CyclicGroup<BondItem> {
        public:
            BondGroupSimple(const int& n_atoms) noexcept
                : CyclicGroup<BondItem>(n_atoms)
            { }

            #ifdef EXT_LOG
                std::string repr_free_list() const {
                    py::list res;
                    for(const auto& item : items)
                        if(item.free)
                            res.append(1);
                        else
                            res.append(0);
                    return py::repr(res).cast<std::string>();
                }

                std::string repr_solved_list() const {
                    py::list res;
                    for(const auto& item : items)
                        if(item.solved)
                            res.append(1);
                        else
                            res.append(0);
                    return py::repr(res).cast<std::string>();
                }
            #endif
    };


    class AtomGroup : public AtomGroupSimple {
        public:
            AtomGroup() noexcept : AtomGroupSimple(0) { }

            AtomGroup(const int& n_atoms) noexcept
                : AtomGroupSimple(n_atoms)
            { }

            void set_vangle_data(const std::vector<double>& new_vangles) {
                if (this->ord != new_vangles.size())
                    throw std::runtime_error(fmt::format("There must be {} valence angles. ({} were provided)",
                                                         this->ord, new_vangles.size()));
                vangles = new_vangles;
            }

            inline void init_tlc_map(std::array<unsigned int, 3> axidx) noexcept {
                tlc_map[0] = axidx[0];
                tlc_map[3] = axidx[1];
                tlc_map[6] = axidx[2];

                tlc_map[1] = this->operator[](axidx[0] + 1).get_num();
                tlc_map[4] = this->operator[](axidx[1] + 1).get_num();
                tlc_map[7] = this->operator[](axidx[2] + 1).get_num();
                
                tlc_map[8] = this->operator[](axidx[0] - 1).get_num();
                tlc_map[2] = this->operator[](axidx[1] - 1).get_num();
                tlc_map[5] = this->operator[](axidx[2] - 1).get_num();
            }

            // #ifdef EXT_LOG
                py::list list_vangles() {
                    py::list l;
                    for (int i = 0; i < this->items.size(); ++i) {
                        py::list newitem;
                        newitem.append(this->items[i].get_num());
                        newitem.append(vangles[i]);
                        l.append(newitem);
                    }
                    return l;
                }

                std::string repr_vangles() {
                    return py::repr(list_vangles()).cast<std::string>();
                }
                
                std::string repr_tlc_vangles() {
                    py::list l;
                    for(const auto& item : tlc_vangles)
                        l.append(item);
                    return py::repr(l).cast<std::string>();
                }
            // #endif

            // Containers for valence angles of the original cycle
            std::vector<double> vangles;

            // Container for valence angles of the effective cycle which will be passed to TLCSolver
            std::array<double, 9> tlc_vangles;

            // Map of indices from TLC to original cycle
            std::array<unsigned int, 9> tlc_map;

        private:

            // TODO Container for cartesian coords
    };
    

    class BondGroup : public BondGroupSimple {
        public:
            BondGroup() noexcept : BondGroupSimple(0) { }

            BondGroup(const int& n_atoms) noexcept
                : BondGroupSimple(n_atoms),
                nfixed(-1), nvariable(-1)
            { }

            void markdown_dihedrals(std::unordered_map<std::string, std::vector<std::pair<unsigned int, unsigned int>>> dof_data, AtomGroup& catoms) {
                nfixed = 0;
                nvariable = 0;
                for(const auto& bond : *this) {
                    if(!bond.free)
                        nfixed++;
                    else if(bond.free && !bond.solved)
                        nvariable++;
                }

                if(nfixed != dof_data["fixed"].size())
                    throw std::runtime_error(fmt::format("Error while parsing fixed dihedrals. {} were provided but {} were parsed",
                                                         dof_data["fixed"].size(), nfixed));
                
                if(nvariable != dof_data["variable"].size())
                    throw std::runtime_error(fmt::format("Error while parsing variable dihedrals. {} were provided but {} were parsed",
                                                         dof_data["variable"].size(), nvariable));
                
                dihedral_permutation.reserve(this->ord);
                #ifdef EXT_LOG
                    if (dihedral_permutation.capacity() != this->ord)
                        throw std::runtime_error(fmt::format("Capacity = {} but requested {} elements", dihedral_permutation.capacity(), this->ord));
                #endif

                std::vector<unsigned int> dih_prototype;
                dih_prototype.reserve(nvariable + nfixed);
                for(const auto& pair : dof_data["variable"])
                    dih_prototype.push_back(pair.first);
                for(const auto& pair : dof_data["fixed"])
                    dih_prototype.push_back(pair.first);

                for(const auto& bond : *this) {
                    if(bond.solved)
                        dihedral_permutation.push_back(-1);
                    else {
                        auto find_iter = std::find(dih_prototype.cbegin(), dih_prototype.cend(), bond.get_num());
                        #ifdef EXT_LOG
                            if (find_iter == dih_prototype.cend())
                                throw std::runtime_error(fmt::format("Index list {} does not contain index {}", 
                                                         py::repr(py::cast(dih_prototype)).cast<std::string>(), bond.get_num()));
                        #endif
                        dihedral_permutation.push_back(static_cast<int>(find_iter - dih_prototype.cbegin()));
                    }
                }
            }

            void set_length_data(const std::vector<double>& new_lengths) {
                if (this->ord != new_lengths.size())
                    throw std::runtime_error(fmt::format("Bond length array must have {} elements. ({} were provided)",
                                                         this->ord, new_lengths.size()));
                this->bonds = new_lengths;
            }
            
            void set_fixed_dihedral_data(const std::vector<double>& new_dihedrals) {
                if (nfixed != new_dihedrals.size())
                    throw std::runtime_error(fmt::format("There must be {} fixed dihedrals. ({} were provided)",
                                                         nfixed, new_dihedrals.size()));
                this->dihedrals = std::vector<double>(nvariable);
                std::copy(new_dihedrals.begin(), new_dihedrals.end(), std::back_inserter(this->dihedrals));
            }
            
            void set_variable_dihedral_data(const std::vector<double>& new_dihedrals) { // OPTIMIZE THIS ONE
                if (nvariable != new_dihedrals.size())
                    throw std::runtime_error(fmt::format("There must be {} variable dihedrals. ({} were provided)",
                                                         nvariable, new_dihedrals.size()));
                std::copy(new_dihedrals.begin(), new_dihedrals.end(), this->dihedrals.begin());
            }
            
            inline double get_dihedral(const int& num) const noexcept {
                return dihedrals[dihedral_permutation[num]];
            }

            // #ifdef EXT_LOG
                py::list list_lengths() {
                    py::list l;
                    for(const auto& item : bonds)
                        l.append(item);
                    return l;
                }

                std::string repr_lengths() {
                    return py::repr(list_lengths()).cast<std::string>();
                }
                
                std::string repr_tlc_lengths() {
                    py::list l;
                    for(const auto& item : tlc_lengths)
                        l.append(item);
                    return py::repr(l).cast<std::string>();
                }

                std::string repr_dihedrals_init() {
                    py::list l;
                    for(int i = 0; i < nvariable; ++i)
                        dihedrals[i] = 0.0;
                    for(const auto& item : dihedrals)
                        l.append(item);
                    return py::repr(l).cast<std::string>();
                }

                py::list list_dihedrals_full() {
                    py::list l;
                    for(const auto& item : dihedrals)
                        l.append(item);
                    return l;
                }
                
                std::string repr_dihedrals_full() {
                    return py::repr(list_dihedrals_full()).cast<std::string>();
                }
                
                std::string repr_dihedral_perm() {
                    py::list l;
                    for(const auto& item : dihedral_permutation)
                        l.append(item);
                    return py::repr(l).cast<std::string>();
                }

                std::string repr_tlc_dihedrals() {
                    py::list l;
                    for(const auto& item : tlc_dihedrals)
                        l.append(item);
                    return py::repr(l).cast<std::string>();
                }
            // #endif

            // Containers for parameters of the original cycle
            // The 'dihedrals' vector is constructed as follows: ---VARIABLE---|---FIXED---
            std::vector<double> bonds;
            std::vector<double> dihedrals;

            // Containers for parameters of the effective cycle which will be passed to TLCSolver
            std::array<double, 9> tlc_lengths;
            std::array<double, 3> tlc_dihedrals;

        private:
            int nfixed, nvariable;

            // Contains indices of bonds in 'dihedrals' vector and -1 for solved bonds:
            std::vector<int> dihedral_permutation;
    };


    template <class Container_type>
    class AtomContainer {
        public:
            AtomContainer() noexcept { }

            inline unsigned int size() const noexcept { return items.size(); }
            inline void set(AtomItem* x, const int& i) noexcept { items[i] = x; }
            
            bool contains(const AtomItem& x) {
                const auto x_idx = x.get_num();
                for(const auto& item : *this)
                    if (item->get_num() == x_idx)
                        return true;
                return false;
            }

            bool contains(const AtomItem* x) {
                const auto x_idx = x->get_num();
                for(const auto& item : *this)
                    if (item->get_num() == x_idx)
                        return true;
                return false;
            }

            #ifdef EXT_LOG
                std::string repr() {
                    std::string res = fmt::format("Atom container with {} items: ", items.size());
                    for(const auto& item : *this)
                        res += fmt::format("{} ", item->get_num());
                    std::cout << res << std::endl;
                    return res;
                }

                std::string repr_list() {
                    py::list l;
                    for(const auto& item : *this)
                        l.append(item->get_num());
                    return py::repr(l).cast<std::string>();
                }
                
                std::string repr_list(const int& n_left) {
                    py::list l;
                    int i = 0;
                    for(const auto& item : *this) {
                        if(i++ >= n_left)
                            break;
                        l.append(item->get_num());
                    }
                    return py::repr(l).cast<std::string>();
                }
            #endif

            inline AtomItem& operator [](int idx) noexcept { return *items[idx]; } // DO NOT USE IT FOR ASSIGNMENT!!!
            inline const AtomItem& operator[](int idx) const noexcept { return *items[idx]; }
            inline typename Container_type::iterator begin() noexcept { return items.begin(); }
            inline typename Container_type::iterator end() noexcept { return items.end(); }
            inline typename Container_type::const_iterator cbegin() const noexcept { return items.cbegin(); }
            inline typename Container_type::const_iterator cend() const noexcept { return items.cend(); }
            
        protected:
            Container_type items;
    };


    class AtomContainerDynamic : public AtomContainer<std::vector<AtomItem*>> {
        public:
            AtomContainerDynamic() noexcept {}

            AtomContainerDynamic(int prealloc_size) noexcept {
                if (prealloc_size > 0) {
                    this->items.reserve(prealloc_size);
                }
            }

            AtomContainerDynamic(const std::vector<unsigned int>& atom_data, AtomGroupSimple& catoms) noexcept {
                for(const auto& atom_idx : atom_data)
                    this->items.push_back(&catoms[atom_idx]);
            }

            AtomContainerDynamic& operator+=(AtomItem& rhs) noexcept {
                this->items.push_back(&rhs);
                return *this;
            }
    };


    template<std::size_t N>
    class AtomContainerStatic : public AtomContainer<std::array<AtomItem*, N>> {
        public:
            AtomContainerStatic() noexcept {}

            template<class AtomGrT>
            AtomContainerStatic(const std::vector<unsigned int>& atom_data, AtomGrT& catoms) {
                int i = 0;
                for(const auto& atom_idx : atom_data) {
                    if(i >= N)
                        throw std::runtime_error("Overflow of static-sized array");
                    this->items[i++] = &catoms[atom_idx];
                }
            }

            template<class AtomGrT>
            AtomContainerStatic(const std::array<unsigned int, N>& atom_data, AtomGrT& catoms) noexcept {
                int i = 0;
                for(const auto& atom_idx : atom_data)
                    this->items[i++] = &catoms[atom_idx];
            }
    };


    template <class Container_type>
    class BondContainer {
        public:
            BondContainer() noexcept { }

            inline unsigned int size() const noexcept { return items.size(); }
            inline void set(BondItem* x, const int& i) noexcept { items[i] = x; }

            bool contains(const BondItem& x) {
                const auto x_idx = x.get_num();
                for(const auto& item : *this)
                    if (item->get_num() == x_idx)
                        return true;
                return false;
            }

            #ifdef EXT_LOG
                std::string repr() {
                    std::string res = fmt::format("Bond container with {} items: ", items.size());
                    for(const auto& item : *this)
                        res += fmt::format("{} ", item->get_num());
                    std::cout << res << std::endl;
                    return res;
                }

                std::string repr_list() {
                    py::list l;
                    for(const auto& item : *this)
                        l.append(item->get_num());
                    return py::repr(l).cast<std::string>();
                }

                std::string repr_list(const int& n_left) {
                    py::list l;
                    int i = 0;
                    for(const auto& item : *this) {
                        if(i++ >= n_left)
                            break;
                        l.append(item->get_num());
                    }
                    return py::repr(l).cast<std::string>();
                }
            #endif

            inline BondItem& operator [](int idx) noexcept { return *items[idx]; } // DO NOT USE IT FOR ASSIGNMENT!!!
            inline const BondItem& operator[](int idx) const noexcept { return *items[idx]; }
            inline typename Container_type::iterator begin() noexcept { return items.begin(); }
            inline typename Container_type::iterator end() noexcept { return items.end(); }
            inline typename Container_type::const_iterator cbegin() const noexcept { return items.cbegin(); }
            inline typename Container_type::const_iterator cend() const noexcept { return items.cend(); }
            
        protected:
            Container_type items;
    };


    class BondContainerDynamic : public BondContainer<std::vector<BondItem*>> {
        public:
            BondContainerDynamic() noexcept {}

            BondContainerDynamic(int prealloc_size) noexcept {
                if (prealloc_size > 0) {
                    this->items.reserve(prealloc_size);
                }
            }

            BondContainerDynamic(const std::vector<unsigned int>& bond_data, BondGroupSimple& cbonds) noexcept {
                for(const auto& bond_idx : bond_data)
                    this->items.push_back(&cbonds[bond_idx]);
            }

            BondContainerDynamic(const std::vector<std::pair<unsigned int, unsigned int>> bond_data, BondGroupSimple& cbonds) {
                for(auto& bond : cbonds) {
                    std::pair<unsigned int, unsigned int> cur_atoms = bond.atom_idxs();
                    std::pair<unsigned int, unsigned int> cur_atoms_rev = bond.atom_idxs_reversed();
                    if   (std::find(bond_data.cbegin(), bond_data.cend(), cur_atoms) != bond_data.cend() || 
                          std::find(bond_data.cbegin(), bond_data.cend(), cur_atoms_rev) != bond_data.cend())
                        this->items.push_back(&bond);
                }
                if (this->items.size() != bond_data.size())
                    throw std::runtime_error(fmt::format("Some bonds weren't found. Required = {}.", py::repr(py::cast(bond_data)).cast<std::string>()));
            }

            BondContainerDynamic& operator+=(BondItem& rhs) noexcept {
                this->items.push_back(&rhs);
                return *this;
            }
    };


    template<std::size_t N>
    class BondContainerStatic : public BondContainer<std::array<BondItem*, N>> {
        public:
            BondContainerStatic() noexcept {}

            BondContainerStatic(const std::vector<unsigned int>& bond_data, BondGroupSimple& cbonds) {
                int i = 0;
                for(const auto& bond_idx : bond_data) {
                    if(i >= N)
                        throw std::runtime_error("Overflow of static-sized array");
                    this->items[i++] = &cbonds[bond_idx];
                }
            }
    };


    class AxialAtomContainer : public AtomContainerStatic<3> {
        public:
            AxialAtomContainer() noexcept : n_set(0) {}
            inline void add_atom(AtomItem& a) noexcept { this->items[n_set++] = &a; }
            inline bool is_full() const noexcept { return n_set == 3; }
        
            unsigned int get_min_metric(const AtomGroupSimple& catoms) const {
                int a = this->items[0]->get_num(),
                    b = this->items[1]->get_num(),
                    c = this->items[2]->get_num();
                #ifdef EXT_LOG
                    if (!(catoms.size() == catoms.get_distance(a, b) + catoms.get_distance(b, c) + catoms.get_distance(c, a) + 3))
                        throw std::runtime_error("Distance calc method is bugged");
                #endif
                auto res = std::min({catoms.get_distance(a, b),
                                 catoms.get_distance(b, c),
                                 catoms.get_distance(c, a)});
                // fmt::print(" Minmetric = {} ", res);
                return res;
            }

            bool is_sensible(const AtomGroupSimple& catoms) const noexcept { return get_min_metric(catoms) > 0; }

            double get_metric(const AtomGroupSimple& catoms) const {
                int a = this->items[0]->get_num(),
                    b = this->items[1]->get_num(),
                    c = this->items[2]->get_num();
                #ifdef EXT_LOG
                    if (!(catoms.size() == catoms.get_distance(a, b) + catoms.get_distance(b, c) + catoms.get_distance(c, a) + 3))
                        throw std::runtime_error("Distance calc method is bugged");
                #endif
                auto res = sqrt(pow(catoms.get_distance(a, b), 2) +
                            pow(catoms.get_distance(b, c), 2) +
                            pow(catoms.get_distance(c, a), 2));
                // fmt::print(" Metric = {} ", res);
                return res;
            }

        private:
            int n_set;
    };


    inline void AtomItem::make_bonds_solved() noexcept {
        left_bond->solved = true;
        right_bond->solved = true;
    }


    inline bool AtomItem::can_be_axial() noexcept {
        return left_bond->free && right_bond->free;
    }


    void connect(AtomGroupSimple&, BondGroupSimple&) noexcept;
    #ifdef EXT_LOG
        std::string tuplelist_to_str(py::list x);
        void check_topology(AtomGroupSimple&, BondGroupSimple&);
        std::string repr_matrix_buffer(const double* buf, const int& size1, const int& size2);

        using Matrix = boost::numeric::ublas::matrix<double>;
        std::string repr_matrix(const Matrix& m);
        std::string repr_row(const Matrix& m, const int& row_idx);
        std::string repr_matrix_buffer_addr(const double* buf, const std::vector<std::size_t>& addrs, const int& size2);
    #endif


    // CAUTION!!! The following implementations are tricky
    inline double vangle_from_dirsA(const double* r1, const double* r2) noexcept {
        double res = cblas_ddot(3, r1, 1, r2, 1) / cblas_dnrm2(3, r2, 1);
        // res = std::copysign(std::min(std::abs(res), 1.0), res);
        res = acos(res);
        return res;
    }

    inline double vangle_from_dirsB(const double* r1, const double* r2) noexcept {
        double res = -cblas_ddot(3, r1, 4, r2, 1) / cblas_dnrm2(3, r2, 1);
        // res = std::copysign(std::min(std::abs(res), 1.0), res);
        res = acos(res);
        return res;
    }

    inline void calc_cross(const double *v_A, const double *v_B, double *c_P) noexcept {
        c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
        c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
        c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
    }
    
    inline void calc_crossL(const double *v_A, const double *v_B, double *c_P) noexcept {
        c_P[0] = v_A[4] * v_B[2] - v_A[8] * v_B[1];
        c_P[1] = -(v_A[0] * v_B[2] - v_A[8] * v_B[0]);
        c_P[2] = v_A[0] * v_B[1] - v_A[4] * v_B[0];
    }

    inline void calc_crossR(const double *v_A, const double *v_B, double *c_P) noexcept {
        c_P[0] = v_A[1] * v_B[8] - v_A[2] * v_B[4];
        c_P[1] = -(v_A[0] * v_B[8] - v_A[2] * v_B[0]);
        c_P[2] = v_A[0] * v_B[4] - v_A[1] * v_B[0];
    }

    inline double dihedral_from_dirs(const double *r1, const double *r2, const double *r3) noexcept {
        double p[3], q[3], s[3];
        TlcUtils::calc_cross(r1, r2, p);
        TlcUtils::calc_crossR(r2, r3, q);
        TlcUtils::calc_crossL(r3, r1, s);
        
        double res = -cblas_ddot(3, &p[0], 1, &q[0], 1);
        res /= cblas_dnrm2(3, &q[0], 1);
        res /= cblas_dnrm2(3, &p[0], 1);
        
        // res = std::copysign(std::min(std::abs(res), 1.0), res);
        res = std::copysign(acos(res), cblas_ddot(3, s, 1, r2, 1));
        return res;
    }

}

#endif // TLCEXTUTILS_H
