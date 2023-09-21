#ifndef UTILS_H
#define UTILS_H

#include "sideheaders.h"

constexpr double H2KC = 627.509474063;
constexpr double DEG2RAD = 0.0174532925199432957692;
constexpr double RAD2DEG = 1 / DEG2RAD;

#define assertm(expr, message) \
    do { \
        if (!(expr)) { \
            throw std::runtime_error(message); \
        } \
    } while (0)
// #define assertm(expr, message) \
//     do { \
//         if (!(expr)) { \
//             throw std::runtime_error(fmt::format("Error: \"{}\" Assertion failed at {}:{}", message, __FILE__, __LINE__)); \
//         } \
//     } while (0)

#ifdef USE_CHECKPOINTS
    #define checkpoint(message) \
        do { \
            fmt::print("Reached checkpoint: \"{}\". File {}:{}", message, __FILE__, __LINE__); std::cout << std::endl; \
        } while (0)
#else
    // Cut out in release version
    #define checkpoint(message)
#endif

#if defined(OVERLAP_DETECTION) || defined(OVERLAP_DETECTION_FINAL)
    #define BUILD_OVERLAP_DETECTION
#endif

namespace py = pybind11;
using boost::typeindex::type_id_with_cvr;

namespace Utils {
    constexpr double STRONG_THRESHOLD = 0.01;
    constexpr double WEAK_THRESHOLD = 0.05;
    #if defined(VALIDATION)
        constexpr double THRESHOLD = STRONG_THRESHOLD;
    #elif defined(WEAK_VALIDATION)
        constexpr double THRESHOLD = WEAK_THRESHOLD;
    #endif

    using umatrix_t = boost::numeric::ublas::matrix<double>;
    using umatrix4d_t = boost::numeric::ublas::fixed_matrix<double, 4, 4>;
    using umatrix3d_t = boost::numeric::ublas::fixed_matrix<double, 3, 3>;

    using uvector_t = boost::numeric::ublas::vector<double>;
    using uvector3d_t = boost::numeric::ublas::fixed_vector<double, 3>;
    using uvector4d_t = boost::numeric::ublas::fixed_vector<double, 4>;

    using dof_t = std::pair<unsigned int, unsigned int>;
    using frame_t = std::tuple<unsigned int, unsigned int, unsigned int>;
    using fulldof_t = std::tuple<unsigned int, unsigned int, unsigned int, unsigned int>;
    
    class PyModules {
        public:
            static inline py::module_* get_nx()
            { return nx; }

            static inline void set_nx(py::module_* new_nx)
            { nx = new_nx; }

            static inline py::module_* get_pyutils()
            { return pyutils; }

            static inline void set_pyutils(py::module_* new_pyutils)
            { pyutils = new_pyutils; }

        private:
            PyModules() {}
            PyModules(PyModules const&) = delete;
            void operator=(PyModules const&) = delete;
            static py::module_* nx;
            static py::module_* pyutils;
    };

    template <class edgeidx_t, class fragidx_t>
    struct EdgeSpec {
        edgeidx_t index;
        fragidx_t parent_frag;
    };

    struct EdgeIndexer { // Implements a "good" way to lookup indices of edges
        using atom_t = int;
        using edgeidx_t = int;
        using fragidx_t = int;
        using spec_t = EdgeSpec<edgeidx_t, fragidx_t>;

        using map_t = std::map<atom_t, spec_t>;
        using container_t = std::vector<map_t>;

        EdgeIndexer() { }

        template <class initial_container_t> // initial_container_t = std::vector<std::pair<atom_t, atom_t>>
        EdgeIndexer(const int& natoms, const initial_container_t& raw_data) {
            data = container_t(natoms);
            
            for(int i = 0; i < raw_data.size(); ++i) {
                auto pair = raw_data[i];
                data[pair.first][pair.second] = {i, -1};
                data[pair.second][pair.first] = {i, -1};
            }
        }

        inline const map_t& atom_data(const atom_t& index) const noexcept {
            return data[index];
        }
        
        inline map_t& set_atom_data(const atom_t& index) noexcept {
            return data[index];
        }

        void validate_frag_indices() const {
            bool okay = true;
            for (int i = 0; i < data.size(); ++i)
                for (const auto& [nb_idx, spec] : data[i])
                    if (spec.parent_frag == -1) {
                        okay = false;
                        fmt::print("Bond {}-{} was not assigned to any fragment\n", i, nb_idx);
                    }
            if (!okay)
                throw std::runtime_error("Validation failed");
        }

        private:
            container_t data;
    };

    template <typename C, typename P> 
    inline void erase_remove_if(C& c, P predicate) noexcept {
        c.erase(std::remove_if(c.begin(), c.end(), predicate), c.end());
    }

    struct InvalidChar {
        bool operator()(char c) const {
            return !isprint(static_cast<unsigned char>(c));
        }
    };

    inline double get_distance_raw(const double* a_xyz, const double* b_xyz) noexcept {
        double res[3];
        cblas_dcopy(3, a_xyz, 1, &res[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &res[0], 1);
        return cblas_dnrm2(3, &res[0], 1);
    }

    inline double get_vangle_raw(const double* a_xyz, const double* b_xyz, const double* c_xyz) noexcept {
        double dirA[3];
        double dirC[3];
        cblas_dcopy(3, a_xyz, 1, &dirA[0], 1);
        cblas_dcopy(3, c_xyz, 1, &dirC[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &dirA[0], 1);
        cblas_daxpy(3, -1.0, b_xyz, 1, &dirC[0], 1);
        
        double normA = cblas_dnrm2(3, &dirA[0], 1);
        double normC = cblas_dnrm2(3, &dirC[0], 1);
        cblas_dscal(3, 1 / normA, &dirA[0], 1);
        cblas_dscal(3, 1 / normC, &dirC[0], 1);
        return acos(cblas_ddot(3, &dirA[0], 1, &dirC[0], 1)) * RAD2DEG;
    }

    inline double get_dihedral_raw(const double* a_xyz, const double* b_xyz, const double* c_xyz, const double* d_xyz) noexcept {
        using namespace boost::numeric::ublas;

        uvector3d_t r_a, r_b, r_c, r_d;
        for (int i = 0; i < 3; ++i) {
            r_a(i) = a_xyz[i];
            r_b(i) = b_xyz[i];
            r_c(i) = c_xyz[i];
            r_d(i) = d_xyz[i];
        }

        uvector3d_t fr1_side = r_a - r_b;
        uvector3d_t fr1_mid = r_c - r_b;
        uvector3d_t fr2_mid = -fr1_mid;
        uvector3d_t fr2_side = r_d - r_c;
        fr1_side -= (inner_prod(fr1_side, fr1_mid) / inner_prod(fr1_mid, fr1_mid)) * fr1_mid;
        fr2_side -= (inner_prod(fr2_side, fr2_mid) / inner_prod(fr2_mid, fr2_mid)) * fr2_mid;
        fr1_side /= norm_2(fr1_side);
        fr2_side /= norm_2(fr2_side);

        auto dotprod = inner_prod(fr1_side, fr2_side);
        if (dotprod >= 1.0)
            dotprod = 0.99999999999;
        else if (dotprod <= -1.0)
            dotprod = -0.99999999999;
        auto ang = acos(dotprod);
        
        umatrix3d_t orient;
        for (int i = 0; i < 3; ++i) {
            orient(i, 0) = fr1_side(i);
            orient(i, 1) = fr1_mid(i);
            orient(i, 2) = fr2_side(i);
        }
        
        auto det =  orient(0, 0) * (orient(1, 1) * orient(2, 2) - orient(2, 1) * orient(1, 2))
                  - orient(0, 1) * (orient(1, 0) * orient(2, 2) - orient(2, 0) * orient(1, 2))
                  + orient(0, 2) * (orient(1, 0) * orient(2, 1) - orient(2, 0) * orient(1, 1));
        if (det < 0)
            ang = -ang;
        return ang;
    }

    inline double standardize_dihedral(double dih) {
        while (dih < -M_PI)
            dih += 2*M_PI;
        while (dih > M_PI)
            dih -= 2*M_PI;
        return dih;
    }

    template <class VType>
    py::list vector_to_list(const VType& v) {
        py::list res;
        for (unsigned int i = 0; i < v.size(); ++ i)
            res.append(v[i]);
        return res;
    }

    
    template <class MType>
    py::list matrix_to_list(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.size1 (); ++ i) {
            py::list temp;
            for (unsigned int j = 0; j < m.size2 (); ++ j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return res;
    }

    template <class VType>
    std::string repr_vector(const VType& v) {
        py::list res;
        for (unsigned int i = 0; i < v.size(); ++ i)
            res.append(v[i]);
        return py::repr(res).cast<std::string>();
    }
    
    template <class MType>
    std::string repr_matrix(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.size1 (); ++ i) {
            py::list temp;
            for (unsigned int j = 0; j < m.size2 (); ++ j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }

    class BasicCoordinates {
        public:
            BasicCoordinates() { }
            
            virtual const std::string repr() const noexcept = 0;
            virtual inline const double* get_atom_raw(const int index) const noexcept = 0;

            inline double get_distance(const int a_idx, const int b_idx) const noexcept
            { return Utils::get_distance_raw(this->get_atom_raw(a_idx),
                                             this->get_atom_raw(b_idx)); }

            inline double get_distance(const std::pair<int, int> idxs) const noexcept
            { return Utils::get_distance_raw(this->get_atom_raw(idxs.first),
                                             this->get_atom_raw(idxs.second)); }
            
            inline double get_vangle(const int a_idx, const int b_idx, const int c_idx) const noexcept
            { return Utils::get_vangle_raw(this->get_atom_raw(a_idx),
                                           this->get_atom_raw(b_idx),
                                           this->get_atom_raw(c_idx)); }
                        
            inline double get_dihedral(const int a_idx, const int b_idx, const int c_idx, const int d_idx) const noexcept
            { return Utils::get_dihedral_raw(this->get_atom_raw(a_idx),
                                             this->get_atom_raw(b_idx),
                                             this->get_atom_raw(c_idx),
                                             this->get_atom_raw(d_idx)); }
    };

    template <class M>
    class MapXyzContainer : public BasicCoordinates {
        protected:
            using parent_t = MapXyzContainer<M>;
            using core_container_t = parent_t;
            using base_t = M;
            base_t xyz_;

        public:
            using xyz_t = typename M::mapped_type;
            using iterator = typename base_t::iterator;
            using const_iterator = typename base_t::const_iterator;

            MapXyzContainer() { }

            template <class C> MapXyzContainer(const C& idx_container) {
                for (const auto& atom : idx_container)
                    xyz_[atom] = xyz_t();
            }

            const std::string repr() const noexcept override {
                py::list atoms, res;
                for (const auto& [idx, _] : xyz_)
                    atoms.append(idx);
                atoms.attr("sort")();

                for (const auto atom : atoms)
                    res.append(vector_to_list(xyz_.at(py::handle(atom).cast<int>())));
                return py::repr(res).cast<std::string>();
            }

            inline xyz_t& operator[](const int& index) noexcept
            { return xyz_[index]; }
            inline xyz_t& at(const int& index) noexcept
            { return xyz_.at(index); }
            inline int count(const int& index) noexcept
            { return xyz_.count(index); }

            inline const double* get_atom_raw(const int index) const noexcept override
            { return &xyz_.at(index)(0); }

            void set_atom(const int& index, std::initializer_list<double> new_coords)
            { std::copy(new_coords.begin(), new_coords.end(), &xyz_[index](0)); }

            inline int size() const noexcept { return xyz_.size(); }

            inline const std::vector<int> get_atoms() const noexcept {
                std::vector<int> res;
                res.reserve(xyz_.size());
                for (const auto& [idx, _] : xyz_)
                    res.push_back(idx);
                return res;
            }

            inline typename base_t::iterator begin() noexcept { return xyz_.begin(); }
            inline typename base_t::iterator end() noexcept { return xyz_.end(); }
            inline typename base_t::const_iterator cbegin() const noexcept { return xyz_.cbegin(); }
            inline typename base_t::const_iterator cend() const noexcept { return xyz_.cend(); }

            inline typename base_t::const_iterator find(const int& index) const noexcept { return xyz_.find(index); }
    };

    template <unsigned int DIM>
    class BoostXyzContainer : public BasicCoordinates {
        protected:
            using parent_t = BoostXyzContainer<DIM>;
            using core_container_t = parent_t;
            using base_t = umatrix_t;
            base_t xyz_;
        
        public:
            using xyz_t = boost::numeric::ublas::matrix_row<base_t>;
            
            BoostXyzContainer() {}
            BoostXyzContainer(unsigned int natoms) : xyz_(base_t(natoms, DIM)) { }
            
            const std::string repr() const noexcept override {
                return repr_matrix(xyz_);
            }

            inline xyz_t operator[](const int& index) noexcept
            { return xyz_t(xyz_, index); }
            inline xyz_t at(const int& index) noexcept
            { return xyz_t(xyz_, index); }

            inline const double* get_atom_raw(const int index) const noexcept override
            { return &xyz_(index, 0); }

            inline base_t& to_boost_format() noexcept 
            { return xyz_; }
            
            inline const base_t& to_boost_format() const noexcept 
            { return xyz_; }

            void set_atom(const int& index, std::initializer_list<double> new_coords)
            { std::copy(new_coords.begin(), new_coords.end(), &xyz_(index, 0)); }

            inline int size() const noexcept { return xyz_.size1(); }

            inline const std::vector<int> get_atoms() const noexcept {
                std::vector<int> res;
                res.reserve(xyz_.size1());
                for (int i = 0; i < xyz_.size1(); ++i)
                    res.push_back(i);
                return res;
            }
    };

    template <class C>
    class Boost4DAdapter : public C {
        protected:
            using parent_t = C;
            using core_container_t = typename C::core_container_t;

        public:
            using xyz_t = boost::numeric::ublas::matrix_vector_slice<typename C::core_container_t::base_t>;
            using fourD_t = typename C::xyz_t;
            
            template <class ...Args>
            Boost4DAdapter(Args &&... args) : parent_t(args...) { }

            // For treating it like a 3D vector
            inline xyz_t operator[](const int& index) noexcept
            { return xyz_t(C::xyz_, boost::numeric::ublas::slice(index, 0, 3), boost::numeric::ublas::slice(0, 1, 3)); }
            inline xyz_t at(const int& index) noexcept
            { return xyz_t(C::xyz_, boost::numeric::ublas::slice(index, 0, 3), boost::numeric::ublas::slice(0, 1, 3)); }
            
            // For access as a 4D vector
            inline fourD_t get_4d(const int& index) noexcept
            { return C::at(index); }
    };

    template<class T>
    inline auto convert_to_py(const T& obj) {
        if constexpr(std::is_base_of_v<py::object, T>)
            return obj;
        else
            return py::cast(obj);
    }

    template <typename T>
    py::list as_list(const T& x) {
        py::list res;
        res.append(convert_to_py(x));
        return res;
    }

    template<typename T, typename = void>
    struct has_value_type : std::false_type {};
    template<typename T>
    struct has_value_type<T, std::void_t<typename T::value_type>> : std::true_type {};
    template<typename T, typename = void>
    struct has_array_type : std::false_type {};
    template<typename T>
    struct has_array_type<T, std::void_t<typename T::array_type>> : std::true_type {};

    template <typename T>
    const std::string repr(T obj) {
        if constexpr(has_value_type<T>::value && has_array_type<T>::value)
        {
            // Matrices
            if constexpr(std::is_base_of_v<umatrix4d_t, T>) {
                return repr_matrix(obj);
            } else if constexpr(std::is_base_of_v<umatrix3d_t, T>) {
                return repr_matrix(obj);
            } else if constexpr(std::is_base_of_v<umatrix_t, T>) {
                return repr_matrix(obj);
            }
            // Vectors
            else if constexpr(std::is_base_of_v<uvector3d_t, T>) {
                return repr_vector(obj);
            } else if constexpr(std::is_base_of_v<uvector4d_t, T>) {
                return repr_vector(obj);
            } else if constexpr(std::is_base_of_v<uvector_t, T>) {
                return repr_vector(obj);
            }
        }
        // Coord containers
        else if constexpr(std::is_base_of_v<BasicCoordinates, T>) {
            return obj.repr();
        } else {
            return py::handle(py::repr(convert_to_py(obj))).cast<std::string>();
        }
    }

    #ifdef KDMOL_LOG
        class Logger {
            private:
                static std::string m_log;
                static std::string m_cur_routine;

            public:
                template <class ...Args>
                static void log(Args &&... args) {
                    auto message = fmt::format(args...);
                    std::cout << message << "\n"; 
                    m_log += message + "\n"; 
                }

                static void log_rstart(const std::string & routine_name)
                { log(" START ROUTINE \"{}\"", routine_name); }

                static void log_rend(const std::string & routine_name)
                { log(" END ROUTINE \"{}\"", routine_name); }

                template<class A, class B>
                static void log_item(const std::string & item, const A & key, const B & value)
                { log(" {}:: {} = {}", item, key, value); }

                static void log_routine(const std::string & routine_name)
                { log_rstart(routine_name); m_cur_routine = routine_name; }

                template<class A, class B>
                static void log_in(const A& key, const B& value)
                { log_item("IN", key, repr(value)); }

                template<class A, class B>
                static void log_check(const A& key, const B& value)
                { log_item("CHECK", key, repr(value)); }

                static void log_finish()
                { log_rend(m_cur_routine); }

                static void clear_log()
                { m_log.clear(); aps_count = 0; }

                // Calls from Python
                static void log_py(const py::str& text)
                { log(text.cast<std::string>()); }

                static void log_rstart_py(const py::str& rname)
                { log_rstart(rname.cast<std::string>()); }

                static void log_rend_py(const py::str& rname)
                { log_rend(rname.cast<std::string>()); }

                static const std::string get_log()
                { return m_log; }
                
                static int aps_count;
                
            private:
                Logger() {}
                Logger(Logger const&) = delete;
                void operator=(Logger const&) = delete;
        };

        namespace LoggerWrapper {
            void log_routine(const std::string & routine_name);
            void log_finish();
            template<class A, class B> void log_in(const A& key, const B& value) { Logger::log_in(key, value); }
            template<class A, class B> void log_check(const A& key, const B& value) { Logger::log_check(key, value); }
        }
    #endif

    inline void calc_cross(const double *v_A, const double *v_B, double *c_P) noexcept {
        c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
        c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
        c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
    }

    template<class T>
    void invert_matrix (const T& input, T& inverse) {
        using namespace boost::numeric::ublas;
        T A(input);
        permutation_matrix<std::size_t> pm(A.size1());
        lu_factorize(A, pm);
        inverse.assign(identity_matrix<typename T::value_type>(A.size1()));
        lu_substitute(A, pm, inverse);
    }

    template <typename T, typename Arg> // Yep, it's definitely me who wrote it
    auto at(T&& obj, Arg&& arg) {
        if constexpr (std::is_pointer_v<std::decay_t<T>>) {
            return (*obj).at(std::forward<Arg>(arg));
        } else {
            return obj.at(std::forward<Arg>(arg));
        }
    }

    template <typename T>
    auto& ref(T&& arg) {
        if constexpr (std::is_pointer_v<std::remove_reference_t<T>>) {
            return *arg;
        } else {
            return arg;
        }
    }

    #if defined(VALIDATION) || defined(WEAK_VALIDATION)
        class GeometryValidator {
            protected:
                using vangle_t = std::tuple<int, int, int>;
                std::vector<std::pair<int, int>> bonds;
                std::vector<int> bond_idxs;
                std::vector<vangle_t> vangles;
                std::vector<int> dof_idxs;
                
            public:
                GeometryValidator() { }

                template <class M, class D>
                GeometryValidator(const M& mol_data, const D& dofs) { 
                    for (const auto& atomA : mol_data.xyz.get_atoms())
                        for (const auto& [atomB, bond_data] : mol_data.edge_indexer.atom_data(atomA))
                            if (atomA > atomB) {
                                bonds.push_back(std::make_pair(atomA, atomB));
                                bond_idxs.push_back(bond_data.index);
                            }
                    bonds.shrink_to_fit();
                    bond_idxs.shrink_to_fit();

                    for (const auto& atomB : mol_data.xyz.get_atoms())
                        for (const auto& [atomA, _] : mol_data.edge_indexer.atom_data(atomB))
                            for (const auto& [atomC, __] : mol_data.edge_indexer.atom_data(atomB))
                                if (atomA > atomC)
                                    vangles.push_back(std::make_tuple(atomA, atomB, atomC));
                    vangles.shrink_to_fit();

                    for (int i = 0; i < dofs.bond_container.size(); ++i)
                        dof_idxs.push_back(i);
                    dof_idxs.shrink_to_fit();
                }

                template <class T, class M>
                bool validate_lengths(T& geom, const M& mol_data) {
                    bool okay = true;
                    for (int i = 0; i < bonds.size(); ++i) {
                        if (abs(geom.get_distance(bonds[i].first, bonds[i].second) - 
                                mol_data.bondlengths[bond_idxs[i]]) > THRESHOLD) {
                            okay = false;
                            // fmt::print("Distance {}-{} is {} instead of {} \n", bonds[i].first, bonds[i].second, 
                            //             geom.get_distance(bonds[i].first, bonds[i].second),
                            //             mol_data.bondlengths[bond_idxs[i]]);
                        }
                    }
                    return okay;
                }

                template <class T, class M>
                bool validate_vangles(T& geom, const M& mol_data) {
                    bool okay = true;
                    for (const auto& [atA, atB, atC] : vangles) {
                        if (abs(geom.get_vangle(atA, atB, atC) -
                                mol_data.polys[atB].get_vangle(atA, atC)) > THRESHOLD) {
                            okay = false;
                            // fmt::print("Vangle {}-{}-{} is {} instead of {} \n", atA, atB, atC,
                            //             geom.get_vangle(atA, atB, atC),
                            //             mol_data.polys[atB].get_vangle(atA, atC));
                        }
                    }
                    return okay;
                }

                template <class T, class D>
                bool validate_dihedrals(T& geom, const D& dofs) {
                    bool okay = true;
                    for (const auto& dofidx : dof_idxs) {
                        int atA = dofs.side_container[dofidx].first,
                            atB = dofs.bond_container[dofidx].first,
                            atC = dofs.bond_container[dofidx].second,
                            atD = dofs.side_container[dofidx].second;
                        if (abs(geom.get_dihedral(atA, atB, atC, atD) -
                                dofs.dofvalue_container[dofidx]) > THRESHOLD) {
                            okay = false;
                            // fmt::print("Dihedral {}-{}-{}-{} is {} instead of {} \n", atA, atB, atC, atD,
                            //             geom.get_dihedral(atA, atB, atC, atD),
                            //             dofs.dofvalue_container[dofidx]);
                        }
                    }
                    return okay;
                }

                template <class T, class M, class D>
                bool validate(T& geom, const M& mol_data, const D& dofs) {
                    bool lengths_okay = validate_lengths(geom, mol_data);
                    bool vangles_okay = validate_vangles(geom, mol_data);
                    bool dihedrals_okay = validate_dihedrals(geom, dofs);
                    return lengths_okay && vangles_okay && dihedrals_okay;
                }
        };

        class LocalGeometryValidator : public GeometryValidator {
            public:
                LocalGeometryValidator() { }

                template <class M, class D, class C>
                LocalGeometryValidator(const M& mol_data, const D& dofs, const C& atom_idxs) { 
                    for (const auto& atomA : mol_data.xyz.get_atoms())
                        if (atom_idxs.find(atomA) != atom_idxs.cend())
                            for (const auto& [atomB, bond_data] : mol_data.edge_indexer.atom_data(atomA))
                                if ((atomA > atomB) && (atom_idxs.find(atomB) != atom_idxs.cend())) {
                                    bonds.push_back(std::make_pair(atomA, atomB));
                                    bond_idxs.push_back(bond_data.index);
                                }
                    bonds.shrink_to_fit();
                    bond_idxs.shrink_to_fit();

                    for (const auto& atomB : mol_data.xyz.get_atoms())
                        if (atom_idxs.find(atomB) != atom_idxs.cend())
                            for (const auto& [atomA, _] : mol_data.edge_indexer.atom_data(atomB))
                                if (atom_idxs.find(atomA) != atom_idxs.cend())
                                    for (const auto& [atomC, __] : mol_data.edge_indexer.atom_data(atomB))
                                        if ((atomA > atomC) && ((atom_idxs.find(atomC) != atom_idxs.cend())))
                                            vangles.push_back(std::make_tuple(atomA, atomB, atomC));
                    vangles.shrink_to_fit();

                    for (int i = 0; i < dofs.bond_container.size(); ++i)
                        if ((atom_idxs.find(dofs.bond_container[i].first) != atom_idxs.cend()) &&
                            (atom_idxs.find(dofs.bond_container[i].second) != atom_idxs.cend()) &&
                            (atom_idxs.find(dofs.side_container[i].first) != atom_idxs.cend()) &&
                            (atom_idxs.find(dofs.side_container[i].second) != atom_idxs.cend()))
                            dof_idxs.push_back(i);
                    dof_idxs.shrink_to_fit();
                }
        };
    #endif


    py::list get_status_feed();
    void add_message_to_feed(const std::string& subject, const std::string& message, const std::vector<int>& atoms, const std::string& filename, const int linenumber);
    void add_message_to_feed(const std::string& json_contents);
    std::string assemble_message(const std::string& subject, const std::string& message, const std::vector<int>& atoms, const std::string& filename, const int linenumber);
    void clear_status_feed();
    extern std::unordered_map<std::string, std::string> warning_codes;

    #define prepare_message(subject, message, atoms) Utils::assemble_message(subject, message, atoms, __FILE__, __LINE__)
    #define new_warning(subject, message, atoms) Utils::add_message_to_feed(subject, message, atoms, __FILE__, __LINE__)

    #ifdef BUILD_OVERLAP_DETECTION
        extern std::unordered_map<std::string, double> VDW_RADII;
        extern double RADIUS_MULT;
        py::dict get_vdw_radii();
        void set_vdw_radii(py::dict new_radii);
        void set_radius_multiplier(double rad_mult);

        class OverlapDetector {
            private:
                using check_pairs_t = std::vector<std::vector<int>> ;
                using atoms_t = std::vector<int>;
                using radii_t = std::vector<double>;
                using failed_pairs_t = std::vector<std::pair<int, int>>;
                check_pairs_t check_pairs;
                atoms_t atoms;
                radii_t radii;
                std::string level;

            public:
                OverlapDetector() { }

                template <class C>
                OverlapDetector(const py::object& graph, C& coords, const std::vector<std::string>& atom_symbols, const std::string& level) :
                    level(level)
                {
                    if constexpr(std::is_base_of_v<BasicCoordinates, C>)
                        atoms = coords.get_atoms();
                    else
                        std::copy(coords.cbegin(), coords.cend(), std::back_inserter(atoms));
                    
                    assert(atoms.size() != 0);
                    const auto nx = Utils::PyModules::get_nx();
                    auto number_of_nodes = py::handle(graph.attr("number_of_nodes")()).cast<int>();

                    check_pairs = std::vector<std::vector<int>>(atoms.size());
                    for (int idxA = 0; idxA < atoms.size(); ++idxA) {
                        const auto nodeA = atoms[idxA];
                        for (int idxB = idxA + 1; idxB < atoms.size(); ++idxB) {
                            const auto nodeB = atoms[idxB];
                            auto shortest_path_length = py::handle(nx->attr("shortest_path_length")(graph, py::cast(nodeA), py::cast(nodeB))).cast<int>();
                            if (shortest_path_length >= 3) {
                                check_pairs[idxA].push_back(idxB);
                            }
                        }
                        check_pairs[idxA].shrink_to_fit();
                    }
                    
                    // Protect against unknown element symbols
                    for (int i = 0; i < atoms.size(); ++i) {
                        const auto atom = atoms[i];
                        auto symbol = atom_symbols[atom];
                        auto it = VDW_RADII.find(symbol);
                        if (it == VDW_RADII.end()) {
                            std::vector<int> atoms;
                            atoms.push_back(i);
                            new_warning(Utils::warning_codes["UNKNOWN_ELEMENT"]+"[important]", fmt::format("The element '{}' is unknown. Setting its VdW radius to that of hydrogen atom.", symbol), atoms);
                            VDW_RADII[symbol] = VDW_RADII["H"];
                        }
                    }

                    radii = radii_t(atoms.size());
                    for (int i = 0; i < atoms.size(); ++i) {
                        const auto atom = atoms[i];
                        auto symbol = atom_symbols[atom];
                        auto radius = VDW_RADII.at(symbol) * RADIUS_MULT;
                        radii[i] = radius;
                    }
                }

                template <class C>
                bool operator()(C& coords) {
                    namespace bg = boost::geometry;
                    namespace bgi = boost::geometry::index;
                    using point_t = bg::model::point<double, 3, bg::cs::cartesian>; // 3D point type
                    using rtree_t = bgi::rtree<point_t, bgi::quadratic<16>>; // R-tree index for efficient spatial queries

                    std::vector<point_t> points;
                    for (const auto atom : atoms) {
                        auto cur_xyz = coords.at(atom);
                        points.push_back(point_t(cur_xyz(0), cur_xyz(1), cur_xyz(2)));
                    }
                    assert(points.size() != 0);

                    // Build an R-tree index for the points
                    rtree_t rtree(points);

                    // Iterate over all pairs of points and check if their distance is less than sum of "radii"
                    failed_pairs_t failed_pairs;
                    for (int idxA = 0; idxA < check_pairs.size(); idxA++) {
                        const auto radA = radii[idxA];
                        for (const auto idxB : check_pairs[idxA]) {
                            const auto radB = radii[idxB];
                            double distance = bg::distance(points[idxA], points[idxB]);
                            if (distance < radA + radB) {
                                failed_pairs.push_back({atoms[idxA], atoms[idxB]});
                            }
                        }
                    }

                    bool okay = failed_pairs.size() == 0;
                    checkpoint("overlap_check");
                    #ifdef KDMOL_LOG
                        int log_okay = (okay) ? 1 : 0;
                        
                        py::list log_coords, log_atoms = py::handle(py::cast(atoms)).cast<py::list>();
                        log_atoms.attr("sort")();
                        for (const auto atom : log_atoms) {
                            log_coords.append(vector_to_list(coords.at(py::handle(atom).cast<int>())));
                        }

                        using namespace Utils::LoggerWrapper;
                        log_routine("overlap_check");
                        log_in("atoms", log_atoms);
                        log_in("coords", log_coords);
                        log_in("level", as_list(level));
                        log_check("result", as_list(log_okay));
                        log_check("failed_bonds_SET2", failed_pairs);
                        log_finish();
                    #endif

                    return failed_pairs.size() == 0;
                }
        };
    #endif


    template<typename T>
    bool sets_are_equal(const std::set<T>& set1, const std::set<T>& set2) {
        if (set1.size() != set2.size()) {
            return false;
        }
        for (const auto& elem : set1) {
            if (set2.find(elem) == set2.end()) {
                return false;
            }
        }
        return true;
    }

    template <class T>
    bool vectors_are_equal(const std::vector<T>& a, const std::vector<T>& b) {
        if (a.size() != b.size())
            return false;
        return std::equal(a.cbegin(), a.cend(), b.cbegin());
    }

    template<typename T>
    std::vector<T> index_from_one(const std::vector<T>& vec) {
        std::vector<T> result;
        result.reserve(vec.size());
        for (const auto& elem : vec) {
            result.push_back(elem + 1);
        }
        return result;
    }

    // Custom exceptions
    class TLCZeroSolutions : public std::runtime_error {
    public:
        TLCZeroSolutions(const std::string& msg) : std::runtime_error(msg) {}
    };
    class TLCFailure : public std::runtime_error {
    public:
        TLCFailure(const std::string& msg) : std::runtime_error(msg) {}
    };
    class GeomFailure : public std::runtime_error {
    public:
        GeomFailure(const std::string& msg) : std::runtime_error(msg) {}
    };
    class OverlapDetected : public std::runtime_error {
    public:
        OverlapDetected(const std::string& msg) : std::runtime_error(msg) {}
    };

    // OPT: constexpr VDW_RADII_MAP
    // constexpr std::pair<std::string_view, double> VDW_RADII_ARR[] = { // ... };
    // constexpr std::size_t VDW_RADII_SIZE = sizeof(VDW_RADII_ARR) / sizeof(VDW_RADII_ARR[0]);
    // constexpr std::unordered_map<std::string_view, double> VDW_RADII_MAP(VDW_RADII_ARR, VDW_RADII_ARR + VDW_RADII_SIZE);

    uvector3d_t gs_rand(uvector3d_t& x, uvector3d_t& y);
    std::vector<std::string> readlines(const std::string& filename);
}

enum aps_return_t : int { success = 0, zero = 1, tlcfail = 2, overlap = 3, validationfail = 4 };

#ifdef KDMOL_LOG
    using Utils::Logger;
    using namespace Utils::LoggerWrapper;
    using Utils::as_list;
#endif
using Utils::repr;

using so_matrix_t = Utils::umatrix3d_t;
using se_matrix_t = Utils::umatrix4d_t;

#ifdef BUILD_OVERLAP_DETECTION
    const char* const OVERLAP_MOLECULE_MESSAGE = "Overlap detected at Molecule-level";
    const char* const OVERLAP_GEOMUNIT_MESSAGE = "Overlap detected at GeomUnit-level";
    const char* const OVERLAP_PROBLEM_MESSAGE = "Overlap detected at Problem-level";
#endif

#endif // UTILS_H
