#ifndef CONFPOOL_H
#define CONFPOOL_H

#include "utils.h"
#include "molproxy.h"
#include "rmsd.h"

struct RmsdSettings {
    double rmsd_cutoff;
    bool mirror_match;
};

class Confpool {
    private:
        using coord_container_t = Utils::BoostXyzContainer<3>;
        using symbol_container_t = std::vector<std::string>;
        using single_isomorphism_t = std::vector<unsigned int>;
        using isomorphism_container_t = std::vector<single_isomorphism_t>;
        using rmsd_calculator_t = RmsdCalculator<isomorphism_container_t>;
 
        symbol_container_t sym_;
        std::vector<coord_container_t> coord_;
        std::vector<std::string> descr_;
        std::unordered_map<std::string, std::vector<double>> keys_;
        std::vector<MolProxy> proxies_;
        isomorphism_container_t isomorphisms;
        single_isomorphism_t simple_reorder;
        unsigned int natoms;
        unsigned int iteration_idx;
        py::object rmsd_graph;
        rmsd_calculator_t rmsd_calculator;
        py::module_ mt_module, utils_module;

        int include(const std::string filename);
        void resize();
        void remove_structure(const int& idx);
        void remove_key(const std::string& keyname);
        void full_check() const;
        void modify_key(const std::string& key_name, const py::function& func);
        void modify_descr(const py::function& func);
        inline py::list key_to_list(const std::string& key) const;
        inline py::list get_atom_symbols() const noexcept { return py::cast(sym_); }

    public:
        Confpool(const Confpool &) = default;
        Confpool(const py::object& copy) {
            if (copy.is(py::none())) {
                natoms = 0;
                rmsd_graph = py::none();
                mt_module = py::module::import("ringo.pyutils.pyxyz.moltopology");
                utils_module = py::module::import("ringo.pyutils.pyxyz.utils");
            } else {
                *this = Confpool(copy.cast<const Confpool&>());
            }
        }
        Confpool& __iter__();
        MolProxy __next__();
        py::object __getitem__(const py::object& key);
        void __setitem__(const py::object& key, const py::object& func);
        py::object __getattr__(const py::str& py_attr);
        void __setattr__(const py::str& py_attr, const py::object& func);
        void __delitem__(const py::object& py_key);
        py::int_ __len__() { return py::cast(coord_.size()); }

        py::int_ include_from_file(py::str py_filename);
        void include_from_xyz_py(const py::array_t<double>& xyz, const py::str& descr);
        void include_from_xyz(const coord_container_t& xyz, const std::string& descr);
        bool include_from_xyz_with_rmsdcheck(const coord_container_t& xyz, const std::string& descr);
        void include_subset(Confpool& other, const py::list& py_idxs);
        void float_from_descr(const py::str& keyname, const py::int_& idx);
        py::object filter(const py::function& py_parser, const py::bool_& inplace);
        py::object upper_cutoff(const py::str& py_keyname, const py::float_& py_cutoff, const py::bool_& inplace);
        py::object lower_cutoff(const py::str& py_keyname, const py::float_& py_cutoff, const py::bool_& inplace);
        py::dict rmsd_filter(const py::float_& py_rmsd_cutoff, const py::kwargs& kwargs);

        py::int_ count(const py::function& py_criterion);

        void generate_connectivity(const int& index, const py::kwargs& kwargs);
        void generate_connectivity_from_graph(py::object molgr, const py::kwargs& kwargs);
        py::int_ generate_isomorphisms();
        void add_rmsd_settings(const RmsdSettings rmsd_settings);

        py::object sort(const py::str& py_keyname, const py::bool_& py_ascend, const py::bool_& inplace);
        void save(const py::str& py_filename) const;
        py::dict as_table() const;
        py::array_t<double> get_rmsd_matrix(const py::kwargs& kwargs);
        py::tuple calc_rmsd(const int& idxA, const int& idxB, const py::kwargs& kwargs);

        inline int size() const noexcept { return coord_.size(); }
        inline const symbol_container_t& symbols_access() const { return sym_; }
        void set_atom_symbols(const symbol_container_t& symbols);
        inline coord_container_t& coord_access(const int& i) { return coord_[i]; }
        inline std::string& descr_access(const int& i) { return descr_[i]; }
        inline std::vector<double>& key_access(const std::string& key) { return keys_.at(key); }
        inline const std::unordered_map<std::string, std::vector<double>>& key_full_access() const { return keys_; }
        inline int get_natoms() const noexcept { return natoms; }
        py::object get_connectivity() { return this->rmsd_graph; }
        inline int get_num_isomorphisms() const noexcept { return isomorphisms.size(); }
        
        inline void prepare_key(const std::string& keyname) noexcept {
            if (keys_.find(keyname) == keys_.end())
                keys_[keyname] = std::vector<double>(coord_.size());
        }

        inline void py_print(const std::string& text) noexcept {
            utils_module.attr("py_print")(py::cast(text));
        }
};


class SharedConfpool {
    public:
        SharedConfpool (Confpool& p) : p(p), conformers_count(0) { }
        
        int size() const {
            return conformers_count.load();
        }

        template <class C>
        void include_from_xyz(const C& new_coords, const std::string& description) {
            std::lock_guard<std::mutex> guard(mutex);
            p.include_from_xyz(new_coords, description);
            conformers_count.fetch_add(1);
        }
        
        template <class C>
        bool include_from_xyz_with_rmsdcheck(const C& new_coords, const std::string& description) {
            std::lock_guard<std::mutex> guard(mutex);
            auto successful = p.include_from_xyz_with_rmsdcheck(new_coords, description);
            if (successful)
                conformers_count.fetch_add(1);
            return successful;
        }

    private:
        Confpool& p;
        std::mutex mutex;
        std::atomic<int> conformers_count;
};

#endif // CONFPOOL_H