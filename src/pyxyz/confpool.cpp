#include <chrono>
#include <limits>

#ifdef WIN_PLATFORM
    #include <windows.h>
    #include <stdio.h>
    #include <wchar.h>
    #include <string.h>
#endif

#include "confpool.h"
#include "molproxy.h"


py::int_ Confpool::include_from_file(py::str py_filename) {
    int res = include(py_filename.cast<std::string>());
    return res;
}

int Confpool::include(const std::string filename) {
    auto mylines = PyxyzUtils::readlines_rofl(filename);
    const auto filename_str = filename;
    int cline = 0;
    int added_count = 0;
    while (cline < mylines.size()) {
        // std::cout << "Casting to int '" << mylines[cline] << "'" << "\n";
        boost::algorithm::trim(mylines[cline]);
        const unsigned int cur_natoms = boost::lexical_cast<int>(mylines[cline].c_str());
        if (natoms == 0) 
            natoms = cur_natoms;
        else if (natoms != cur_natoms) 
            throw std::runtime_error("Wrong numer of atoms");
        
        auto description = mylines[cline + 1];

        auto geom = coord_container_t(natoms);
        symbol_container_t atom_types;
        for (unsigned int i = 2; i < natoms + 2; ++i)
        {
            if (cline + i >= mylines.size())
                throw std::runtime_error(fmt::format("Unexpected number of atoms (expected {}, nlines={}). Check {}", natoms, mylines.size(), filename_str));
            
            std::vector<std::string> parts;
            boost::algorithm::trim(mylines[cline + i]);
            boost::split(parts, mylines[cline + i], boost::is_any_of(" "), boost::token_compress_on);
            if(parts.size() != 4)
                throw std::runtime_error("Unexpected number of parts in line. Check " + filename_str);
            
            atom_types.push_back(parts[0]);
            geom.set_atom(i - 2, {boost::lexical_cast<double>(parts[1].c_str()),
                                  boost::lexical_cast<double>(parts[2].c_str()),
                                  boost::lexical_cast<double>(parts[3].c_str())});
        }
        
        if (sym_.size() == 0)
            sym_ = atom_types;
        else if (sym_ != atom_types)
            throw std::runtime_error("Unexpected atom types. Check " + filename_str);
        
        coord_.push_back(geom);
        descr_.push_back(description);
        cline += 2 + natoms;
        added_count++;
    }
    resize();
    return added_count;
}

py::object Confpool::filter(const py::function& py_criterion, const py::bool_& inplace) {
    full_check();

    if (!inplace) {
        py::dict res;
        auto obj = Confpool(*this);
        res["DelCount"] = obj.filter(py_criterion, true);
        res["Object"] = obj;
        return res;
    }

    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (!py_criterion(proxies_[i]).cast<bool>()) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return py::cast(del_count);
}

py::int_ Confpool::count(const py::function& py_criterion) {
    full_check();

    unsigned int res = 0;
    for(int i = coord_.size() - 1; i >= 0; --i)
        if (py_criterion(proxies_[i]).cast<bool>())
            res += 1;
    return res;
}

void Confpool::remove_key(const std::string& keyname) {
    keys_.erase(keyname);
}

void Confpool::remove_structure(const int& i) {
    auto c_it = coord_.begin();
    std::advance(c_it, i);
    coord_.erase(c_it);

    auto d_it = descr_.begin();
    std::advance(d_it, i);
    descr_.erase(d_it);

    for(auto& pair : keys_) { // pair = key, std::vector<double>
        auto it = pair.second.begin();
        std::advance(it, i);
        pair.second.erase(it);
    }
}

void Confpool::resize() {
    if (coord_.size() != descr_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and descr): {} vs. {}", coord_.size(), descr_.size()));
    
    for(auto& pair : keys_) { // pair = key, std::vector<double>
        if (pair.second.size() > coord_.size())
            throw std::runtime_error(fmt::format("Mismatch of container sizes (key container > coord): {} < {}", pair.second.size(), coord_.size()));
        else if (pair.second.size() < coord_.size())
            pair.second.resize(coord_.size());
    }
    
    if (proxies_.size() != coord_.size()) {
        proxies_.resize(coord_.size());
        for (int i = 0; i < coord_.size(); ++i)
            proxies_[i] = MolProxy(this, i);
    }
}

py::object Confpool::upper_cutoff(const py::str& py_keyname, const py::float_& py_cutoff, const py::bool_& inplace) {
    full_check();

    if (!inplace) {
        py::dict res;
        auto obj = Confpool(*this);
        res["DelCount"] = obj.upper_cutoff(py_keyname, py_cutoff, true);
        res["Object"] = obj;
        return res;
    }

    const auto cutoff = py_cutoff.cast<double>();
    if (cutoff <= 0.0)
        throw std::runtime_error(fmt::format("Cutoff value must be > 0. {} given.", cutoff));

    const auto& key_data = keys_[py_keyname.cast<std::string>()];

    auto minimal_value = key_data[0];
    for (const auto& current_val : key_data)
        if (current_val < minimal_value)
            minimal_value = current_val;
    
    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (key_data[i] - minimal_value > cutoff) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return py::cast(del_count);
}

py::object Confpool::lower_cutoff(const py::str& py_keyname, const py::float_& py_cutoff, const py::bool_& inplace) {
    full_check();

    if (!inplace) {
        py::dict res;
        auto obj = Confpool(*this);
        res["DelCount"] = obj.lower_cutoff(py_keyname, py_cutoff, true);
        res["Object"] = obj;
        return res;
    }

    const auto cutoff = py_cutoff.cast<double>();
    if (cutoff <= 0.0)
        throw std::runtime_error(fmt::format("Cutoff value must be > 0. {} given.", cutoff));

    const auto& key_data = keys_[py_keyname.cast<std::string>()];

    auto maximal_value = key_data[0];
    for (const auto& current_val : key_data)
        if (maximal_value < current_val)
            maximal_value = current_val;
    
    unsigned int del_count = 0;
    for(int i = coord_.size() - 1; i >= 0; --i) {
        if (maximal_value - key_data[i] > cutoff) {
            remove_structure(i);
            del_count += 1;
        }
    }
    resize();
    return py::cast(del_count);
}

void Confpool::save(const py::str& py_filename) const {
    full_check();

    auto natoms_str = boost::lexical_cast<std::string>(natoms);
    std::vector<std::string> reslines;
    for(auto i = 0; i < coord_.size(); ++i) {
        reslines.push_back(natoms_str);
        reslines.push_back(descr_[i]);

        for(auto j = 0; j < natoms; ++j) {
            const auto* coords = coord_[i].get_atom_raw(j);
            reslines.push_back(fmt::format("{:>2}  {:12.8f}  {:12.8f}  {:12.8f}", sym_[j], coords[0], coords[1], coords[2]));
        }
    }

    auto joined = boost::algorithm::join(reslines, "\n");
    #ifdef WIN_PLATFORM
        std::ofstream out(PyxyzUtils::ConvertFromUtf16ToUtf8(py_filename.cast<std::wstring>())); //filename.c_str(), std::ios::in | std::ios::binary
    #else
        std::ofstream out(py_filename.cast<std::string>());
    #endif
    out << joined << "\n";
    out.close();
}

py::dict Confpool::as_table() const {
    py::dict res;
    for(const auto& pair : keys_) { // pair = key, std::vector<double>
        res[py::cast(pair.first)] = py::list(py::cast(pair.second));
    }
    return res;
}

void Confpool::full_check() const {
    if (coord_.size() != descr_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and descr): {} vs. {}", coord_.size(), descr_.size()));
    
    for(auto& pair : keys_) { // pair = key, std::vector<double>
        if (pair.second.size() != coord_.size())
            throw std::runtime_error(fmt::format("Mismatch of container sizes (key='{}' container != coord): {} != {}", pair.first, pair.second.size(), coord_.size()));
    }

    if (coord_.size() != proxies_.size())
        throw std::runtime_error(fmt::format("Mismatch of container sizes (coord and proxy container): {} vs. {}", coord_.size(), proxies_.size()));
    
    for (int i = 0; i < coord_.size(); ++i)
        if (proxies_[i].get_index() != i)
            throw std::runtime_error(fmt::format("MolProxy #{} has number {} (must be equal)", i, proxies_[i].get_index()));
}

py::object Confpool::sort(const py::str& py_keyname, const py::bool_& py_ascend, const py::bool_& inplace) {
    full_check();

    if (!inplace) {
        auto obj = Confpool(*this);
        obj.sort(py_keyname, py_ascend, true);
        return py::cast(obj);
    }

    // bool ascend = true;
    // if (kwargs.attr("__contains__")("ascending").cast<bool>())
    //     ascend = kwargs["ascending"].cast<bool>();
    bool ascend = py_ascend.cast<bool>();

    const auto keyname = py_keyname.cast<std::string>();
    auto p = PyxyzUtils::sort_permutation(keys_[keyname], ascend);

    PyxyzUtils::apply_permutation_in_place(descr_, p);
    PyxyzUtils::apply_permutation_in_place(coord_, p);
    for(auto& pair : keys_) // pair = key, std::vector<double>
        PyxyzUtils::apply_permutation_in_place(pair.second, p);
    return py::none();
}

py::dict Confpool::rmsd_filter(const py::float_& py_rmsd_cutoff, const py::kwargs& kwargs) {
    full_check();

    bool inplace = true;
    if (kwargs.attr("__contains__")("inplace").cast<bool>())
        inplace = kwargs["inplace"].cast<bool>();
    if (!inplace) {
        kwargs["inplace"] = true;
        auto obj = Confpool(*this);
        py::dict res = obj.rmsd_filter(py_rmsd_cutoff, kwargs);
        res["Object"] = py::cast(obj);
        return res;
    }

    double energy_threshold = 0.0;
    double* energies;
    bool use_energy_threshold = false;
    if (kwargs.attr("__contains__")("energy_threshold").cast<bool>()) {
        use_energy_threshold = true;
        energy_threshold = kwargs["energy_threshold"].cast<double>();
        auto energy_key = kwargs["energy_key"].cast<std::string>();
        energies = &(keys_[energy_key][0]);
    }

    bool mirror_match = true;
    if (kwargs.attr("__contains__")("mirror_match").cast<bool>()) {
        mirror_match = kwargs["mirror_match"].cast<bool>();
    }
    
    bool print_status = false;
    if (kwargs.attr("__contains__")("print_status").cast<bool>()) {
        print_status = kwargs["print_status"].cast<bool>();
    }

    bool use_rmsd_matrix = false;
    py::buffer_info rmsd_buffer;
    py::array_t<double> rmsd_matrix;
    double* rmsd_matrix_raw = nullptr;
    if (kwargs.attr("__contains__")("rmsd_matrix").cast<bool>()) {
        use_rmsd_matrix = true;
        rmsd_matrix = kwargs["rmsd_matrix"].cast<py::array_t<double>>();
        rmsd_buffer = rmsd_matrix.request();
        if ((rmsd_matrix.shape(0) != coord_.size()) || (rmsd_matrix.shape(1) != coord_.size())) {
            throw std::runtime_error(fmt::format("Mismatch between dimensions of RMSD matrix {}x{} and Nstructures={}", 
                                rmsd_matrix.shape(0), rmsd_matrix.shape(1), coord_.size()));
        }
        rmsd_matrix_raw = static_cast<double *>(rmsd_buffer.ptr);
    }
    
    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();
    unsigned int del_count = 0;
    const auto cutoff = py_rmsd_cutoff.cast<double>();
    
    if (isomorphisms.size() == 0)
        throw std::runtime_error("No isomorphisms were generated yet. Run generate_connectivity and generate_isomorphisms");
    RmsdCalculator<isomorphism_container_t> rmsd(isomorphisms, simple_reorder, cutoff, mirror_match);
    double minimal_rmsd = std::numeric_limits<double>::max();
    int min_pairA = -1, min_pairB = -1;
    int rmsd_calc_count = 0, skip_count = 0;
    int start_size = coord_.size(),
        full_set = start_size * (start_size - 1) / 2;
    
    std::vector<int> struct_idxs;
    if (use_rmsd_matrix) {
        for (int i = 0; i < coord_.size(); ++i)
            struct_idxs.push_back(i);
        struct_idxs.shrink_to_fit();
    }

    for (int i = 0; i < coord_.size(); i++) {
        if (print_status) {
            double progress = (i + del_count - 1) * (i + del_count) * 100;
            progress /= 2 * full_set;
            py_print(fmt::format("Progress {:5.2f}%, step = {}/{}", progress, i + del_count, start_size));
        }
        double cur_min_rmsd = minimal_rmsd;
        int cur_min_pair = -1;
        auto& curgeom = coord_[i].to_boost_format();
        bool removed = false;
        for (int j = i - 1; j >= 0; j--) {
            if (use_energy_threshold && (abs(energies[i] - energies[j]) > energy_threshold)) {
                skip_count++;
                continue;
            }

            double rmsd_value;
            if (use_rmsd_matrix) {
                rmsd_value = rmsd_matrix.at(struct_idxs[i], struct_idxs[j]);
            } else {
                auto& testgeom = coord_[j].to_boost_format();
                rmsd_value = rmsd.calc(curgeom, testgeom);
            }
            rmsd_calc_count++;
            if (rmsd_value < cutoff) {
                remove_structure(i);
                
                if (use_rmsd_matrix) {
                    auto c_it = struct_idxs.begin();
                    std::advance(c_it, i);
                    struct_idxs.erase(c_it);
                }

                del_count++;
                removed = true;
                i--;
                break;
            } else if (rmsd_value < cur_min_rmsd) {
                cur_min_pair = j;
                cur_min_rmsd = rmsd_value;
            }
        }

        if (!removed && (cur_min_rmsd < minimal_rmsd)) {
            if (cur_min_pair == -1)
                throw std::runtime_error("A bug detected. 306");
            min_pairA = i;
            min_pairB = cur_min_pair;
            minimal_rmsd = cur_min_rmsd;
        }
    }
    auto stop_time = high_resolution_clock::now();
    duration<double, std::milli> fp_ms = stop_time - start_time;

    if ((min_pairA == -1) && (coord_.size() > 1))
        throw std::runtime_error("A bug detected. 314");
    if ((min_pairB == -1) && (coord_.size() > 1))
        throw std::runtime_error("A bug detected. 316");
    py::dict res;
    res["DelCount"] = del_count;
    res["MinRMSD_pairA"] = min_pairA;
    res["MinRMSD_pairB"] = min_pairB;
    res["MinRMSD"] = minimal_rmsd;
    res["TimeElapsed"] = fp_ms.count();
    res["NRMSDCalcs"] = rmsd_calc_count;
    res["NSkipped"] = skip_count;
    resize();
    return res;
}

py::object Confpool::__getitem__(const py::object& key) {
    if (py::isinstance<py::int_>(key)) {
        if (key.cast<int>() >= coord_.size())
            throw py::index_error(fmt::format("Index {} is invalid. Size of the pool is {}", key.cast<int>(), coord_.size()));
        return py::cast(proxies_[key.cast<int>()]);
    } else if (py::isinstance<py::str>(key))
        return key_to_list(key.cast<std::string>());
    else
        throw std::runtime_error(fmt::format("Expected either an integer (conformer index) or a string (a key). Got a {}", py::repr(key).cast<std::string>()));
}

void Confpool::__setitem__(const py::object& key, const py::object& func) {
    if (!py::isinstance<py::str>(key))
            throw std::runtime_error(fmt::format("Expected str as key. Got a {}", py::repr(key).cast<std::string>()));
    if (!py::isinstance<py::function>(func))
            throw std::runtime_error(fmt::format("Expected a function for keyvalue modification. Got a {}", py::repr(func).cast<std::string>()));
    modify_key(key.cast<std::string>(), func.cast<py::function>());
}

void Confpool::modify_key(const std::string& keyname, const py::function& func) {
    full_check();
    prepare_key(keyname);
    for(int i = 0; i < coord_.size(); ++i)
        keys_[keyname][i] = func(proxies_[i]).cast<double>();
}

void Confpool::modify_descr(const py::function& func) {
    full_check();
    for(int i = 0; i < coord_.size(); ++i)
        descr_[i] = func(proxies_[i]).cast<std::string>();
}

py::object Confpool::__getattr__(const py::str& py_attr) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "size")
        return py::cast(coord_.size());
    else if (attr == "atom_symbols")
        return get_atom_symbols();
    else
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
}

void Confpool::__setattr__(const py::str& py_attr, const py::object& value) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "descr") {
        if (!py::isinstance<py::function>(value))
            throw std::runtime_error(fmt::format("Expected a function for description modification. Got a {}", py::repr(value).cast<std::string>()));
        modify_descr(value.cast<py::function>());
    } else if (attr == "atom_symbols") {
        if (!py::isinstance<py::list>(value))
            throw std::runtime_error(fmt::format("Expected a list of strs. Got a {}", py::repr(value).cast<std::string>()));
        sym_ = value.cast<symbol_container_t>();
    } else {
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
    }
}

void Confpool::__delitem__(const py::object& py_key) {
    if (py::isinstance<py::int_>(py_key)) {
        remove_structure(py_key.cast<int>());
        resize();
    } else if (py::isinstance<py::str>(py_key)) {
        remove_key(py_key.cast<std::string>());
    } else {
        throw std::runtime_error(fmt::format("Expected either keyname or structure index for deletion. Got {}", py::repr(py_key).cast<std::string>()));
    }
    full_check();
}

inline py::list Confpool::key_to_list(const std::string& key) const {
    py::list res;
    for (const auto& item : keys_.at(key))
        res.append(item);
    return res;
}

void Confpool::generate_connectivity(const int& index, const py::kwargs& kwargs) {
    rmsd_graph = mt_module.attr("generate_connectivity")(this, index, kwargs).cast<py::object>();
}

void Confpool::generate_connectivity_from_graph(py::object molgr, const py::kwargs& kwargs) {
    rmsd_graph = mt_module.attr("generate_connectivity_from_graph")(this, molgr, kwargs).cast<py::object>();
}

py::int_ Confpool::generate_isomorphisms() {
    isomorphisms = mt_module.attr("generate_isomorphisms")(rmsd_graph).cast<isomorphism_container_t>();
    simple_reorder = py::list(rmsd_graph.attr("nodes")()).cast<single_isomorphism_t>();
    return isomorphisms.size(); // Number of isomorphisms
}

py::array_t<double> Confpool::get_rmsd_matrix(const py::kwargs& kwargs) {
    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();

    int nstructs = coord_.size();
    const size_t size = nstructs * nstructs;
    double *coord_array = new double[size];
    
    bool mirror_match = true;
    if (kwargs.attr("__contains__")("mirror_match").cast<bool>()) {
        mirror_match = kwargs["mirror_match"].cast<bool>();
    }

    bool print_status = false;
    if (kwargs.attr("__contains__")("print_status").cast<bool>()) {
        print_status = kwargs["print_status"].cast<bool>();
    }

    bool print_stats = false;
    if (kwargs.attr("__contains__")("print_stats").cast<bool>()) {
        print_stats = kwargs["print_stats"].cast<bool>();
    }

    RmsdCalculator<isomorphism_container_t> rmsd(isomorphisms, simple_reorder, -1.0, mirror_match);
    int rmsd_calc_count = 0;
    int start_size = coord_.size(),
        full_set = start_size * (start_size - 1) / 2;
    for (int i = 0; i < coord_.size(); i++) {
        if (print_status) {
            double progress = (i - 1) * i * 100;
            progress /= 2 * full_set;
            py_print(fmt::format("Progress {:5.2f}%, step = {}/{}", progress, i, start_size));
        }
        auto& curgeom = coord_[i].to_boost_format();
        for (int j = i; j >= 0; j--) {
            auto& testgeom = coord_[j].to_boost_format();
            auto rmsd_value = rmsd.calc(curgeom, testgeom);
            coord_array[i * nstructs + j] = rmsd_value;
            coord_array[j * nstructs + i] = rmsd_value;
            rmsd_calc_count++;
        }
    }
    auto stop_time = high_resolution_clock::now();
    duration<double, std::milli> fp_ms = stop_time - start_time;
    if (print_stats)
        py_print(fmt::format("RMSD matrix statistics:\nTime elapsed = {} ms\nPer RMSD calc = {} ms", fp_ms.count(), fp_ms.count() / rmsd_calc_count));
    
    py::capsule free_when_done(coord_array, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        delete[] foo;
    });
    
    return py::array_t<double>(
        {nstructs, nstructs}, // shape
        {nstructs*8, 8}, // C-style contiguous strides for double
        coord_array, // the data pointer
        free_when_done);
}

void Confpool::include_from_xyz_py(const py::array_t<double>& xyz, const py::str& descr) {
    py::buffer_info input_buf = xyz.request();
    if (input_buf.ndim != 2)
        throw std::runtime_error("numpy.ndarray dims must be 2!");

    if (natoms == 0) {
        natoms = input_buf.shape[0];
        sym_ = symbol_container_t(natoms);
        py_print(fmt::format("The object is initialized with {} atoms and empty element symbols.", natoms));
    }

    if ((input_buf.shape[0] != natoms) || (input_buf.shape[1] != 3))
        throw std::runtime_error(fmt::format("Expected dimensions ({}, {}). Got ({}, {})", natoms, 3, input_buf.shape[0], input_buf.shape[1]));

    full_check();
    auto geom = coord_container_t(natoms);
    double* ptr1 = (double*)input_buf.ptr;
    for (int i = 0; i < input_buf.shape[0]; i++)
    {
        geom.set_atom(i, {ptr1[i * input_buf.shape[1]],
                          ptr1[i * input_buf.shape[1] + 1],
                          ptr1[i * input_buf.shape[1] + 2]});
    }
    coord_.push_back(geom);
    descr_.push_back(descr.cast<std::string>());
    resize();
}

void Confpool::set_atom_symbols(const symbol_container_t& symbols) {
    if (natoms == 0)
        natoms = symbols.size();
    else
        assertm(natoms == symbols.size(), "Mismatch between natoms and symbols.size()");
    sym_ = symbols;
}

void Confpool::include_from_xyz(const coord_container_t& xyz, const std::string& descr) {
    const auto& raw_boost_matrix = xyz.to_boost_format();
    if (natoms == 0) { // TODO !!! SKIP THESE CHECKS !!!
        natoms = raw_boost_matrix.size1();
        sym_ = symbol_container_t(natoms);
        py_print(fmt::format("The object is initialized with {} atoms and empty element symbols.", natoms));
    }

    if ((raw_boost_matrix.size1() != natoms) || (raw_boost_matrix.size2() != 3))
        throw std::runtime_error(fmt::format("Expected dimensions ({}, {}). Got ({}, {})", natoms, 3, raw_boost_matrix.size1(), raw_boost_matrix.size2()));

    full_check();
    coord_.push_back(xyz);
    PyxyzUtils::make_centered(coord_[coord_.size() - 1].to_boost_format());
    descr_.push_back(descr);
    resize();
}

void Confpool::add_rmsd_settings(const RmsdSettings rmsd_settings) {
    if (isomorphisms.size() == 0)
        throw std::runtime_error("No isomorphisms were generated yet. Run generate_connectivity and generate_isomorphisms");
    rmsd_calculator = rmsd_calculator_t(isomorphisms, simple_reorder, rmsd_settings.rmsd_cutoff, rmsd_settings.mirror_match);
}

bool Confpool::include_from_xyz_with_rmsdcheck(const coord_container_t& xyz, const std::string& descr) {
    const auto& raw_boost_matrix = xyz.to_boost_format();
    if (natoms == 0) {
        natoms = raw_boost_matrix.size1();
        sym_ = symbol_container_t(natoms);
        py_print(fmt::format("The object is initialized with {} atoms and empty element symbols.", natoms));
    }

    if ((raw_boost_matrix.size1() != natoms) || (raw_boost_matrix.size2() != 3))
        throw std::runtime_error(fmt::format("Expected dimensions ({}, {}). Got ({}, {})", natoms, 3, raw_boost_matrix.size1(), raw_boost_matrix.size2()));

    if (isomorphisms.size() == 0)
        throw std::runtime_error("No isomorphisms were generated yet. Run generate_connectivity and generate_isomorphisms");
    
    if (rmsd_calculator.get_natoms() == -1)
        throw std::runtime_error("Settings for confsearch were not provided.");

    for (int i = 0; i < coord_.size(); i++) {
        auto& curgeom = coord_[i].to_boost_format();

        double rmsd_value = rmsd_calculator.calc(curgeom, raw_boost_matrix);
        if (rmsd_value < rmsd_calculator.get_cutoff())
            return false;
    }

    full_check();
    coord_.push_back(xyz);
    PyxyzUtils::make_centered(coord_[coord_.size() - 1].to_boost_format());
    descr_.push_back(descr);
    resize();
    return true;
}

void Confpool::include_subset(Confpool& other, const py::list& py_idxs) {
    full_check();
    auto idxs = py_idxs.cast<std::vector<int>>();
    if (natoms == 0) {
        natoms = other.get_natoms();
        sym_ = other.symbols_access();
    }

    for(const auto& pair : other.key_full_access())
        prepare_key(pair.first);

    for (const auto& idx : idxs) {
        coord_.push_back(other.coord_access(idx));
        descr_.push_back(other.descr_access(idx));
    }
    
    for (const auto& idx : idxs)
        for(const auto& pair : other.key_full_access())
            keys_[pair.first].push_back(pair.second[idx]);
    resize();
}

void Confpool::float_from_descr(const py::str& keyname, const py::int_& py_idx) {
    full_check();
    prepare_key(keyname);
    auto idx = py_idx.cast<int>() - 1;
    if (idx + 1 < 0)
        throw std::runtime_error(fmt::format("Index must be 1 or greater. {} was given.", idx + 1));

    std::string float_regex = R"###([-+]?(?:\d*\.\d+|\d+))###";
    std::string int_regex = R"(^\d+$)";
    for(int i = 0; i < coord_.size(); ++i) {
        int count = 0;
        double myvalue;
        py::list float_matches = utils_module.attr("get_matches")(float_regex, descr_[i]);
        for (const auto& float_part : float_matches) {
            if (utils_module.attr("get_matches")(int_regex, float_part).cast<py::list>().size() == 0) {
                if (count == idx) {
                    myvalue = boost::lexical_cast<double>(float_part.cast<std::string>());
                    break;
                }
                count += 1;
            }
        }
        keys_[keyname][i] = myvalue;
    }
}

py::tuple Confpool::calc_rmsd(const int& idxA, const int& idxB, const py::kwargs& kwargs) {
    bool mirror_match = true;
    if (kwargs.attr("__contains__")("mirror_match").cast<bool>()) {
        mirror_match = kwargs["mirror_match"].cast<bool>();
    }

    if (isomorphisms.size() == 0)
        throw std::runtime_error("No isomorphisms were generated yet. Run generate_connectivity and generate_isomorphisms");
    
    RmsdCalculator<isomorphism_container_t> rmsd(isomorphisms, simple_reorder, -1.0, mirror_match);
    auto& curgeom = coord_[idxA].to_boost_format();
    auto& testgeom = coord_[idxB].to_boost_format();
    return rmsd.manual_calc(curgeom, testgeom);
}

Confpool& Confpool::__iter__() {
    iteration_idx = 0;
    return *this;
}

MolProxy Confpool::__next__() {
    if (iteration_idx < coord_.size()) {
        auto& res = proxies_[iteration_idx];
        iteration_idx++;
        return res;
    } else {
        throw py::stop_iteration();
    }
}
