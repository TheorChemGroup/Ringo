#include "confpool.h"
#include "molproxy.h"


MolProxy::MolProxy(Confpool* base, const int& idx) :
    base_(base), idx_(idx)
{}

py::float_ MolProxy::__getitem__(const py::object& key) {
    if (!py::isinstance<py::str>(key))
        throw std::runtime_error(fmt::format("Expected string as a key. Got {}", py::repr(key).cast<std::string>()));
    return base_->key_access(key.cast<std::string>())[idx_];
}

void MolProxy::__setitem__(const py::object& key, const py::object& value) {
    if (!py::isinstance<py::str>(key))
        throw std::runtime_error(fmt::format("Expected string as a key. Got {}", py::repr(key).cast<std::string>()));
    if (!py::isinstance<py::float_>(value))
        throw std::runtime_error(fmt::format("Keys must refer to floats. Got {}", py::repr(value).cast<std::string>()));
    
    const auto keyname = key.cast<std::string>();
    base_->prepare_key(keyname);
    base_->key_access(keyname)[idx_] = value.cast<double>();
}

py::float_ MolProxy::l(const int a, const int b)
{ return base_->coord_access(idx_).get_distance(a - 1, b - 1); }

py::float_ MolProxy::v(const int a, const int b, const int c)
{ return base_->coord_access(idx_).get_vangle(a - 1, b - 1, c - 1); }

py::float_ MolProxy::z(const int a, const int b, const int c, const int d)
{ return base_->coord_access(idx_).get_dihedral(a - 1, b - 1, c - 1, d - 1); }

py::tuple MolProxy::rmsd(const MolProxy& other, const py::kwargs& kwargs) {
    if (this->base_ != other.base_)
        throw std::runtime_error("Structures appear to be from different Confpool objects");
    return this->base_->calc_rmsd(this->idx_, other.idx_, kwargs);
}

inline py::array_t<double> MolProxy::get_xyz() {
    const size_t size = base_->get_natoms() * 3;
    double *coord_array = new double[size];
    const auto& coords = base_->coord_access(idx_).to_boost_format();
    for (size_t i = 0; i < base_->get_natoms(); i++)
        for (size_t j = 0; j < 3; j++)
            coord_array[i * 3 + j] = coords(i, j);

    py::capsule free_when_done(coord_array, [](void *f) {
        double *foo = reinterpret_cast<double *>(f);
        delete[] foo;
    });

    return py::array_t<double>(
        {static_cast<int>(base_->get_natoms()), 3}, // shape
        {3*8, 8}, // C-style contiguous strides for double
        coord_array, // the data pointer
        free_when_done);
}

inline const std::string& MolProxy::get_descr()
{ return base_->descr_access(idx_); }

inline void MolProxy::set_descr(const std::string& new_descr)
{ base_->descr_access(idx_) = new_descr; }

py::object MolProxy::__getattr__(const py::str& py_attr) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "xyz")
        return get_xyz();
    else if (attr == "descr")
        return py::cast(get_descr());
    else if (attr == "idx")
        return py::cast(get_index());
    else
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
}

void MolProxy::__setattr__(const py::str& py_attr, const py::object& value) {
    const auto attr = py_attr.cast<std::string>();
    if (attr == "descr") {
        if (!py::isinstance<py::str>(value))
            throw std::runtime_error(fmt::format("Description is expected to be a string. Got {}", py::repr(value).cast<std::string>()));
        set_descr(value.cast<std::string>());
    } else {
        throw std::runtime_error(fmt::format("Unknown attr {}", attr));
    }
}

py::object MolProxy::get_connectivity() { return base_->get_connectivity(); }
