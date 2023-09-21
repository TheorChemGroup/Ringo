#ifndef MOLPROXY_H
#define MOLPROXY_H

class Confpool;

class MolProxy {
    public:
        MolProxy() {}
        MolProxy(Confpool* base, const int& idx);
        py::float_ __getitem__(const py::object& key);
        void __setitem__(const py::object& key, const py::object& value);
        py::object __getattr__(const py::str& py_attr);
        void __setattr__(const py::str& py_attr, const py::object& value);
        
        py::float_ l(const int a, const int b);
        py::float_ v(const int a, const int b, const int c);
        py::float_ z(const int a, const int b, const int c, const int d);

        py::tuple rmsd(const MolProxy& other, const py::kwargs& kwargs);
        py::object get_connectivity();

        int get_index() const { return idx_; }

    private:
        inline const std::string& get_descr();
        inline py::array_t<double> get_xyz();
        
        inline void set_descr(const std::string& new_descr);

        Confpool* base_;
        int idx_;
};


#endif // MOLPROXY_H
