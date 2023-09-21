#include "utils.h"

namespace TlcUtils {
    AtomItem AtomItem::operator++(int) noexcept { 
        *this = this->right_bond->get_right_atom();
        return *this;
    }
    
    AtomItem AtomItem::operator--(int) noexcept {
        *this = this->left_bond->get_left_atom();
        return *this;
    }

    AtomItem AtomItem::operator+(const int& step) const noexcept {
        int i = 0;
        AtomItem x = *this;
        while(i++ < step)
            x++;
        return x;
    }

    AtomItem& AtomItem::operator-=(const int& step) noexcept {
        *this = *this + step;
        return *this;
    }

    AtomItem AtomItem::operator-(const int& step) const noexcept {
        int i = 0;
        AtomItem x = *this;
        while(i++ < step)
            x--;
        return x;
    }

    AtomItem& AtomItem::operator+=(const int& step) noexcept {
        *this = *this - step;
        return *this;
    }

    BondItem BondItem::operator++(int) noexcept {
        *this = this->right_atom->get_right_bond();
        return *this;
    }
    
    BondItem BondItem::operator--(int) noexcept {
        *this = this->left_atom->get_left_bond();
        return *this;
    }

    AtomContainerDynamic AtomGroupSimple::collect_proper_atoms() noexcept {
        AtomContainerDynamic res;
        for(auto& atom : items)
            if (atom.can_be_axial())
                res += atom;
        return res;
    }

    void connect(AtomGroupSimple& atoms, BondGroupSimple& bonds) noexcept {
        int n_atoms = atoms.size();
        for(int i = 0; i < n_atoms; ++i) {
            atoms[i].set_bonds(&bonds[i - 1], &bonds[i]);
            bonds[i].set_atoms(&atoms[i], &atoms[i + 1]);
        }
    }

    #ifdef EXT_LOG
        std::string AtomItem::repr() {
            py::list res;
            res.append(left_bond->get_num());
            res.append(right_bond->get_num());
            return py::repr(res).cast<std::string>();
        }

        py::list AtomItem::pyrepr() {
            py::list res;
            res.append(left_bond->get_num());
            res.append(right_bond->get_num());
            return res;
        }

        std::string BondItem::repr() {
            py::list res;
            res.append(left_atom->get_num());
            res.append(right_atom->get_num());
            return py::repr(res).cast<std::string>();
        }

        py::list BondItem::pyrepr() {
            py::list res;
            res.append(left_atom->get_num());
            res.append(right_atom->get_num());
            return res;
        }

        void check_topology(AtomGroupSimple& atoms, BondGroupSimple& bonds) {
            int n_atoms = atoms.size();
            bool ok = true;
            for(int i = 0; i < n_atoms; ++i) {
                auto atomrepr = atoms[i].pyrepr();
                // fmt::print("Atom {} has bonds {}\n", i, atoms[i].repr());
                if ((atomrepr[0].cast<int>() != bonds[i - 1].get_num()) || (atomrepr[1].cast<int>() != bonds[i].get_num())) {
                    ok = false;
                    // fmt::print("Wrong!\n");
                }
                
                auto bondrepr = bonds[i].pyrepr();
                // fmt::print("Bond {} has atoms {}\n", i, bonds[i].repr());
                if ((bondrepr[0].cast<int>() != atoms[i].get_num()) || (bondrepr[1].cast<int>() != atoms[i + 1].get_num())) {
                    ok = false;
                    // fmt::print("Wrong!\n");
                }
            }

            if(!ok)
                throw std::runtime_error("Topology is corrupted");
        }

        std::string tuplelist_to_str(py::list x) {
            py::list ll;
            for(auto& item : x) {
                if(!py::isinstance<py::tuple>(item))
                    throw std::runtime_error("Not a list of tuples");
                ll.append(item.cast<py::list>());
            }
            return py::repr(ll).cast<std::string>();
        }

        std::string repr_matrix_buffer(const double* buf, const int& size1, const int& size2) {
            py::list res;
            for (unsigned int i = 0; i < size1; ++i) {
                py::list temp;
                for (unsigned int j = 0; j < size2; ++j)
                    temp.append(buf[i * size2 + j]);
                res.append(temp);
            }
            return py::repr(res).cast<std::string>();
        }
        
        std::string repr_matrix_buffer_addr(const double* buf, const std::vector<std::size_t>& addrs, const int& size2) {
            py::list res;
            for (const auto& addr : addrs) {
                py::list temp;
                for (unsigned int j = 0; j < size2; ++j)
                    temp.append(buf[addr + j]);
                res.append(temp);
            }
            return py::repr(res).cast<std::string>();
        }

        std::string repr_matrix(const Matrix& m) {
            py::list res;
            for (unsigned int i = 0; i < m.size1 (); ++ i) {
                py::list temp;
                for (unsigned int j = 0; j < m.size2 (); ++ j)
                    temp.append(py::cast(m(i, j)));
                res.append(temp);
            }
            return py::repr(res).cast<std::string>();
        }

        std::string repr_row(const Matrix& m, const int& row_idx) {
            py::list res;
            for (unsigned int j = 0; j < m.size2 (); ++ j)
                res.append(py::cast(m(row_idx, j)));
            return py::repr(res).cast<std::string>();
        }
    #endif
}
