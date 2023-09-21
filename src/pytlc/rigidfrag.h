#ifndef RIGIDFRAG_H
#define RIGIDFRAG_H

#include "utils.h"


using TlcUtils::AtomItem;
using TlcUtils::BondItem;
using TlcUtils::AtomGroup;
using TlcUtils::BondGroup;
using TlcUtils::AtomContainerStatic;
using TlcUtils::AtomContainerDynamic;
using TlcUtils::BondContainerStatic;
using TlcUtils::BondContainerDynamic;

enum class RFCase { atom, bond, general };
using Matrix = boost::numeric::ublas::matrix<double>;
using SEMatrix = boost::numeric::ublas::fixed_matrix<double, 4, 4>;
using XYZMatrix = boost::numeric::ublas::fixed_matrix<double, 4, 3>;
using FrameMatrix = boost::numeric::ublas::fixed_matrix<double, 3, 3, boost::numeric::ublas::column_major>;
using IMatrix = boost::numeric::ublas::identity_matrix<double>;

class RigidFrag {
    public:
        RigidFrag() noexcept {}
        void initialize(AtomItem left_atom, AtomItem right_atom, AtomGroup& catoms, BondGroup& cbonds);

        inline void write_tlc_data_init(AtomGroup& catoms, BondGroup& cbonds) noexcept {
            cbonds.tlc_lengths[3 * num] = cbonds.bonds[this->ebonds[0].get_num()];
            cbonds.tlc_lengths[3 * num + 2] = cbonds.bonds[this->ebonds[1].get_num()];
            catoms.tlc_vangles[3 * num] = catoms.vangles[this->eatoms[0].get_num()] * RAD2DEG;
        }

        void write_tlc_data_work(AtomGroup& catoms, BondGroup& cbonds) noexcept;
        void reconstruct(double* buffer, AtomGroup& catoms, BondGroup& cbonds) noexcept;

        void set_num(int new_num) noexcept { this->num = new_num; }

    private:
        inline void insert_iatom_xyz(Matrix& r, const Matrix& m, const int& idx) noexcept {
            r(idx, 0) = m(0, 3);
            r(idx, 1) = m(1, 3);
            r(idx, 2) = m(2, 3);
        }

        inline void build_frame_Lopt(FrameMatrix& m) noexcept {
            double *x = &m(0, 0), *y = &m(0, 1), *z = &m(0, 2);
            // x /= norm(x)
            cblas_daxpy(3, -cblas_ddot(3, x, 1, y, 1), x, 1, y, 1);
            cblas_dscal(3, 1 / cblas_dnrm2(3, y, 1), y, 1);
            TlcUtils::calc_cross(x, y, z);
        }
        
        inline void build_frame(double m[]) noexcept {
            double *x = &m[0], *y = &m[3], *z = &m[6];
            cblas_dscal(3, 1 / cblas_dnrm2(3, x, 1), x, 1);
            cblas_daxpy(3, -cblas_ddot(3, x, 1, y, 1), x, 1, y, 1);
            cblas_dscal(3, 1 / cblas_dnrm2(3, y, 1), y, 1);
            TlcUtils::calc_cross(x, y, z);
        }

        unsigned int num;
        RFCase mycase;
        AtomContainerDynamic iatoms;
        AtomContainerStatic<2> eatoms, satoms;
        BondContainerDynamic ibonds;
        BondContainerStatic<2> ebonds, sbonds;

        Matrix iatoms_xyz;
        std::size_t eatom_addr, satomL_addr, satomR_addr, firstiatom_addr, segA_length, segB_length;
        bool do_Bseg;
};


class RigidFragContainer {
    public:
        RigidFragContainer() noexcept {
            for(int i = 0; i < 3; ++i)
                items[i].set_num(i);
        }

        inline unsigned int size() const noexcept { return items.size(); }

        RigidFrag& operator [](int idx) noexcept { return items[idx]; }
        const RigidFrag& operator[](int idx) const noexcept { return items[idx]; }
        inline std::array<RigidFrag, 3>::iterator begin() noexcept { return items.begin(); }
        inline std::array<RigidFrag, 3>::iterator end() noexcept { return items.end(); }
        inline std::array<RigidFrag, 3>::const_iterator cbegin() const noexcept { return items.cbegin(); }
        inline std::array<RigidFrag, 3>::const_iterator cend() const noexcept { return items.cend(); }
            
    private:
        std::array<RigidFrag, 3> items;
};

#endif //RIGIDFRAG_H
