#ifndef RMSD_H
#define RMSD_H

#include <limits>

#include "utils.h"


template<class IC>
class RmsdCalculator {
    public:
        using isomorphism_container_t = IC;
        using single_isomorphism_t = typename IC::value_type;

        using Matrix = boost::numeric::ublas::matrix<double>;
        template <std::size_t N, std::size_t M> using FMatrix = boost::numeric::ublas::fixed_matrix<double, N, M>;
        template <std::size_t N> using FVector = boost::numeric::ublas::fixed_vector<double, N>;
        using ElemVector = boost::numeric::ublas::vector<int>;
        using XYZVector = boost::numeric::ublas::fixed_vector<int, 3>;
        using CoordVector = boost::numeric::ublas::fixed_vector<double, 3>;
        using MatrixSlice = boost::numeric::ublas::matrix_slice<Matrix>;

        using IndArray = boost::numeric::ublas::indirect_array<boost::numeric::ublas::vector<int>>;
        using Permutation = boost::numeric::ublas::permutation_matrix<int>;
        using MatrixView = boost::numeric::ublas::matrix_indirect<Matrix, IndArray>;
        template <class T> using VectorView = boost::numeric::ublas::vector_indirect<boost::numeric::ublas::vector<T>, IndArray>;
    
    private:
        int n_atoms;
        double cutoff;
        Matrix p_coord, q_coord;
        isomorphism_container_t* isomorphisms;
        single_isomorphism_t* simple_reorder;
        bool mirror_match;

    public:
        RmsdCalculator () : n_atoms(-1) { }

        RmsdCalculator(isomorphism_container_t& isomorphisms, single_isomorphism_t& simple_reorder, const double cutoff, const bool mirror_match)
            : n_atoms(simple_reorder.size()),
              cutoff(cutoff),
              p_coord(Matrix(n_atoms, 3)),
              q_coord(Matrix(n_atoms, 3)),
              isomorphisms(&isomorphisms),
              simple_reorder(&simple_reorder),
              mirror_match(mirror_match)
        { }

        inline double get_natoms() const noexcept { return n_atoms; }
        inline double get_cutoff() const noexcept { return cutoff; }

        double calc(const Matrix& p_all, const Matrix& q_all) {
            double best_rmsd = std::numeric_limits<double>::max();
            for(const auto& isomorphism : (*isomorphisms)) {
                for(int i = 0; i < n_atoms; ++i) { // Must apply all rearrangements right here!
                    for(int j = 0; j < 3; ++j) {
                        p_coord(i, j) = p_all((*simple_reorder).at(i), j);
                        q_coord(i, j) = q_all(isomorphism[i], j);
                    }
                }

                PyxyzUtils::make_centered(p_coord);
                PyxyzUtils::make_centered(q_coord);

                FMatrix<3, 3> C, V, W;
                FVector<3> S;
                noalias(C) = prod(trans(p_coord), q_coord);
                
                Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> C_ei(&C(0, 0));
                Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(C_ei);
                
                auto U_ei = svd.matrixU();
                for (unsigned int i = 0; i < 3; ++i)
                    for (unsigned int j = 0; j < 3; ++j)
                        V(i, j) = U_ei(i, j);
                auto V_ei = svd.matrixV();
                for (unsigned int i = 0; i < 3; ++i)
                    for (unsigned int j = 0; j < 3; ++j)
                        W(i, j) = V_ei(j, i);

                if (!mirror_match && (PyxyzUtils::det3x3(V) * PyxyzUtils::det3x3(W) < 0.0)) {
                    for(int i = 0; i < 3; i++)
                        V(i, 2) = -V(i, 2);
                }
                FMatrix<3, 3> U;
                noalias(U) = prod(V, W);
                p_coord = prod(p_coord, U);

                noalias(p_coord) -= q_coord;
                for (auto i = &p_coord(0, 0); i < &p_coord(n_atoms - 1, 2) + 1; ++i)
                    *i *= *i; // Nice
                
                double res = sqrt(std::accumulate(&p_coord(0, 0), &p_coord(n_atoms - 1, 2) + 1, 0.0) / n_atoms);
                if (res < best_rmsd)
                    best_rmsd = res;
                if (best_rmsd < cutoff)
                    break;
            }
            // std::cout << "RMSD = " << best_rmsd << std::endl;
            return best_rmsd;
        }

        py::tuple manual_calc(Matrix& p_all, Matrix& q_all) {
            double best_rmsd = std::numeric_limits<double>::max();
            FMatrix<3, 3> U_best;
            CoordVector q_center;
            for(const auto& isomorphism : (*isomorphisms)) {
                for(int i = 0; i < n_atoms; ++i) { // Must apply all rearrangements right here!
                    for(int j = 0; j < 3; ++j) {
                        p_coord(i, j) = p_all((*simple_reorder)[i], j);
                        q_coord(i, j) = q_all(isomorphism[i], j);
                    }
                }

                PyxyzUtils::make_centered(p_coord);
                PyxyzUtils::make_centered(q_coord, &q_center(0));

                FMatrix<3, 3> C, V, W;
                FVector<3> S;
                noalias(C) = prod(trans(p_coord), q_coord);
                
                Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> C_ei(&C(0, 0));
                Eigen::JacobiSVD<Eigen::Matrix3d, Eigen::ComputeThinU | Eigen::ComputeThinV> svd(C_ei);
                
                auto U_ei = svd.matrixU();
                for (unsigned int i = 0; i < 3; ++i)
                    for (unsigned int j = 0; j < 3; ++j)
                        V(i, j) = U_ei(i, j);
                auto V_ei = svd.matrixV();
                for (unsigned int i = 0; i < 3; ++i)
                    for (unsigned int j = 0; j < 3; ++j)
                        W(i, j) = V_ei(j, i);

                if (!mirror_match && (PyxyzUtils::det3x3(V) * PyxyzUtils::det3x3(W) < 0.0)) {
                    for(int i = 0; i < 3; i++)
                        V(i, 2) = -V(i, 2);
                }
                FMatrix<3, 3> U;
                noalias(U) = prod(V, W);
                p_coord = prod(p_coord, U);

                noalias(p_coord) -= q_coord;
                for (auto i = &p_coord(0, 0); i < &p_coord(n_atoms - 1, 2) + 1; ++i)
                    *i *= *i; // Nice
                
                double res = sqrt(std::accumulate(&p_coord(0, 0), &p_coord(n_atoms - 1, 2) + 1, 0.0) / n_atoms);
                if (res < best_rmsd) {
                    best_rmsd = res;
                    noalias(U_best) = U;
                }
                if (best_rmsd < cutoff)
                    break;
            }
            auto U_best_np = PyxyzUtils::ublas_matrix_to_np(U_best);
            auto q_center_np = PyxyzUtils::ublas_vector_to_np(q_center);
            return py::make_tuple(best_rmsd, U_best_np, q_center_np);
        }
};

#endif // RMSD_H
