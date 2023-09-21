#ifndef TLCMAINSOLVER_H
#define TLCMAINSOLVER_H

#include "utils.h"


class TLCMain {
    public:
        TLCMain();
        py::int_ solve(const py::list& bonds, const py::list& vangles, const py::list& tangles);
        int solve_internal(const std::array<double, 9>& bonds, const std::array<double, 9>& vangles, const std::array<double, 3>& tangles);
        void set_solution(py::array_t<double> solution, const int& sol_idx);
        #ifdef TLC_LOG
            py::str get_log() const
            { return py::cast(Utils::Logger::get_log()); }
        #endif

        inline void solution_to_buffer(double* buffer, const int& sol_idx) noexcept {
            const int ind[][2] = {{1,0},{2,0},{0,1},{1,1},{2,1},{0,2},{1,2},{2,2},{0,0}};
            for(int i = 0; i < 9; ++i) {
                if(ind[i][0] == 0)
                    cblas_dcopy(3, &r_soln_n[sol_idx][ind[i][1]][0], 1, &buffer[i * 3], 1);
                else if(ind[i][0] == 1)
                    cblas_dcopy(3, &r_soln_a[sol_idx][ind[i][1]][0], 1, &buffer[i * 3], 1);
                else // == 2
                    cblas_dcopy(3, &r_soln_c[sol_idx][ind[i][1]][0], 1, &buffer[i * 3], 1);
            }
        }

        inline void solution_to_buffer(double* buffer, const int& sol_idx, const std::array<unsigned int, 9>& map) noexcept {
            const int ind[][2] = {{1,0},{2,0},{0,1},{1,1},{2,1},{0,2},{1,2},{2,2},{0,0}};
            for(int i = 0; i < 9; ++i) {
                if(ind[i][0] == 0)
                    cblas_dcopy(3, &r_soln_n[sol_idx][ind[i][1]][0], 1, &buffer[map[i] * 3], 1);
                else if(ind[i][0] == 1)
                    cblas_dcopy(3, &r_soln_a[sol_idx][ind[i][1]][0], 1, &buffer[map[i] * 3], 1);
                else // == 2
                    cblas_dcopy(3, &r_soln_c[sol_idx][ind[i][1]][0], 1, &buffer[map[i] * 3], 1);
            }
        }
        bool some_removed;
        
    private:
        void initialize();
        int start_tlc();
        bool get_input_angles(const double r_n1[], const double r_a1[], const double r_a3[], const double r_c3[]);
        bool test_two_cone_existence_soln(const double, const double, const double, const double);
        void get_poly_coeff();
        void coord_from_poly_roots(int& n_soln, double roots[16], const double r_n1[3], const double r_a1[3], const double r_a3[3], const double r_c3[3]);

        void poly_mul_sub2(const double[5][5], const double[5][5], const double[5][5], const double[5][5], const int[2], const int[2], const int[2], const int[2], double[5][5], int[2]);
        void poly_mul2(const double [5][5], const double [5][5], const int [2], const int [2], double [5][5], int [2]);
        void poly_sub2(const double [5][5], const double [5][5], const int [2], const int [2], double [5][5], int [2]);
        void poly_mul_sub1(const double u1[17], const double u2[17], const double u3[17], const double u4[17], const int p1, const int p2, const int p3, const int p4, double u5[17], int& p5);
        void poly_mul1(const double u1[17], const double u2[17], const int p1, const int p2, double u3[17], int& p3);
        void poly_sub1(const double u1[17], const double u2[17], const int p1, const int p2, double u3[17], int& p3);

        double calc_t2(double t0);
        double calc_t1(double t0, double t2);
        void quaternion(const double axis[], const double& ang, double res[]);
        void rotation_matrix(const double q[], double U[][3]);
        double calc_bnd_ang(const double r1[], const double r2[]);
        double calc_dih_ang(const double r1[], const double r2[], const double r3[]);
        inline double get_vangle(double (&a)[3], double (&b)[3], double (&c)[3]) noexcept;

        
        double solutions[16][9][3];
        int nsolutions;
        double m_bonds[9];
        double m_vangles[9];
        double m_tangles[3];
        double len0[6], b_ang0[7], t_ang0[2];
        double C[3][3][3];

        double delta[3], xi[3], eta[3], alpha[3], theta[3]; // ATTENTION Some weird indexing in delta
        double aa13_min_sqr, aa13_max_sqr;
        
        double cos_alpha[3], sin_alpha[3], cos_theta[3], sin_theta[3];
        double cos_delta[3], sin_delta[3];
        double cos_xi[3], cos_eta[3], sin_xi[3], sin_eta[3];

        double r_a1a3[3], r_a1n1[3], r_a3c3[3];
        double b_a1a3[3], b_a1n1[3], b_a3c3[3];
        double len_na[3], len_ac[3], len_aa[3];
        double poly_coeff[17], Q[5][17], R[3][17];
        double r_soln_n[16][3][3], r_soln_a[16][3][3], r_soln_c[16][3][3];
};

#endif // TLCMAINSOLVER_H