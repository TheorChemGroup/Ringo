#include "tlcmain.h"

extern "C" void initialize_sturm_(double *tol_secant, int *max_iter_sturm, int *max_iter_secant);
extern "C" void solve_sturm_(int *p_order, int *n_root, double *poly_coeffs, double *roots);

#define MAX_SOLN 16


// This specifies the acceptablee errors in valence angles in TLC solution
#ifdef VALIDATION
    constexpr double THRESHOLD = Utils::STRONG_THRESHOLD / 100;
#else
    constexpr double THRESHOLD = Utils::WEAK_THRESHOLD / 5;
#endif

namespace TlcUtils {
    template <std::size_t N>
    inline void copy_array(const std::array<double, N> src, double (&dest)[N]) noexcept{
        std::copy(src.cbegin(), src.cend(), std::begin(dest));
    }

    void cross(const double v_A[3], const double v_B[3], double c_P[3]) {
        c_P[0] = v_A[1] * v_B[2] - v_A[2] * v_B[1];
        c_P[1] = -(v_A[0] * v_B[2] - v_A[2] * v_B[0]);
        c_P[2] = v_A[0] * v_B[1] - v_A[1] * v_B[0];
    }

    #ifdef EXT_LOG
        template <std::size_t N>
        std::string repr_array(const double (&a)[N])
        {
            py::list res;
            for (double i : a)
                res.append(i * RAD2DEG);
            return py::repr(res).cast<std::string>();
        }
    #endif

    #ifdef TLC_LOG
        template <typename T, std::size_t N>
        std::string print_array(const T(&a)[N])
        {
            std::string res;
            for (std::size_t i = 0; i < N; ++i)
            {
                res += boost::lexical_cast<std::string>(a[i]);
                res += " ";
            }
            return res;
        }

        template <typename T>
        std::string print_array(const T* a, const int size)
        {
            std::string res;
            for (std::size_t i = 0; i < size; ++i)
            {
                res += boost::lexical_cast<std::string>(a[i]);
                res += " ";
            }
            return res;
        }
    #endif
}


TLCMain::TLCMain()
{
    double tol_secant = 1.0e-15; // TODO: Optimize these parameters??
    int max_iter_sturm = 100, max_iter_secant = 20;
    initialize_sturm_(&tol_secant, &max_iter_sturm, &max_iter_secant);
}

py::int_ TLCMain::solve(const py::list& bonds, const py::list& vangles, const py::list& tangles) {
    for(int i = 0; i < 9; ++i) {
        m_bonds[i] = bonds[i].cast<double>();
        m_vangles[i] = vangles[i].cast<double>();
    }
    
    #ifdef TLC_LOG
        log("[CTOR] NEW CALL");
        std::string bonds_data = "";
        std::string vangle_data = "";
        for(int i = 0; i < 9; ++i){
            bonds_data += boost::lexical_cast<std::string>(m_bonds[i]);
            bonds_data += " ";
            vangle_data += boost::lexical_cast<std::string>(m_vangles[i]);
            vangle_data += " ";
        }
        log("INIT:: bonds = " + bonds_data);
        log("INIT:: vangles = " + vangle_data);
    #endif
    
    for(int i = 0; i < 3; ++i)
        m_tangles[i] = tangles[i].cast<double>();

    #ifdef TLC_LOG
        std::string tangle_data = "";
        for(int i = 0; i < 3; ++i) {
            tangle_data += boost::lexical_cast<std::string>(m_tangles[i]);
            tangle_data += " ";
        }
        log("INIT:: tangles = " + tangle_data);
    #endif
    return py::cast(start_tlc());
}

int TLCMain::solve_internal(const std::array<double, 9>& bonds, const std::array<double, 9>& vangles, const std::array<double, 3>& tangles) {
    TlcUtils::copy_array(bonds, m_bonds);
    TlcUtils::copy_array(vangles, m_vangles);
    TlcUtils::copy_array(tangles, m_tangles);
    // std::array<double, 3> xxx;
    // for(int i = 0; i < 3; ++i)
    //     xxx[i] = m_tangles[i];
    // throw std::runtime_error(fmt::format("source = {} dest = {}", py::repr(py::cast(tangles)).cast<std::string>(), py::repr(py::cast(xxx)).cast<std::string>()));
    return start_tlc();
}

int TLCMain::start_tlc() {
    double a1 = m_vangles[7]/180*M_PI;
	double a2 = m_vangles[8]/180*M_PI;
	double t = -m_tangles[2]/180*M_PI;
	double b6 = m_bonds[6];
	double b7 = m_bonds[7];
	double b8 = m_bonds[8];

	double at4[] = {0,0,0};
	double at3[] = {b6,0,0};
	double at2[] = {b7*cos(a1)-b8*(cos(a1)*cos(a2)-cos(t)*sin(a1)*sin(a2)),b7*sin(a1)+b8*(-cos(t)*cos(a1)*sin(a2)-cos(a2)*sin(a1)),-sin(t)*sin(a2)*b8};
	double at1[] = {b7*cos(a1),b7*sin(a1),0};
	double r_n[3][3] = {{at1[0],at1[1],at1[2]},{0,0,0},{0,0,0}};
	double r_a[3][3] = {{at2[0],at2[1],at2[2]},{0,0,0},{at3[0],at3[1],at3[2]}};
	double r_c[3][3] = {{0,0,0},{0,0,0},{at4[0],at4[1],at4[2]}};

	for(int i = 0; i < 6; ++i)
		len0[i] = m_bonds[i];
	
	for(int i = 0; i < 7; ++i)
		b_ang0[i] = m_vangles[i] / 180 * M_PI;
    
	for(int i = 0; i < 2; ++i)
		t_ang0[i] = -m_tangles[i] / 180 * M_PI;
    
    some_removed = false;

    #ifdef TLC_LOG
        log("START ROUTINE \"runtlc\"");
        log(fmt::format("IN:: b_len = {}", TlcUtils::print_array(len0)));
        log(fmt::format("IN:: b_ang = {}", TlcUtils::print_array(b_ang0)));
        log(fmt::format("IN:: t_ang = {}", TlcUtils::print_array(t_ang0)));
        log(fmt::format("IN:: r0_n(:, 1) = {}", TlcUtils::print_array(r_n[0])));
        log(fmt::format("IN:: r0_n(:, 2) = {}", TlcUtils::print_array(r_n[1])));
        log(fmt::format("IN:: r0_n(:, 3) = {}", TlcUtils::print_array(r_n[2])));
        log(fmt::format("IN:: r0_a(:, 1) = {}", TlcUtils::print_array(r_a[0])));
        log(fmt::format("IN:: r0_a(:, 2) = {}", TlcUtils::print_array(r_a[1])));
        log(fmt::format("IN:: r0_a(:, 3) = {}", TlcUtils::print_array(r_a[2])));
        log(fmt::format("IN:: r0_c(:, 1) = {}", TlcUtils::print_array(r_c[0])));
        log(fmt::format("IN:: r0_c(:, 2) = {}", TlcUtils::print_array(r_c[1])));
        log(fmt::format("IN:: r0_c(:, 3) = {}", TlcUtils::print_array(r_c[2])));
        log("END ROUTINE \"runtlc\"");
    #endif

    initialize();
    if(!get_input_angles(r_n[0], r_a[0], r_a[2], r_c[2]))
        return 0;
    
    get_poly_coeff();
    int deg_pol = MAX_SOLN, n_soln;
    double roots[MAX_SOLN];
    solve_sturm_(&deg_pol, &n_soln, poly_coeff, roots);

    #ifdef TLC_LOG
        log("START ROUTINE \"sturm_res\"");
        log(fmt::format("IN:: deg_pol = {}", deg_pol));
        log(fmt::format("IN:: poly_coeff = {}", TlcUtils::print_array(poly_coeff)));
        log(fmt::format("CHECK:: roots = {}", TlcUtils::print_array(roots)));
        log(fmt::format("CHECK:: n_soln = {}", n_soln));
        log("END ROUTINE \"sturm_res\"");
    #endif

    if(n_soln > 0)
        coord_from_poly_roots(n_soln, roots, r_n[0], r_a[0], r_a[2], r_c[2]); // Might decrease the n_soln value
    return n_soln;
}

void TLCMain::initialize() {
    double rr_c1[] = {0.0, 0.0, 0.0};
    double axis[] = {1.0, 0.0, 0.0};
    for(int i = 0; i < 2; ++i) {
        double rr_a1[] = {cos(b_ang0[3*i+1])*len0[3*i], sin(b_ang0[3*i+1])*len0[3*i], 0.0};
        double rr_n2[] = {len0[3*i+1], 0.0, 0.0};
        double rr_c1a1[3];
        for(int j = 0; j < 3; ++j)
            rr_c1a1[j] = rr_a1[j] - rr_c1[j];
        double rr_n2a2_ref[] = {-cos(b_ang0[3*i+2])*len0[3*i+2], sin(b_ang0[3*i+2])*len0[3*i+2], 0.0};

        double p[4];
        quaternion(axis, t_ang0[i] * 0.25, p);
        
        double Us[3][3];
        rotation_matrix(p, Us); // std::cout << fmt::format("Double check p = {}", TlcUtils::print_array(p)) << std::endl;

        double rr_a2[3];
        cblas_dcopy(3, &rr_n2[0], 1, &rr_a2[0], 1);
        // rr_a2(:) = matmul(Us, rr_n2a2_ref) + rr_n2(:)
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &rr_n2a2_ref[0], 1, 1.0, &rr_a2[0], 1);
        #ifdef TLC_LOG
            log("START ROUTINE \"matmul\"");
            log(fmt::format("IN:: rr_n2 = {}", TlcUtils::print_array(rr_n2)));
            log(fmt::format("IN:: rr_n2a2_ref = {}", TlcUtils::print_array(rr_n2a2_ref)));
            log(fmt::format("CHECK:: rr_a2 = {}", TlcUtils::print_array(rr_a2)));
            log("END ROUTINE \"matmul\"");
        #endif

        double dr[3];
        for(int j = 0; j < 3; ++j)
            dr[j] = rr_a2[j] - rr_a1[j];
        len_aa[i + 1] = cblas_dnrm2(3, &dr[0], 1);

        double bb_c1a1[3], bb_a1a2[3], bb_a2n2[3];
        for(int j = 0; j < 3; ++j) {
            bb_c1a1[j] = rr_c1a1[j] / len0[3 * i];
            bb_a1a2[j] = dr[j] / len_aa[i + 1];
            bb_a2n2[j] = (rr_n2[j] - rr_a2[j])/len0[3*i+2];
        }

        // for(int j = 0; j < 3; ++j)
        //     bb_a1a2[j] = -bb_a1a2[j];
        xi[i + 1] = M_PI - calc_bnd_ang(bb_a1a2, bb_a2n2); // Better?
        // for(int j = 0; j < 3; ++j) {
        //     bb_a1a2[j] = -bb_a1a2[j];
        //     bb_c1a1[j] = -bb_c1a1[j];
        // }
        eta[i] = M_PI - calc_bnd_ang(bb_a1a2, bb_c1a1);
        // for(int j = 0; j < 3; ++j)
        //     bb_c1a1[j] = -bb_c1a1[j];
        delta[i] = M_PI - calc_dih_ang(bb_c1a1, bb_a1a2, bb_a2n2);
    }

    double a_min = b_ang0[3] - (xi[1] + eta[1]);
    double a_max = std::min(b_ang0[3] + (xi[1] + eta[1]), M_PI);
    aa13_min_sqr = pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0*len_aa[1]*len_aa[2]*cos(a_min);
    aa13_max_sqr = pow(len_aa[1], 2) + pow(len_aa[2], 2) - 2.0*len_aa[1]*len_aa[2]*cos(a_max);

    #ifdef TLC_LOG
        log("START ROUTINE \"initialize_loop_closure\"");
        log(fmt::format("IN:: b_len = {}", TlcUtils::print_array(len0, 6)));
        log(fmt::format("IN:: b_ang = {}", TlcUtils::print_array(b_ang0, 7)));
        log(fmt::format("IN:: t_ang = {}", TlcUtils::print_array(t_ang0, 2)));
        log(fmt::format("CHECK:: len_aa(2) = {}", len_aa[1]));
        log(fmt::format("CHECK:: len_aa(3) = {}", len_aa[2]));
        log(fmt::format("CHECK:: a_min = {}", a_min));
        log(fmt::format("CHECK:: a_max = {}", a_max));
        log(fmt::format("CHECK:: aa13_min_sqr = {}", aa13_min_sqr));
        log(fmt::format("CHECK:: aa13_max_sqr = {}", aa13_max_sqr));
        log("END ROUTINE \"initialize_loop_closure\"");
    #endif
}

bool TLCMain::get_input_angles(const double r_n1[], const double r_a1[], const double r_a3[], const double r_c3[]) {
    for(int i = 0; i < 3; ++i)
        r_a1a3[i] = r_a3[i] - r_a1[i];
    
    len_aa[0] = cblas_dnrm2(3, &r_a1a3[0], 1);
    double dr_sqr = pow(len_aa[0], 2); // Get rid of it
    
    if ((dr_sqr < aa13_min_sqr) || (dr_sqr > aa13_max_sqr))
        return 0; // No solutions

    for(int i = 0; i < 3; ++i)
        r_a1n1[i] = r_n1[i] - r_a1[i];
    len_na[0] = cblas_dnrm2(3, &r_a1n1[0], 1);
    len_na[1] = len0[2];
    len_na[2] = len0[5];

    for(int i = 0; i < 3; ++i)
        r_a3c3[i] = r_c3[i] - r_a3[i];
    len_ac[0] = len0[0];
    len_ac[1] = len0[3];
    len_ac[2] = cblas_dnrm2(3, &r_a3c3[0], 1);

    for(int i = 0; i < 3; ++i) {
        b_a1n1[i] = r_a1n1[i]/len_na[0];
        b_a3c3[i] = r_a3c3[i]/len_ac[2];
        b_a1a3[i] = r_a1a3[i]/len_aa[0];
    }

    #ifdef TLC_LOG
    log("START ROUTINE \"get_input_angles_temp\"");
    log(fmt::format("IN:: r_n1 = {}", TlcUtils::print_array(r_n1, 3)));
    log(fmt::format("IN:: r_a1 = {}", TlcUtils::print_array(r_a1, 3)));
    log(fmt::format("IN:: r_a3 = {}", TlcUtils::print_array(r_a3, 3)));
    log(fmt::format("IN:: r_c3 = {}", TlcUtils::print_array(r_c3, 3)));
    log(fmt::format("CHECK:: r_a1n1 = {}", TlcUtils::print_array(r_a1n1)));
    log(fmt::format("CHECK:: r_a3c3 = {}", TlcUtils::print_array(r_a3c3)));
    log(fmt::format("CHECK:: r_a1a3 = {}", TlcUtils::print_array(r_a1a3)));
    log(fmt::format("CHECK:: b_a1n1 = {}", TlcUtils::print_array(b_a1n1)));
    log(fmt::format("CHECK:: b_a3c3 = {}", TlcUtils::print_array(b_a3c3)));
    log(fmt::format("CHECK:: b_a1a3 = {}", TlcUtils::print_array(b_a1a3)));
    log(fmt::format("CHECK:: len_na = {}", TlcUtils::print_array(len_na)));
    log(fmt::format("CHECK:: len_ac = {}", TlcUtils::print_array(len_ac)));
    log(fmt::format("CHECK:: len_aa = {}", TlcUtils::print_array(len_aa)));
    log("END ROUTINE \"get_input_angles_temp\"");
    #endif

    for(int i = 0; i < 3; ++i)
        b_a1n1[i] = -b_a1n1[i];
    delta[2] = calc_dih_ang(b_a1n1, b_a1a3, b_a3c3);
    for(int i = 0; i < 3; ++i) {
        b_a1n1[i] = -b_a1n1[i];
        b_a1a3[i] = -b_a1a3[i];
    }
    xi[0] = calc_bnd_ang(b_a1a3, b_a1n1);
    for(int i = 0; i < 3; ++i)
        b_a1a3[i] = -b_a1a3[i];
    eta[2] = calc_bnd_ang(b_a1a3, b_a3c3);

    theta[0] = b_ang0[0];
    theta[1] = b_ang0[3];
    theta[2] = b_ang0[6];
    for(int i = 0; i < 3; ++i) {
        cos_delta[i] = cos(delta[i]);
        sin_delta[i] = sin(delta[i]);
        cos_xi[i] = cos(xi[i]);
        sin_xi[i] = sin(xi[i]);
        cos_eta[i] = cos(eta[i]);
        sin_eta[i] = sin(eta[i]);
        cos_theta[i] = cos(theta[i]);
    }

    cos_alpha[0] = -(pow(len_aa[0], 2) + pow(len_aa[1], 2) - pow(len_aa[2], 2))/(2.0*len_aa[0]*len_aa[1]);
    alpha[0] = acos(cos_alpha[0]);
    sin_alpha[0] = sin(alpha[0]);
    cos_alpha[1] = (pow(len_aa[1], 2) + pow(len_aa[2], 2) - pow(len_aa[0], 2))/(2.0*len_aa[1]*len_aa[2]);
    alpha[1] = acos(cos_alpha[1]);
    sin_alpha[1] = sin(alpha[1]);
    alpha[2] = M_PI - alpha[0] + alpha[1];
    cos_alpha[2] = cos(alpha[2]);
    sin_alpha[2] = sin(alpha[2]);

    #ifdef TLC_LOG
    log("START ROUTINE \"get_input_angles\"");
    log(fmt::format("IN:: r_n1 = {}", TlcUtils::print_array(r_n1, 3)));
    log(fmt::format("IN:: r_a1 = {}", TlcUtils::print_array(r_a1, 3)));
    log(fmt::format("IN:: r_a3 = {}", TlcUtils::print_array(r_a3, 3)));
    log(fmt::format("IN:: r_c3 = {}", TlcUtils::print_array(r_c3, 3)));
    log(fmt::format("CHECK:: xi = {}", TlcUtils::print_array(xi)));
    log(fmt::format("CHECK:: eta = {}", TlcUtils::print_array(eta)));
    log(fmt::format("CHECK:: delta = {}", TlcUtils::print_array(delta)));
    log(fmt::format("CHECK:: theta = {}", TlcUtils::print_array(theta)));
    log(fmt::format("CHECK:: alpha = {}", TlcUtils::print_array(alpha)));
    log("END ROUTINE \"get_input_angles\"");
    #endif

    for(int i = 0; i < 3; ++i) {
        if(!test_two_cone_existence_soln(theta[i], xi[i], eta[i], alpha[i]))
            return false;
    }
    return true;
}

void TLCMain::get_poly_coeff() {
    double A[5], B[9][3];
    double A21, A22, A31, A32, A41, A42;
    for(int i = 0; i < 3; ++i) {
        A[0] = cos_alpha[i]*cos_xi[i]*cos_eta[i] - cos_theta[i];
        A[1] = -sin_alpha[i]*cos_xi[i]*sin_eta[i];
        A[2] = sin_alpha[i]*sin_xi[i]*cos_eta[i];
        A[3] = sin_xi[i]*sin_eta[i];
        A[4] = A[3]*cos_alpha[i];

        #ifdef TLC_LOG
        log("START ROUTINE \"get_poly_coeff_pos0\"");
        log(fmt::format("IN:: i = {}", i));
        log(fmt::format("CHECK:: A0 = {}", A[0]));
        log(fmt::format("CHECK:: A1 = {}", A[1]));
        log(fmt::format("CHECK:: A2 = {}", A[2]));
        log(fmt::format("CHECK:: A3 = {}", A[3]));
        log(fmt::format("CHECK:: A4 = {}", A[4]));
        log("END ROUTINE \"get_poly_coeff_pos0\"");
        #endif
        
        int j = i - 1;
        if(i == 0)
            j = 2;

        A21 = A[2]*cos_delta[j];
        A22 = A[2]*sin_delta[j];
        A31 = A[3]*cos_delta[j];
        A32 = A[3]*sin_delta[j];
        A42 = A[4]*sin_delta[j];
        A41 = A[4]*cos_delta[j];

        #ifdef TLC_LOG
        log("START ROUTINE \"get_poly_coeff_pos00\"");
        log(fmt::format("IN:: i = {}", i));
        log(fmt::format("CHECK:: A21 = {}", A21));
        log(fmt::format("CHECK:: A22 = {}", A22));
        log(fmt::format("CHECK:: A31 = {}", A31));
        log(fmt::format("CHECK:: A32 = {}", A32));
        log(fmt::format("CHECK:: A42 = {}", A42));
        log(fmt::format("CHECK:: A41 = {}", A41));
        log("END ROUTINE \"get_poly_coeff_pos00\"");
        #endif

        B[0][i] = A[0] + A22 + A31;
        B[1][i] = 2.0*(A[1] + A42);
        B[2][i] = 2.0*(A32 - A21);
        B[3][i] = -4.0*A41;
        B[4][i] = A[0] + A22 - A31;
        B[5][i] = A[0] - A22 - A31;
        B[6][i] = -2.0*(A21 + A32);
        B[7][i] = 2.0*(A[1] - A42);
        B[8][i] = A[0] - A22 + A31;
    }

    #ifdef TLC_LOG
    log("START ROUTINE \"get_poly_coeff_pos1\"");
    log(fmt::format("CHECK:: B0 = {}", TlcUtils::print_array(B[0])));
    log(fmt::format("CHECK:: B1 = {}", TlcUtils::print_array(B[1])));
    log(fmt::format("CHECK:: B2 = {}", TlcUtils::print_array(B[2])));
    log(fmt::format("CHECK:: B3 = {}", TlcUtils::print_array(B[3])));
    log(fmt::format("CHECK:: B4 = {}", TlcUtils::print_array(B[4])));
    log(fmt::format("CHECK:: B5 = {}", TlcUtils::print_array(B[5])));
    log(fmt::format("CHECK:: B6 = {}", TlcUtils::print_array(B[6])));
    log(fmt::format("CHECK:: B7 = {}", TlcUtils::print_array(B[7])));
    log(fmt::format("CHECK:: B8 = {}", TlcUtils::print_array(B[8])));
    log("END ROUTINE \"get_poly_coeff_pos1\"");
    #endif

    memset(C, 0, sizeof(C));
    C[0][0][0] = B[0][0]; C[0][0][1] = B[2][0]; C[0][0][2] = B[5][0];
    C[1][0][0] = B[1][0]; C[1][0][1] = B[3][0]; C[1][0][2] = B[7][0];
    C[2][0][0] = B[4][0]; C[2][0][1] = B[6][0]; C[2][0][2] = B[8][0];
    for(int i = 1; i < 3; ++i) {
        C[0][i][0] = B[0][i]; C[0][i][1] = B[1][i]; C[0][i][2] = B[4][i];
        C[1][i][0] = B[2][i]; C[1][i][1] = B[3][i]; C[1][i][2] = B[6][i];
        C[2][i][0] = B[5][i]; C[2][i][1] = B[7][i]; C[2][i][2] = B[8][i];
    }

    double u11[5][5], u12[5][5], u13[5][5], u31[5][5], u32[5][5], u33[5][5];
    for(int i = 0; i < 5; ++i) {
        for(int j = 0; j < 5; ++j) {
            u11[i][j] = 0.0;
            u12[i][j] = 0.0;
            u13[i][j] = 0.0;
            u31[i][j] = 0.0;
            u32[i][j] = 0.0;
            u33[i][j] = 0.0;
        }
    }

    for(int i = 0; i < 3; ++i) {
        u11[0][i] = C[0][0][i];
        u12[0][i] = C[1][0][i];
        u13[0][i] = C[2][0][i];
        u31[i][0] = C[0][1][i];
        u32[i][0] = C[1][1][i];
        u33[i][0] = C[2][1][i];
    }

    #ifdef TLC_LOG
    log("START ROUTINE \"get_poly_coeff_pos2\"");
    log(fmt::format("CHECK:: u11(:, 0) = {}", TlcUtils::print_array(u11[0])));
    log(fmt::format("CHECK:: u12(:, 0) = {}", TlcUtils::print_array(u12[0])));
    log(fmt::format("CHECK:: u13(:, 0) = {}", TlcUtils::print_array(u13[0])));
    log(fmt::format("CHECK:: u31(:, 0) = {}", TlcUtils::print_array(u31[0])));
    log(fmt::format("CHECK:: u32(:, 0) = {}", TlcUtils::print_array(u32[0])));
    log(fmt::format("CHECK:: u33(:, 0) = {}", TlcUtils::print_array(u33[0])));

    log(fmt::format("CHECK:: u11(:, 1) = {}", TlcUtils::print_array(u11[1])));
    log(fmt::format("CHECK:: u12(:, 1) = {}", TlcUtils::print_array(u12[1])));
    log(fmt::format("CHECK:: u13(:, 1) = {}", TlcUtils::print_array(u13[1])));
    log(fmt::format("CHECK:: u31(:, 1) = {}", TlcUtils::print_array(u31[1])));
    log(fmt::format("CHECK:: u32(:, 1) = {}", TlcUtils::print_array(u32[1])));
    log(fmt::format("CHECK:: u33(:, 1) = {}", TlcUtils::print_array(u33[1])));

    log(fmt::format("CHECK:: u11(:, 2) = {}", TlcUtils::print_array(u11[2])));
    log(fmt::format("CHECK:: u12(:, 2) = {}", TlcUtils::print_array(u12[2])));
    log(fmt::format("CHECK:: u13(:, 2) = {}", TlcUtils::print_array(u13[2])));
    log(fmt::format("CHECK:: u31(:, 2) = {}", TlcUtils::print_array(u31[2])));
    log(fmt::format("CHECK:: u32(:, 2) = {}", TlcUtils::print_array(u32[2])));
    log(fmt::format("CHECK:: u33(:, 2) = {}", TlcUtils::print_array(u33[2])));

    log(fmt::format("CHECK:: u11(:, 3) = {}", TlcUtils::print_array(u11[3])));
    log(fmt::format("CHECK:: u12(:, 3) = {}", TlcUtils::print_array(u12[3])));
    log(fmt::format("CHECK:: u13(:, 3) = {}", TlcUtils::print_array(u13[3])));
    log(fmt::format("CHECK:: u31(:, 3) = {}", TlcUtils::print_array(u31[3])));
    log(fmt::format("CHECK:: u32(:, 3) = {}", TlcUtils::print_array(u32[3])));
    log(fmt::format("CHECK:: u33(:, 3) = {}", TlcUtils::print_array(u33[3])));

    log(fmt::format("CHECK:: u11(:, 4) = {}", TlcUtils::print_array(u11[4])));
    log(fmt::format("CHECK:: u12(:, 4) = {}", TlcUtils::print_array(u12[4])));
    log(fmt::format("CHECK:: u13(:, 4) = {}", TlcUtils::print_array(u13[4])));
    log(fmt::format("CHECK:: u31(:, 4) = {}", TlcUtils::print_array(u31[4])));
    log(fmt::format("CHECK:: u32(:, 4) = {}", TlcUtils::print_array(u32[4])));
    log(fmt::format("CHECK:: u33(:, 4) = {}", TlcUtils::print_array(u33[4])));
    log("END ROUTINE \"get_poly_coeff_pos2\"");
    #endif

    int p1[2] = { 2, 0 };
    int p3[2] = { 0, 2 };
    int p_um1[2], p_um2[2], p_um3[2], p_um4[2], p_um5[2], p_um6[2], p_Q[2];
    double um1[5][5], um2[5][5], um3[5][5], um4[5][5], um5[5][5], um6[5][5], q_tmp[5][5];
    
    poly_mul_sub2(u32, u32, u31, u33, p3, p3, p3, p3, um1, p_um1);
    poly_mul_sub2(u12, u32, u11, u33, p1, p3, p1, p3, um2, p_um2);
    poly_mul_sub2(u12, u33, u13, u32, p1, p3, p1, p3, um3, p_um3);
    poly_mul_sub2(u11, u33, u31, u13, p1, p3, p3, p1, um4, p_um4);
    poly_mul_sub2(u13, um1, u33, um2, p1, p_um1, p3, p_um2, um5, p_um5);
    poly_mul_sub2(u13, um4, u12, um3, p1, p_um4, p1, p_um3, um6, p_um6);
    poly_mul_sub2(u11, um5, u31, um6, p1, p_um5, p3, p_um6, q_tmp, p_Q);
    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 17; ++j)
            Q[i][j] = 0.0;

    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
            Q[i][j] = q_tmp[i][j];
    
    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
            Q[i][j] = q_tmp[i][j];
    
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 17; ++j)
            R[i][j] = 0.0;
    
    for(int i = 0; i < 3; ++i) {
        R[0][i] = C[0][2][i];
        R[1][i] = C[1][2][i];
        R[2][i] = C[2][2][i];
    }

    double f1[17], f2[17], f3[17], f4[17], f5[17], f6[17], f7[17], f8[17], f9[17], f10[17], f11[17], f12[17], f13[17], f14[17], f15[17], f16[17], f17[17], f18[17], f19[17], f20[17], f21[17], f22[17], f23[17], f24[17], f25[17], f26[17];
    int p2 = 2, p4 = 4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7, p_f8, p_f9, p_f10, p_f11, p_f12, p_f13, p_f14, p_f15, p_f16, p_f17, p_f18, p_f19, p_f20, p_f21, p_f22, p_f23, p_f24, p_f25, p_f26, p_final;

    #ifdef TLC_LOG
    log("START ROUTINE \"get_poly_coeff_pos3\"");
    log(fmt::format("CHECK:: R(:, 0) = {}", TlcUtils::print_array(R[0])));
    log(fmt::format("CHECK:: R(:, 1) = {}", TlcUtils::print_array(R[1])));
    log(fmt::format("CHECK:: R(:, 2) = {}", TlcUtils::print_array(R[2])));
    log("END ROUTINE \"get_poly_coeff_pos3\"");
    #endif

    poly_mul_sub1(R[1], R[1], R[0], R[2], p2, p2, p2, p2, f1, p_f1);
    poly_mul1(R[1], R[2], p2, p2, f2, p_f2);
    poly_mul_sub1(R[1], f1, R[0], f2, p2, p_f1, p2, p_f2, f3, p_f3);
    poly_mul1(R[2], f1, p2, p_f1, f4, p_f4);
    poly_mul_sub1(R[1], f3, R[0], f4, p2, p_f3, p2, p_f4, f5, p_f5);

    poly_mul_sub1(Q[1], R[1], Q[0], R[2], p4, p2, p4, p2, f6, p_f6);
    poly_mul_sub1(Q[2], f1, R[2], f6, p4, p_f1, p2, p_f6, f7, p_f7);
    poly_mul_sub1(Q[3], f3, R[2], f7, p4, p_f3, p2, p_f7, f8, p_f8);
    poly_mul_sub1(Q[4], f5, R[2], f8, p4, p_f5, p2, p_f8, f9, p_f9);

    poly_mul_sub1(Q[3], R[1], Q[4], R[0], p4, p2, p4, p2, f10, p_f10);
    poly_mul_sub1(Q[2], f1, R[0], f10, p4, p_f1, p2, p_f10, f11, p_f11);
    poly_mul_sub1(Q[1], f3, R[0], f11, p4, p_f3, p2, p_f11, f12, p_f12);
    
    poly_mul_sub1(Q[2], R[1], Q[1], R[2], p4, p2, p4, p2, f13, p_f13);
    poly_mul_sub1(Q[3], f1, R[2], f13, p4, p_f1, p2, p_f13, f14, p_f14);
    poly_mul_sub1(Q[3], R[1], Q[2], R[2], p4, p2, p4, p2, f15, p_f15);
    poly_mul_sub1(Q[4], f1, R[2], f15, p4, p_f1, p2, p_f15, f16, p_f16);
    poly_mul_sub1(Q[1], f14, Q[0], f16, p4, p_f14, p4, p_f16, f17, p_f17);
    
    poly_mul_sub1(Q[2], R[2], Q[3], R[1], p4, p2, p4, p2, f18, p_f18);
    poly_mul_sub1(Q[1], R[2], Q[3], R[0], p4, p2, p4, p2, f19, p_f19);
    poly_mul_sub1(Q[3], f19, Q[2], f18, p4, p_f19, p4, p_f18, f20, p_f20);
    poly_mul_sub1(Q[1], R[1], Q[2], R[0], p4, p2, p4, p2, f21, p_f21);
    poly_mul1(Q[4], f21, p4, p_f21, f22, p_f22);
    poly_sub1(f20, f22, p_f20, p_f22, f23, p_f23);
    poly_mul1(R[0], f23, p2, p_f23, f24, p_f24);
    poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25);
    poly_mul_sub1(Q[4], f12, R[2], f25, p4, p_f12, p2, p_f25, f26, p_f26);
    poly_mul_sub1(Q[0], f9, R[0], f26, p4, p_f9, p2, p_f26, poly_coeff, p_final);

    if (poly_coeff[16] < 0.0)
        for(int i = 0; i < 17; ++i)
            poly_coeff[i] = -poly_coeff[i];
    
    #ifdef TLC_LOG
    log("START ROUTINE \"get_poly_coeff_final\"");
    log(fmt::format("CHECK:: poly_coeff = {}", TlcUtils::print_array(poly_coeff)));
    log("END ROUTINE \"get_poly_coeff_final\"");
    #endif
}

void TLCMain::poly_mul_sub2(const double u1[5][5], const double u2[5][5], const double u3[5][5], const double u4[5][5], const int p1[2], const int p2[2], const int p3[2], const int p4[2], double u5[5][5], int p5[2]) {
    double d1[5][5], d2[5][5];
    int pd1[2], pd2[2];
    poly_mul2(u1, u2, p1, p2, d1, pd1);
    poly_mul2(u3, u4, p3, p4, d2, pd2);
    poly_sub2(d1, d2, pd1, pd2, u5, p5);
}

void TLCMain::poly_mul2(const double u1[5][5], const double u2[5][5], const int p1[2], const int p2[2], double u3[5][5], int p3[2]) {
    for(int i = 0; i < 2; ++i)
        p3[i] = p1[i] + p2[i];
    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
            u3[i][j] = 0.0;
    int p11, p12, p21, p22;
    p11 = p1[0];
    p12 = p1[1];
    p21 = p2[0];
    p22 = p2[1];
    double u1ij;
    int i3, j3;
    for(int i1 = 0; i1 < p12+1; ++i1) {
        for(int j1 = 0; j1 < p11+1; ++j1) {
            u1ij = u1[i1][j1];
            for(int i2 = 0; i2 < p22+1; ++i2) {
                i3 = i1 + i2;
                for(int j2 = 0; j2 < p21+1; ++j2) {
                    j3 = j1 + j2;
                    u3[i3][j3] += u1ij*u2[i2][j2];
                }
            }
        }
    }
    
    #ifdef TLC_LOG
    log("START ROUTINE \"poly_mul2\"");
    log(fmt::format("IN:: u1(:, 0) = {}", TlcUtils::print_array(u1[0])));
    log(fmt::format("IN:: u1(:, 1) = {}", TlcUtils::print_array(u1[1])));
    log(fmt::format("IN:: u1(:, 2) = {}", TlcUtils::print_array(u1[2])));
    log(fmt::format("IN:: u1(:, 3) = {}", TlcUtils::print_array(u1[3])));
    log(fmt::format("IN:: u1(:, 4) = {}", TlcUtils::print_array(u1[4])));
    log(fmt::format("IN:: u2(:, 0) = {}", TlcUtils::print_array(u2[0])));
    log(fmt::format("IN:: u2(:, 1) = {}", TlcUtils::print_array(u2[1])));
    log(fmt::format("IN:: u2(:, 2) = {}", TlcUtils::print_array(u2[2])));
    log(fmt::format("IN:: u2(:, 3) = {}", TlcUtils::print_array(u2[3])));
    log(fmt::format("IN:: u2(:, 4) = {}", TlcUtils::print_array(u2[4])));
    log(fmt::format("IN:: p1 = {}", TlcUtils::print_array(p1, 2)));
    log(fmt::format("IN:: p2 = {}", TlcUtils::print_array(p2, 2)));
    log(fmt::format("CHECK:: u3(:, 0) = {}", TlcUtils::print_array(u3[0])));
    log(fmt::format("CHECK:: u3(:, 1) = {}", TlcUtils::print_array(u3[1])));
    log(fmt::format("CHECK:: u3(:, 2) = {}", TlcUtils::print_array(u3[2])));
    log(fmt::format("CHECK:: u3(:, 3) = {}", TlcUtils::print_array(u3[3])));
    log(fmt::format("CHECK:: u3(:, 4) = {}", TlcUtils::print_array(u3[4])));
    log(fmt::format("CHECK:: p3 = {}", TlcUtils::print_array(p3, 2)));
    log("END ROUTINE \"poly_mul2\"");
    #endif
}

void TLCMain::poly_sub2(const double u1[5][5], const double u2[5][5], const int p1[2], const int p2[2], double u3[5][5], int p3[2]) {
    for(int i = 0; i < 5; ++i)
        for(int j = 0; j < 5; ++j)
            u3[i][j] = 0.0;
    int p11, p12, p21, p22;
    p11 = p1[0];
    p12 = p1[1];
    p21 = p2[0];
    p22 = p2[1];
    p3[0] = std::max(p11, p21);
    p3[1] = std::max(p12, p22);
    #ifdef TLC_LOG
    log("START ROUTINE \"poly_sub2\"");
    #endif
    bool i1_ok, i2_ok;
    for(int i = 0; i < p3[1] + 1; ++i) {
        i1_ok = (i > p12);
        i2_ok = (i > p22);
        for(int j = 0; j < p3[0] + 1; ++j) {
            if (i2_ok || (j > p21)) {
                #ifdef TLC_LOG
                if((j == 4) && (i == 2))
                    log(fmt::format("CHECK:: case = 0"));
                #endif
                u3[i][j] = u1[i][j];
            } else if (i1_ok || (j > p11)) {
                #ifdef TLC_LOG
                if((j == 4) && (i == 2))
                    log(fmt::format("CHECK:: case = 1"));
                #endif
                u3[i][j] = -u2[i][j];
            } else {
                #ifdef TLC_LOG
                if((j == 4) && (i == 2))
                    log(fmt::format("CHECK:: case = 2"));
                #endif
                u3[i][j] = u1[i][j] - u2[i][j];
            }
        }
    }

    #ifdef TLC_LOG
    log(fmt::format("IN:: u1(:, 0) = {}", TlcUtils::print_array(u1[0])));
    log(fmt::format("IN:: u1(:, 1) = {}", TlcUtils::print_array(u1[1])));
    log(fmt::format("IN:: u1(:, 2) = {}", TlcUtils::print_array(u1[2])));
    log(fmt::format("IN:: u1(:, 3) = {}", TlcUtils::print_array(u1[3])));
    log(fmt::format("IN:: u1(:, 4) = {}", TlcUtils::print_array(u1[4])));
    log(fmt::format("IN:: u2(:, 0) = {}", TlcUtils::print_array(u2[0])));
    log(fmt::format("IN:: u2(:, 1) = {}", TlcUtils::print_array(u2[1])));
    log(fmt::format("IN:: u2(:, 2) = {}", TlcUtils::print_array(u2[2])));
    log(fmt::format("IN:: u2(:, 3) = {}", TlcUtils::print_array(u2[3])));
    log(fmt::format("IN:: u2(:, 4) = {}", TlcUtils::print_array(u2[4])));
    log(fmt::format("IN:: p1 = {}", TlcUtils::print_array(p1, 2)));
    log(fmt::format("IN:: p2 = {}", TlcUtils::print_array(p2, 2)));
    log(fmt::format("CHECK:: u3(:, 0) = {}", TlcUtils::print_array(u3[0])));
    log(fmt::format("CHECK:: u3(:, 1) = {}", TlcUtils::print_array(u3[1])));
    log(fmt::format("CHECK:: u3(:, 2) = {}", TlcUtils::print_array(u3[2])));
    log(fmt::format("CHECK:: u3(:, 3) = {}", TlcUtils::print_array(u3[3])));
    log(fmt::format("CHECK:: u3(:, 4) = {}", TlcUtils::print_array(u3[4])));
    log(fmt::format("CHECK:: p3 = {}", TlcUtils::print_array(p3, 2)));
    log("END ROUTINE \"poly_sub2\"");
    #endif
}

void TLCMain::poly_mul_sub1(const double u1[17], const double u2[17], const double u3[17], const double u4[17], const int p1, const int p2, const int p3, const int p4, double u5[17], int& p5) {
    double d1[17], d2[17];
    int pd1, pd2;
    poly_mul1(u1, u2, p1, p2, d1, pd1);
    poly_mul1(u3, u4, p3, p4, d2, pd2);
    poly_sub1(d1, d2, pd1, pd2, u5, p5);
}

void TLCMain::poly_mul1(const double u1[17], const double u2[17], const int p1, const int p2, double u3[17], int& p3) {
    for(int i = 0; i < 17; ++i)
        u3[i] = 0.0;
    double u1i;
    p3 = p1 + p2;
    
    for(int i1 = 0; i1 < p1 + 1; ++i1) {
        u1i = u1[i1];
        for(int i2 = 0; i2 < p2 + 1; ++i2) {
            int i3 = i1 + i2;
            u3[i3] += u1i * u2[i2];
        }
    }

    #ifdef TLC_LOG
    log("START ROUTINE \"poly_mul1\"");
    log(fmt::format("IN:: u1(:) = {}", TlcUtils::print_array(u1, 17)));
    log(fmt::format("IN:: u2(:) = {}", TlcUtils::print_array(u2, 17)));
    log(fmt::format("IN:: p1 = {}", p1));
    log(fmt::format("IN:: p2 = {}", p2));
    log(fmt::format("CHECK:: u3(:) = {}", TlcUtils::print_array(u3, 17)));
    log(fmt::format("CHECK:: p3 = {}", p3));
    log("END ROUTINE \"poly_mul1\"");
    #endif
}

void TLCMain::poly_sub1(const double u1[17], const double u2[17], const int p1, const int p2, double u3[17], int& p3) {
    for(int i = 0; i < 17; ++i)
        u3[i] = 0.0;
    p3 = std::max(p1, p2);

    for(int i = 0; i < p3 + 1; ++i) {
        if(i > p2)
            u3[i] = u1[i];
        else if(i > p1)
            u3[i] = -u2[i];
        else
            u3[i] = u1[i] - u2[i];
    }

    #ifdef TLC_LOG
    log("START ROUTINE \"poly_sub1\"");
    log(fmt::format("IN:: u1(:) = {}", TlcUtils::print_array(u1, 17)));
    log(fmt::format("IN:: u2(:) = {}", TlcUtils::print_array(u2, 17)));
    log(fmt::format("IN:: p1 = {}", p1));
    log(fmt::format("IN:: p2 = {}", p2));
    log(fmt::format("CHECK:: u3(:) = {}", TlcUtils::print_array(u3, 17)));
    log(fmt::format("CHECK:: p3 = {}", p3));
    log("END ROUTINE \"poly_sub1\"");
    #endif
}

void TLCMain::coord_from_poly_roots(int& n_soln, double roots[16], const double r_n1[3], const double r_a1[3], const double r_a3[3], const double r_c3[3]) {
    if(n_soln <= 0)
        return; // Better throw runtime_error

    double ex[3], ey[3], ez[3];
    cblas_dcopy(3, &b_a1a3[0], 1, &ex[0], 1);
    TlcUtils::cross(r_a1n1, ex, ez);
    cblas_dscal(3, 1/cblas_dnrm2(3, &ez[0], 1), &ez[0], 1);
    TlcUtils::cross(ez, ex, ey);

    #ifdef TLC_LOG
        log("START ROUTINE \"coord_pos1\"");
        log(fmt::format("IN:: r_n1 = {}", TlcUtils::print_array(r_n1, 3)));
        log(fmt::format("IN:: r_a1 = {}", TlcUtils::print_array(r_a1, 3)));
        log(fmt::format("IN:: r_a3 = {}", TlcUtils::print_array(r_a3, 3)));
        log(fmt::format("IN:: r_c3 = {}", TlcUtils::print_array(r_c3, 3)));
        log(fmt::format("CHECK:: ex = {}", TlcUtils::print_array(ex, 3)));
        log(fmt::format("CHECK:: ey = {}", TlcUtils::print_array(ey, 3)));
        log(fmt::format("CHECK:: ez = {}", TlcUtils::print_array(ez, 3)));
        log("END ROUTINE \"coord_pos1\"");
    #endif

    double b_a1a2[3], b_a3a2[3];
    for(int i = 0; i < 3; ++i) {
        b_a1a2[i] = -cos_alpha[0]*ex[i] + sin_alpha[0]*ey[i];
        b_a3a2[i] = cos_alpha[2]*ex[i] + sin_alpha[2]*ey[i];
    }
    double p_s[3][3], s1[3][3], s2[3][3], p_t[3][3], t1[3][3], t2[3][3];
    for(int i = 0; i < 3; ++i) {
        p_s[0][i] = -ex[i];
        s1[0][i] = ez[i];
        s2[0][i] = ey[i];
        p_t[0][i] = b_a1a2[i];
        t1[0][i] = ez[i] ;
        t2[0][i] = sin_alpha[0]*ex[i] + cos_alpha[0]*ey[i];

        p_s[1][i] = -b_a1a2[i];
        s1[1][i] = -ez[i];
        s2[1][i] = t2[0][i];
        p_t[1][i] = -b_a3a2[i];
        t1[1][i] = -ez[i];
        t2[1][i] = sin_alpha[2]*ex[i] - cos_alpha[2]*ey[i];

        p_s[2][i] = b_a3a2[i];
        s2[2][i] = t2[1][i];
        s1[2][i] = ez[i];
        p_t[2][i] = ex[i];
        t1[2][i] = ez[i];
        t2[2][i] = -ey[i];
    }

    double p_s_c[3][3], s1_s[3][3], s2_s[3][3], p_t_c[3][3], t1_s[3][3], t2_s[3][3];
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            p_s_c[i][j] = p_s[i][j]*cos_xi[i];
            s1_s[i][j] = s1[i][j]*sin_xi[i];
            s2_s[i][j] = s2[i][j]*sin_xi[i];
            p_t_c[i][j] = p_t[i][j]*cos_eta[i];
            t1_s[i][j] = t1[i][j]*sin_eta[i];
            t2_s[i][j] = t2[i][j]*sin_eta[i];
        }
    }

    #ifdef TLC_LOG
        log("START ROUTINE \"coord_pos2\"");
        for(int i = 0; i < 3; ++i) {
            log(fmt::format("CHECK:: p_s_c(:,{}) = {}", i + 1, TlcUtils::print_array(p_s_c[i], 3)));
            log(fmt::format("CHECK:: s1_s(:,{}) = {}", i + 1, TlcUtils::print_array(s1_s[i], 3)));
            log(fmt::format("CHECK:: s2_s(:,{}) = {}", i + 1, TlcUtils::print_array(s2_s[i], 3)));
            log(fmt::format("CHECK:: p_t_c(:,{}) = {}", i + 1, TlcUtils::print_array(p_t_c[i], 3)));
            log(fmt::format("CHECK:: t1_s(:,{}) = {}", i + 1, TlcUtils::print_array(t1_s[i], 3)));
            log(fmt::format("CHECK:: t2_s(:,{}) = {}", i + 1, TlcUtils::print_array(t2_s[i], 3)));
        }
        log("END ROUTINE \"coord_pos2\"");
    #endif

    double r_tmp[3];
    for(int i = 0; i < 3; ++i)
        r_tmp[i] = (r_a1n1[i]/len_na[0] - p_s_c[0][i])/sin_xi[0];
    double angle = calc_bnd_ang(s1[0], r_tmp);

    // using Eigen::Vector3d;
    // Eigen::Map<const Eigen::Vector3d> r_tmp_v(r_tmp);
    // Eigen::Map<const Eigen::Vector3d> s21_v(s2[0]);
    // double sig1_init = std::copysign(angle, r_tmp_v.dot(s21_v));
    double sig1_init = std::copysign(angle, cblas_ddot(3, &r_tmp[0], 1, &s2[0][0], 1));

    double r_n[3][3], r_a[3][3], r_c[3][3], r0[3];
    for(unsigned int i = 0; i < 3; ++i) {
        r_a[0][i] = r_a1[i];
        r_a[1][i] = r_a1[i] + len_aa[1] * b_a1a2[i];
        r_a[2][i] = r_a3[i];
        r0[i] = r_a1[i];
    }

    double half_tan[3], cos_tau[3], sin_tau[3];

    double minus_ex[3];
    cblas_dcopy(3, &ex[0], 1, &minus_ex[0], 1);
    cblas_dscal(3, -1, &minus_ex[0], 1);
    unsigned int soln_idx = 0;
    for(unsigned int i_soln = 0; i_soln < n_soln; ++i_soln) {
        half_tan[2] = roots[i_soln];
        half_tan[1] = calc_t2(half_tan[2]);
        half_tan[0] = calc_t1(half_tan[2], half_tan[1]);
        for(int i = 0; i < 3; ++i) {
            double ht = half_tan[i];
            double tmp = 1.0 + ht*ht;
            int j = 0; // if i == 2
            if(i < 2) // if i == 0 or 1
                j = i + 1;
            cos_tau[j] = (1.0 - ht*ht)/tmp;
            sin_tau[j] = 2.0*ht/tmp;
        }

        #ifdef TLC_LOG
            log("START ROUTINE \"coord_pos21\"");
            log(fmt::format("IN:: i_soln = {}", i_soln + 1));
            log(fmt::format("CHECK:: half_tan = {}", TlcUtils::print_array(half_tan)));
            log(fmt::format("CHECK:: cos_tau = {}", TlcUtils::print_array(cos_tau)));
            log(fmt::format("CHECK:: sin_tau = {}", TlcUtils::print_array(sin_tau)));
            log("END ROUTINE \"coord_pos21\"");
        #endif
        
        double cos_sig[3], sin_sig[3];
        for(int i = 0; i < 3; ++i) {
            int j = 0;
            if(i < 2)
                j = i + 1;
            cos_sig[j] = cos_delta[i]*cos_tau[j] + sin_delta[i]*sin_tau[j];
            sin_sig[j] = sin_delta[i]*cos_tau[j] - cos_delta[i]*sin_tau[j];
        }

        double temp; // r_s[3], r_t[3]
        for(int i = 0; i < 3; ++i) {
            int k = 0;
            if (i < 2)
                k = i + 1;
            for(int j = 0; j < 3; ++j) {
                // r_s[j] = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
                // r_t[j] = p_t_c[i][j] + cos_tau[k]*t1_s[i][j] + sin_tau[k]*t2_s[i][j];
                temp = p_s_c[i][j] + cos_sig[i]*s1_s[i][j] + sin_sig[i]*s2_s[i][j];
                r_n[i][j] = temp * len_na[i] + r_a[i][j];
                temp = p_t_c[i][j] + cos_tau[k]*t1_s[i][j] + sin_tau[k]*t2_s[i][j];
                r_c[i][j] = temp * len_ac[i] + r_a[i][j];
            }
        }

        #ifdef TLC_LOG
            log("START ROUTINE \"coord_pos3\"");
            log(fmt::format("IN:: i_soln = {}", i_soln + 1));
            log(fmt::format("IN:: len_na = {}", TlcUtils::print_array(len_na)));
            log(fmt::format("IN:: len_ac = {}", TlcUtils::print_array(len_ac)));
            log(fmt::format("IN:: cos_sig = {}", TlcUtils::print_array(cos_sig)));
            log(fmt::format("IN:: sin_sig = {}", TlcUtils::print_array(sin_sig)));
            log(fmt::format("IN:: cos_tau = {}", TlcUtils::print_array(cos_tau)));
            log(fmt::format("IN:: sin_tau = {}", TlcUtils::print_array(sin_tau)));
            for(int i = 0; i < 3; ++i) {
                log(fmt::format("CHECK:: r_n(:,{}) = {}", i + 1, TlcUtils::print_array(r_n[i], 3)));
                log(fmt::format("CHECK:: r_c(:,{}) = {}", i + 1, TlcUtils::print_array(r_c[i], 3)));
            }
            log("END ROUTINE \"coord_pos3\"");
        #endif

        double sig1 = atan2(sin_sig[0], cos_sig[0]);
        double p[4], Us[3][3];
        quaternion(minus_ex, -(sig1 - sig1_init) * 0.25, p);
        rotation_matrix(p, Us);

        #ifdef TLC_LOG
            log("START ROUTINE \"coord_pos4\"");
            log(fmt::format("IN:: i_soln = {}", i_soln + 1));
            log(fmt::format("CHECK:: sig1 = {}", sig1));
            log(fmt::format("CHECK:: sig1_init = {}", sig1_init));
            log(fmt::format("CHECK:: sin_sig(1) = {}", sin_sig[0]));
            log(fmt::format("CHECK:: cos_sig(1) = {}", cos_sig[0]));
            log(fmt::format("CHECK:: p = {}", TlcUtils::print_array(p)));
            for(int i = 0; i < 3; ++i) {
                log(fmt::format("CHECK:: Us(:,{}) = {}", i + 1, TlcUtils::print_array(Us[i], 3)));
            }
            log("END ROUTINE \"coord_pos4\"");
        #endif

        cblas_dcopy(3, &r_n1[0], 1, &r_soln_n[i_soln][0][0], 1);
        cblas_dcopy(3, &r_a1[0], 1, &r_soln_a[i_soln][0][0], 1);

        double r0_add[3];
        cblas_dcopy(3, &r0[0], 1, &r0_add[0], 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, -1.0, &Us[0][0], 3, &r0[0], 1, 1.0, &r0_add[0], 1);
        
        cblas_dcopy(3, &r0_add[0], 1, &r_soln_c[soln_idx][0][0], 1);
        cblas_dcopy(3, &r0_add[0], 1, &r_soln_n[soln_idx][1][0], 1);
        cblas_dcopy(3, &r0_add[0], 1, &r_soln_a[soln_idx][1][0], 1);
        cblas_dcopy(3, &r0_add[0], 1, &r_soln_c[soln_idx][1][0], 1);
        cblas_dcopy(3, &r0_add[0], 1, &r_soln_n[soln_idx][2][0], 1);

        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &r_c[0][0], 1, 1.0, &r_soln_c[soln_idx][0][0], 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &r_n[1][0], 1, 1.0, &r_soln_n[soln_idx][1][0], 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &r_a[1][0], 1, 1.0, &r_soln_a[soln_idx][1][0], 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &r_c[1][0], 1, 1.0, &r_soln_c[soln_idx][1][0], 1);
        cblas_dgemv(CblasColMajor, CblasNoTrans, 3, 3, 1.0, &Us[0][0], 3, &r_n[2][0], 1, 1.0, &r_soln_n[soln_idx][2][0], 1);

        cblas_dcopy(3, &r_a3[0], 1, &r_soln_a[soln_idx][2][0], 1);
        cblas_dcopy(3, &r_c3[0], 1, &r_soln_c[soln_idx][2][0], 1);

        #ifdef TLC_LOG
            log("START ROUTINE \"coord_final\"");
            log(fmt::format("IN:: i_soln = {}", i_soln + 1));
            log(fmt::format("CHECK:: r_n1 = {}", TlcUtils::print_array(r_n1, 3)));
            log(fmt::format("CHECK:: r_a1 = {}", TlcUtils::print_array(r_a1, 3)));
            log(fmt::format("CHECK:: r_a3 = {}", TlcUtils::print_array(r_a3, 3)));
            log(fmt::format("CHECK:: r_c3 = {}", TlcUtils::print_array(r_c3, 3)));
            for(int i = 0; i < 3; ++i) {
                log(fmt::format("CHECK:: r_soln_a(:,{}) = {}", i + 1, TlcUtils::print_array(r_soln_a[soln_idx][i], 3)));
                log(fmt::format("CHECK:: r_soln_n(:,{}) = {}", i + 1, TlcUtils::print_array(r_soln_n[soln_idx][i], 3)));
                log(fmt::format("CHECK:: r_soln_c(:,{}) = {}", i + 1, TlcUtils::print_array(r_soln_c[soln_idx][i], 3)));
            }
            log("END ROUTINE \"coord_final\"");
        #endif

        const double vangleA = get_vangle(r_soln_n[soln_idx][0], r_soln_a[soln_idx][0], r_soln_c[soln_idx][0]);
        const double vangleB = get_vangle(r_soln_n[soln_idx][1], r_soln_a[soln_idx][1], r_soln_c[soln_idx][1]);
        const double vangleC = get_vangle(r_soln_n[soln_idx][2], r_soln_a[soln_idx][2], r_soln_c[soln_idx][2]);
        if ((abs(vangleA - this->theta[0]) < THRESHOLD) &&
            (abs(vangleB - this->theta[1]) < THRESHOLD) &&
            (abs(vangleC - this->theta[2]) < THRESHOLD)) {
            soln_idx++;
        } else {
            some_removed = true;
        }
    }
    n_soln = soln_idx;

    #ifdef TLC_LOG
        using FMatrix = boost::numeric::ublas::fixed_matrix<double, 9, 3>;
        FMatrix conf(9, 3);
        for (int i = 0; i < n_soln; ++i) {
            solution_to_buffer(&conf(0, 0), i);

            log("START ROUTINE \"xyzcheck\"");
            log(fmt::format(" IN:: bonds = {}", TlcUtils::print_array(m_bonds)));
            log(fmt::format(" IN:: vangles = {}", TlcUtils::print_array(m_vangles)));
            log(fmt::format(" IN:: tangles = {}", TlcUtils::print_array(m_tangles)));
            log(fmt::format(" IN:: sol_idx = {}", i));
            log(fmt::format(" CHECK:: xyz = {}", TlcUtils::repr_matrix(conf)));
            log("END ROUTINE \"xyzcheck\"");
        }
    #endif
}

inline double TLCMain::get_vangle(double (&a)[3], double (&b)[3], double (&c)[3]) noexcept {
    double dirA[3], dirB[3];
    cblas_dcopy(3, &a[0], 1, &dirA[0], 1);
    cblas_dcopy(3, &c[0], 1, &dirB[0], 1);
    cblas_daxpy(3, -1.0, &b[0], 1, &dirA[0], 1);
    cblas_daxpy(3, -1.0, &b[0], 1, &dirB[0], 1);
    cblas_dscal(3, 1 / cblas_dnrm2(3, &dirA[0], 1), &dirA[0], 1);
    cblas_dscal(3, 1 / cblas_dnrm2(3, &dirB[0], 1), &dirB[0], 1);
    return calc_bnd_ang(dirA, dirB);
}

double TLCMain::calc_t2(double t0) {
    double t0_2 = t0*t0;
    double t0_3 = t0_2*t0;
    double t0_4 = t0_3*t0;

    double A0 = Q[0][0] + Q[0][1]*t0 + Q[0][2]*t0_2 + Q[0][3]*t0_3 + Q[0][4]*t0_4;
    double A1 = Q[1][0] + Q[1][1]*t0 + Q[1][2]*t0_2 + Q[1][3]*t0_3 + Q[1][4]*t0_4;
    double A2 = Q[2][0] + Q[2][1]*t0 + Q[2][2]*t0_2 + Q[2][3]*t0_3 + Q[2][4]*t0_4;
    double A3 = Q[3][0] + Q[3][1]*t0 + Q[3][2]*t0_2 + Q[3][3]*t0_3 + Q[3][4]*t0_4;
    double A4 = Q[4][0] + Q[4][1]*t0 + Q[4][2]*t0_2 + Q[4][3]*t0_3 + Q[4][4]*t0_4;

    double B0 = R[0][0] + R[0][1]*t0 + R[0][2]*t0_2;
    double B1 = R[1][0] + R[1][1]*t0 + R[1][2]*t0_2;
    double B2 = R[2][0] + R[2][1]*t0 + R[2][2]*t0_2;

    double B2_2 = B2*B2;
    double B2_3 = B2_2*B2;

    double K0 = A2*B2 - A4*B0;
    double K1 = A3*B2 - A4*B1;
    double K2 = A1*B2_2 - K1*B0;
    double K3 = K0*B2 - K1*B1;
    double res = (K3*B0 - A0*B2_3)/(K2*B2 - K3*B1);
    // if (abs(K2*B2 - K3*B1) < 1E-8)
    //     // throw std::runtime_error(fmt::format("res = {}; den = {}", res, (K2*B2 - K3*B1)));
    // } else {
        #ifdef TLC_LOG
            log("START ROUTINE \"calc_t2\"");
            log(fmt::format("IN:: t0 = {}", t0));
            log(fmt::format("CHECK:: calc_t2 = {}", res));
            log("END ROUTINE \"calc_t2\"");
        #endif
    // }
    return res;
}

double TLCMain::calc_t1(double t0, double t2) {
    double t0_2 = t0*t0;
    double t2_2 = t2*t2;

    double U11 = C[0][0][0] + C[0][0][1]*t0 + C[0][0][2]*t0_2;
    double U12 = C[1][0][0] + C[1][0][1]*t0 + C[1][0][2]*t0_2;
    double U13 = C[2][0][0] + C[2][0][1]*t0 + C[2][0][2]*t0_2;
    double U31 = C[0][1][0] + C[0][1][1]*t2 + C[0][1][2]*t2_2;
    double U32 = C[1][1][0] + C[1][1][1]*t2 + C[1][1][2]*t2_2;
    double U33 = C[2][1][0] + C[2][1][1]*t2 + C[2][1][2]*t2_2;

    double res = (U31*U13-U11*U33)/(U12*U33-U13*U32);
    #ifdef TLC_LOG
        log("START ROUTINE \"calc_t1\"");
        log(fmt::format("IN:: t0 = {}", t0));
        log(fmt::format("IN:: t2 = {}", t2));
        log(fmt::format("CHECK:: calc_t1 = {}", res));
        log("END ROUTINE \"calc_t1\"");
    #endif
    return res;
}

bool TLCMain::test_two_cone_existence_soln(const double tt, const double kx, const double et, const double ap) {
    int n_soln = MAX_SOLN;
    if (abs(ap - tt) > kx + et) {
        n_soln = 0;
    }

    #ifdef TLC_LOG
        log("START ROUTINE \"test_two_cone_existence_soln\"");
        log(fmt::format("IN:: tt = {}", tt));
        log(fmt::format("IN:: kx = {}", kx));
        log(fmt::format("IN:: et = {}", et));
        log(fmt::format("IN:: ap = {}", ap));
        log(fmt::format("CHECK:: n_soln = {}", n_soln));
        log("END ROUTINE \"test_two_cone_existence_soln\"");
    #endif
    return n_soln > 0;
}

void TLCMain::quaternion(const double axis[], const double& ang, double res[]) {
    double tan_w = tan(ang);
    double tan_sqr = pow(tan_w, 2);
    double cosine = (1.0 - tan_sqr)/(1.0 + tan_sqr);
    double sine = 2.0*tan_w/(1.0 + tan_sqr);
    res[0] = cosine;
    for(int i = 0; i < 3; ++i)
        res[i + 1] = axis[i] * sine;

    #ifdef TLC_LOG
        log("START ROUTINE \"quaternion\"");
        log(fmt::format("IN:: axis = {}", TlcUtils::print_array(axis, 3)));
        log(fmt::format("IN:: quarter_ang = {}", ang));
        log(fmt::format("CHECK:: p = {}", TlcUtils::print_array(res, 4)));
        log("END ROUTINE \"quaternion\"");
    #endif
}

void TLCMain::rotation_matrix(const double q[], double U[][3]) {
    double b[4];
    for(int i = 0; i < 4; ++i)
        b[i] = 2 * q[i];
    
    double Q[4][4];
    Q[0][0] = b[0]*q[0]-1.0;
    Q[0][2] = b[0]*q[2];
    Q[0][3] = b[0]*q[3];
    Q[1][1] = b[1]*q[1];
    Q[1][2] = b[1]*q[2];
    Q[1][3] = b[1]*q[3];
    Q[0][1] = b[0]*q[1];
    Q[2][2] = b[2]*q[2];
    Q[2][3] = b[2]*q[3];
    Q[3][3] = b[3]*q[3];

    U[0][0] = Q[0][0]+Q[1][1]; U[1][0] = Q[1][2]-Q[0][3]; U[2][0] = Q[1][3]+Q[0][2];
    U[0][1] = Q[1][2]+Q[0][3]; U[1][1] = Q[0][0]+Q[2][2]; U[2][1] = Q[2][3]-Q[0][1];
    U[0][2] = Q[1][3]-Q[0][2]; U[1][2] = Q[2][3]+Q[0][1]; U[2][2] = Q[0][0]+Q[3][3];

    #ifdef TLC_LOG
        log("START ROUTINE \"rotation_matrix\"");
        log(fmt::format("IN:: q = {}", TlcUtils::print_array(q, 4)));
        log(fmt::format("CHECK:: U(:, 1) = {}", TlcUtils::print_array(U[0], 3)));
        log(fmt::format("CHECK:: U(:, 2) = {}", TlcUtils::print_array(U[1], 3)));
        log(fmt::format("CHECK:: U(:, 3) = {}", TlcUtils::print_array(U[2], 3)));
        log("END ROUTINE \"rotation_matrix\"");
    #endif
}

double TLCMain::calc_bnd_ang(const double r1[], const double r2[]) {
    double res = cblas_ddot(3, &r1[0], 1, &r2[0], 1);
    res = std::copysign(std::min(std::abs(res), 1.0), res);
    res = acos(res);

    // #ifdef TLC_LOG
    //     log("START ROUTINE \"calc_bnd_ang\"");
    //     log(fmt::format("IN:: r1 = {}", TlcUtils::print_array(r1, 3)));
    //     log(fmt::format("IN:: r2 = {}", TlcUtils::print_array(r2, 3)));
    //     log(fmt::format("CHECK:: angle = {}", res));
    //     log("END ROUTINE \"calc_bnd_ang\"");
    // #endif
    return res;
}

double TLCMain::calc_dih_ang(const double r1[], const double r2[], const double r3[]) {
    double p[3], q[3], s[3];
    TlcUtils::cross(r1, r2, p);
    TlcUtils::cross(r2, r3, q);
    TlcUtils::cross(r3, r1, s);
    
    double res = cblas_ddot(3, &p[0], 1, &q[0], 1);
    res /= cblas_dnrm2(3, &q[0], 1);
    res /= cblas_dnrm2(3, &p[0], 1);
    
    res = std::copysign(std::min(std::abs(res), 1.0), res);
    res = std::copysign(acos(res), cblas_ddot(3, &s[0], 1, &r2[0], 1));
    
    #ifdef TLC_LOG
        log("START ROUTINE \"calc_dih_ang\"");
        log(fmt::format("IN:: r1 = {}", TlcUtils::print_array(r1, 3)));
        log(fmt::format("IN:: r2 = {}", TlcUtils::print_array(r2, 3)));
        log(fmt::format("IN:: r3 = {}", TlcUtils::print_array(r3, 3)));
        log(fmt::format("CHECK:: angle = {}", res));
        log("END ROUTINE \"calc_dih_ang\"");
    #endif
    return res;
}

void TLCMain::set_solution(py::array_t<double> solution, const int& sol_idx) {
    py::buffer_info buffer = solution.request();
    double *buffer_raw = static_cast<double*>(buffer.ptr);

    if (buffer.ndim != 2)
        throw std::runtime_error("Input should be 2-D NumPy array");

    if ((solution.shape()[0] != 9) || (solution.shape()[1] != 3))
        throw std::runtime_error(fmt::format("Dimensions of np.array are not right: ({}; {}) instead of (9; 3)", solution.shape()[0], solution.shape()[1]));
    solution_to_buffer(buffer_raw, sol_idx);
}
