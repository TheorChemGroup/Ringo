#ifndef PYXYZ_UTILS_H
#define PYXYZ_UTILS_H

#include "../utils.h"

# define DOUBLE_THRESHOLD 0.9999999999


using BoostRow = boost::numeric::ublas::matrix_row<boost::numeric::ublas::matrix<double>>;

namespace PyxyzUtils {
    template <typename C, typename P> 
    void erase_remove_if(C& c, P predicate) {
        c.erase(std::remove_if(c.begin(), c.end(), predicate), c.end());
    }

    struct InvalidChar {
        bool operator()(char c) const {
            return !isprint(static_cast<unsigned char>(c));
        }
    };

    template <class S>
    std::vector<std::string> readlines_rofl(S filename) {
        #ifdef WIN_PLATFORM
            std::ifstream infile(ConvertFromUtf16ToUtf8(py::handle(filename).cast<std::wstring>())); // py::handle(filename).cast<std::wstring>().c_str(), std::ios::in | std::ios::binary
        #else
            std::ifstream infile(filename);
        #endif
        std::string str;
        std::vector<std::string> reslines;

        if(!infile.is_open())
            throw std::runtime_error("File not found " + filename);

        while (std::getline(infile, str))
        {
            // if(str.size() > 0) {
			erase_remove_if(str, InvalidChar());
			reslines.push_back(str);
            // }
        }
        return reslines;
    }

    #ifdef WIN_PLATFORM
        std::string ConvertFromUtf16ToUtf8(const std::wstring& wstr);
    #endif

    template<class MType>
    void make_centered(MType& xyzs) {
        using boost::numeric::ublas::matrix_column;
        using boost::numeric::ublas::sum;
        auto size = xyzs.size1();
        for(int i = 0; i < 3; ++i) {
            matrix_column<MType> c(xyzs, i);
            double mean = sum(c)/size;
            for(auto& coord : c)
                coord -= mean;
        }
    }
    
    template<class MType>
    void make_centered(MType& xyzs, double* center) {
        using boost::numeric::ublas::matrix_column;
        using boost::numeric::ublas::sum;
        auto size = xyzs.size1();
        for(int i = 0; i < 3; ++i) {
            matrix_column<MType> c(xyzs, i);
            double mean = sum(c)/size;
            center[i] = mean;
            for(auto& coord : c)
                coord -= mean;
        }
    }

    template<class MType>
    py::array_t<double> ublas_matrix_to_np(const MType& m) {
        const size_t size = m.size1() * m.size2();
        double *data_array = new double[size];
        for (unsigned int i = 0; i < m.size1(); ++i)
            for (unsigned int j = 0; j < m.size2(); ++j)
                data_array[m.size2() * i + j] = m(i, j);
        py::capsule free_when_done(data_array, [](void *f) {
            double *foo = reinterpret_cast<double *>(f);
            delete[] foo;
        });
        
        return py::array_t<double>(
            std::vector<long unsigned int>{m.size1(), m.size2()},
            data_array,
            free_when_done);
    }
    
    template<class VType>
    py::array_t<double> ublas_vector_to_np(const VType& v) {
        const size_t size = v.size();
        double *data_array = new double[size];
        for (unsigned int i = 0; i < v.size(); ++i)
            data_array[i] = v(i);
        py::capsule free_when_done(data_array, [](void *f) {
            double *foo = reinterpret_cast<double *>(f);
            delete[] foo;
        });
        
        return py::array_t<double>(
            {v.size()}, // shape
            {8}, // C-style contiguous strides for double
            data_array, // the data pointer
            free_when_done);
    }

    template<class MType>
    void fill_xyz_matrix(MType& xyz_matr, const py::list& py_xyzs) {
        for(int i = 0; i < xyz_matr.size1(); ++i) {
            auto xyz = py_xyzs[i].cast< std::array<double, 3> >();
            xyz_matr(i, 0) = xyz[0];
            xyz_matr(i, 1) = xyz[1];
            xyz_matr(i, 2) = xyz[2];
        }
    }

    template<class MType>
    inline double det3x3(const MType& m) {
        return m(0, 0)*m(1, 1)*m(2, 2)+
               m(1, 0)*m(2, 1)*m(0, 2)+
               m(2, 0)*m(1, 2)*m(0, 1)-
               m(2, 0)*m(1, 1)*m(0, 2)-
               m(1, 0)*m(0, 1)*m(2, 2)-
               m(2, 1)*m(1, 2)*m(0, 0);
    }

    template <typename T>
    std::vector<std::size_t> sort_permutation(
        const std::vector<T>& vec, bool ascend)
    {
        std::vector<std::size_t> p(vec.size());
        std::iota(p.begin(), p.end(), 0);
        if (ascend)
            std::sort(p.begin(), p.end(),
                [&](std::size_t i, std::size_t j){ return vec[i] < vec[j]; });
        else
            std::sort(p.begin(), p.end(),
                [&](std::size_t i, std::size_t j){ return vec[i] > vec[j]; });
        return p;
    }

    template <typename T>
    void apply_permutation_in_place(
        std::vector<T>& vec,
        const std::vector<std::size_t>& p)
    {
        std::vector<bool> done(vec.size());
        for (std::size_t i = 0; i < vec.size(); ++i)
        {
            if (done[i])
            {
                continue;
            }
            done[i] = true;
            std::size_t prev_j = i;
            std::size_t j = p[i];
            while (i != j)
            {
                std::swap(vec[prev_j], vec[j]);
                done[j] = true;
                prev_j = j;
                j = p[j];
            }
        }
    }

    template <typename MType>
    std::string repr_matrix(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.size1 (); ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m.size2 (); ++j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_eigen_matrix(const MType& m) {
        py::list res;
        for (unsigned int i = 0; i < m.rows (); ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m.cols (); ++j)
                temp.append(py::cast(m(i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_gsl_matrix(const MType* m) {
        py::list res;
        for (unsigned int i = 0; i < m->size1; ++i) {
            py::list temp;
            for (unsigned int j = 0; j < m->size2; ++j)
                temp.append(py::cast(gsl_matrix_get(m, i, j)));
            res.append(temp);
        }
        return py::repr(res).cast<std::string>();
    }
    
    template <typename MType>
    std::string repr_gsl_vector(const MType* m) {
        py::list res;
        for (unsigned int i = 0; i < m->size; ++i)
            res.append(py::cast(gsl_vector_get(m, i)));
        return py::repr(res).cast<std::string>();
    }
    
    template <typename VType>
    std::string repr_vector(const VType& v) {
        py::list res;
        for (unsigned int i = 0; i < v.size(); ++ i)
            res.append(v[i]);
        return py::repr(res).cast<std::string>();
    }
    
    template <typename Type>
    std::string repr_sorted(const Type& v) {
        py::list res;
        for (const auto& item : v)
            res.append(item);
        res.attr("sort")();
        return py::repr(res).cast<std::string>();
    }
    std::string repr_matrix_buffer(const double* buf, const int& size1, const int& size2);


    const std::unordered_map<std::string, int> NAMES_ELEMENT = {
        {"H", 1},
        {"He", 2},
        {"Li", 3},
        {"Be", 4},
        {"B", 5},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9},
        {"Ne", 10},
        {"Na", 11},
        {"Mg", 12},
        {"Al", 13},
        {"Si", 14},
        {"P", 15},
        {"S", 16},
        {"Cl", 17},
        {"Ar", 18},
        {"K", 19},
        {"Ca", 20},
        {"Sc", 21},
        {"Ti", 22},
        {"V", 23},
        {"Cr", 24},
        {"Mn", 25},
        {"Fe", 26},
        {"Co", 27},
        {"Ni", 28},
        {"Cu", 29},
        {"Zn", 30},
        {"Ga", 31},
        {"Ge", 32},
        {"As", 33},
        {"Se", 34},
        {"Br", 35},
        {"Kr", 36},
        {"Rb", 37},
        {"Sr", 38},
        {"Y", 39},
        {"Zr", 40},
        {"Nb", 41},
        {"Mo", 42},
        {"Tc", 43},
        {"Ru", 44},
        {"Rh", 45},
        {"Pd", 46},
        {"Ag", 47},
        {"Cd", 48},
        {"In", 49},
        {"Sn", 50},
        {"Sb", 51},
        {"Te", 52},
        {"I", 53},
        {"Xe", 54},
        {"Cs", 55},
        {"Ba", 56},
        {"La", 57},
        {"Ce", 58},
        {"Pr", 59},
        {"Nd", 60},
        {"Pm", 61},
        {"Sm", 62},
        {"Eu", 63},
        {"Gd", 64},
        {"Tb", 65},
        {"Dy", 66},
        {"Ho", 67},
        {"Er", 68},
        {"Tm", 69},
        {"Yb", 70},
        {"Lu", 71},
        {"Hf", 72},
        {"Ta", 73},
        {"W", 74},
        {"Re", 75},
        {"Os", 76},
        {"Ir", 77},
        {"Pt", 78},
        {"Au", 79},
        {"Hg", 80},
        {"Tl", 81},
        {"Pb", 82},
        {"Bi", 83},
        {"Po", 84},
        {"At", 85},
        {"Rn", 86},
        {"Fr", 87},
        {"Ra", 88},
        {"Ac", 89},
        {"Th", 90},
        {"Pa", 91},
        {"U", 92},
        {"Np", 93},
        {"Pu", 94},
        {"Am", 95},
        {"Cm", 96},
        {"Bk", 97},
        {"Cf", 98},
        {"Es", 99},
        {"Fm", 100},
        {"Md", 101},
        {"No", 102},
        {"Lr", 103},
        {"Rf", 104},
        {"Db", 105},
        {"Sg", 106},
        {"Bh", 107},
        {"Hs", 108},
        {"Mt", 109},
        {"Ds", 110},
        {"Rg", 111},
        {"Cn", 112},
        {"Uuq", 114},
        {"Uuh", 116}
    };

    const std::vector<int> generate_atom_ints(std::vector<std::string> inp);
}

#endif // PYXYZ_UTILS_H
