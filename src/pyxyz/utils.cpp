// #include <iostream>
// #include <fstream>
// #include <math.h>
// #include <numeric>

#ifdef WIN_PLATFORM
    #include <windows.h>
    #include <stdio.h>
    #include <wchar.h>
    #include <string.h>
#endif

// #include <pybind11/pybind11.h>
// #include <pybind11/stl.h>
// #include <pybind11/numpy.h>

// #include <boost/numeric/ublas/matrix.hpp>
// #include <boost/numeric/ublas/matrix_proxy.hpp>

// #include <Eigen/Dense>
// #include <gsl/gsl_cblas.h>
// #include <fmt/xchar.h> // <<<<<<<<<<<<<<<<<<<

#include "utils.h"

namespace PyxyzUtils {
    #ifdef WIN_PLATFORM
        std::string ConvertFromUtf16ToUtf8(const std::wstring& wstr) {
            std::string convertedString;
            int requiredSize = WideCharToMultiByte(CP_ACP, 0, wstr.c_str(), -1, 0, 0, 0, 0);
            if (requiredSize > 0) {
                std::vector<char> buffer(requiredSize);
                WideCharToMultiByte(CP_ACP, 0, wstr.c_str(), -1, &buffer[0], requiredSize, 0, 0);
                convertedString.assign(buffer.begin(), buffer.end() - 1);
            } 
            return convertedString;
        }
    #endif
    
    // template <class S>
    // std::vector<std::string> readlines(S filename) {
    //     fmt::print("Prepare to suck {}\n", filename.cast<std::string>()); std::cout << std::endl;
    //     #ifdef WIN_PLATFORM
    //         std::ifstream infile(ConvertFromUtf16ToUtf8(filename.cast<std::wstring>())); // filename.cast<std::wstring>().c_str(), std::ios::in | std::ios::binary
    //     #else
    //         std::ifstream infile(filename.cast<std::string>());
    //     #endif
    //     fmt::print("Passed {}\n", filename.cast<std::string>()); std::cout << std::endl;
    //     std::string str;
    //     std::vector<std::string> reslines;

    //     if(!infile)
    //         throw std::runtime_error("File not found " + filename.cast<std::string>());
    //         // throw std::runtime_error(fmt::format(L"File {} not found", filename));

    //     while (std::getline(infile, str))
    //     {
    //         // if(str.size() > 0) {
	// 		erase_remove_if(str, InvalidChar());
	// 		reslines.push_back(str);
    //         // }
    //     }
    //     return reslines;
    // }

    const std::vector<int> generate_atom_ints(std::vector<std::string> inp) {
        std::vector<int> res;
        res.reserve(inp.size());
        for(const auto& item : inp)
            res.push_back(NAMES_ELEMENT.at(item));
        return res;
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
}