#ifndef SIDEHEADERS_H
#define SIDEHEADERS_H

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <numeric>
#include <fstream>
#include <algorithm>
#include <memory>
#include <tuple>
#include <string>
#include <map>
#include <initializer_list>
#include <cmath>
#include <stdexcept>   // for std::exception
#include <cassert>
#include <chrono>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/type_index.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <random>
#include "pcg/pcg_random.hpp"

#include "fmt/core.h"
#include "fmt/args.h"

#include <Eigen/Dense>
// #include <openblas/cblas.h>
#include <gsl/gsl_cblas.h>
#include <nlohmann/json.hpp>

// For parallel runs:
#include <thread>
#include <mutex>

#endif
