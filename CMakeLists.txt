cmake_minimum_required(VERSION 3.4...3.18)
project(ringo)


set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)
set(CMAKE_CXX_STANDARD 17)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)

add_subdirectory(pybind11)
add_subdirectory(fmt)
add_subdirectory(eigen)
add_subdirectory(json)

link_directories(fmt/build)
include_directories(fmt/include)
include_directories(eigen/Eigen)
include_directories(json/include)

# find_package(GSL REQUIRED)
# include_directories(${GSL_INCLUDE_DIRS})

# include(CheckIncludeFileCXX)
# CHECK_INCLUDE_FILE_CXX(openblas/cblas.h HAVE_OPENBLAS_CBLAS_H)
# if (NOT HAVE_OPENBLAS_CBLAS_H)
#     message(FATAL_ERROR "ERROR: <openblas/cblas.h> is unavailable!!!")
# endif()

if (DEBUG)
    # message(FATAL_ERROR "DEBUG flag is not allowed in this project!")
    add_definitions(-DKDMOL_DEBUG=1)
    if (UNITTEST)
        add_definitions(-DKDMOL_LOG=1)
        add_definitions(-DVALIDATION=1)
    endif()
    if (TLC_UNITTEST)
        add_definitions(-DEXT_LOG=1)
        # add_definitions(-DTLC_LOG=1)
    endif()
endif()
if (ENDTOEND)
    # message(FATAL_ERROR "ENDTOEND flag is not allowed in this project!")
    add_definitions(-DVALIDATION=1)
endif()
if (NOT ENDTOENDWEAK)
    # message(FATAL_ERROR "ENDTOENDWEAK flag must be defined in this project!")
endif()
if (ENDTOENDWEAK)
add_definitions(-DWEAK_VALIDATION=1)
endif()

if (OVERLAP_DETECTION)
    add_definitions(-DOVERLAP_DETECTION=1)
endif()
if (OVERLAP_DETECTION_FINAL)
    add_definitions(-DOVERLAP_DETECTION_FINAL=1)
endif()

if (OVERLAP_MULTIPLIER)
    add_definitions(-DOVERLAP_MULTIPLIER=${OVERLAP_MULTIPLIER})
endif()

if (RMSD_FILTER)
    add_definitions(-DRMSD_FILTER=1)
    add_definitions(-DBUILD_PYXYZ=1)
    option(BUILD_PYXYZ "Building Pyxyz" ON)
endif()

if (USE_CHECKPOINTS)
    add_definitions(-DUSE_CHECKPOINTS=1)
    endif()
add_definitions(-DRINGO_BUILDFLAGS="\"${BUILDFLAGS}\"")

find_package (Eigen3 REQUIRED NO_MODULE)
find_package(Boost COMPONENTS system filesystem REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")

set(PROJ_SRC_FILES src/main.cpp src/molecule.cpp src/utils.cpp
        src/pytlc/sturm.c src/pytlc/ringext.cpp src/pytlc/tlcmain.cpp
        src/pytlc/utils.cpp src/pytlc/rigidfrag.cpp src/solvers.cpp)
if (BUILD_PYXYZ)
    list(APPEND PROJ_SRC_FILES src/pyxyz/molproxy.cpp src/pyxyz/utils.cpp src/pyxyz/confpool.cpp)
endif()
pybind11_add_module(ringo_base ${PROJ_SRC_FILES})
target_include_directories(ringo_base PRIVATE src/)
target_compile_options(ringo_base PRIVATE -pedantic)
target_link_options(ringo_base PRIVATE -pedantic)

# if (NOT DEBUG)
    target_compile_options(ringo_base PRIVATE -pedantic -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG) # -fsanitize=address
    target_link_options(ringo_base PRIVATE -pedantic -O3 -D NDEBUG -DBOOST_UBLAS_NDEBUG) #  -fsanitize=address
# endif()
target_link_libraries(ringo_base PRIVATE fmt::fmt Eigen3::Eigen -lgslcblas) # ${OPENBLAS_LINKING}

