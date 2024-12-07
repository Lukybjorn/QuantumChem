cmake_minimum_required(VERSION 3.28)
project(QuantumChem LANGUAGES CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

include(FetchContent)
# Fetch Eigen
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
    GIT_SHALLOW ON
)
FetchContent_MakeAvailable(Eigen)



add_executable(main src/main.cpp
        src/libmath/overlap.h
        src/libmath/overlap.cpp
        src/libmath/Gaussian.h
        src/libmath/kinetic.cpp
        src/libmath/kinetic.h
        src/libmath/nucAttraction.cpp
        src/libmath/nucAttraction.h
        src/libmath/twoElectronRepulsion.cpp
        src/libmath/twoElectronRepulsion.h
        src/libgeom/Atom.h
)
target_link_libraries(main PRIVATE Eigen3::Eigen)
target_compile_features(main PRIVATE cxx_std_20)
