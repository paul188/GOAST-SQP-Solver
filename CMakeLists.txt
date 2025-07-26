cmake_minimum_required(VERSION 3.5)

# project name
project(GOAST VERSION 0.1)
#set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

################################################################################
# general settings
################################################################################

set(pybind11_DIR "/opt/software/easybuild-INTEL/software/pybind11/2.9.2-GCCcore-11.3.0/share/cmake/pybind11")

# where to find the modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules)

################################################################################
# Library CMake include
################################################################################
include(cmake/GOAST.cmake)

include_directories($ENV{CPATH})
link_directories($ENV{LIBRARY_PATH})

# This pulled in the imported target Eigen3::Eigen
include_directories(
    /opt/software/easybuild-INTEL/software/Eigen/3.4.0-GCCcore-13.3.0/include/eigen3
    /usr/include/eigen3
    /usr/local/include/eigen3
)


target_link_libraries(GOAST INTERFACE /home/s24pjoha_hpc/openmesh/OpenMesh-11.0.0/lib/libOpenMeshCore.a)
target_link_libraries(GOAST INTERFACE /home/s24pjoha_hpc/openmesh/OpenMesh-11.0.0/lib/libOpenMeshTools.a)

#set(CMAKE_PREFIX_PATH "/opt/software/easybuild-INTEL/software/Eigen/3.4.0-GCCcore-11.3.0/share/eigen3/cmake" ${CMAKE_PREFIX_PATH})
#find_package(Eigen3 REQUIRED)

# Then link Eigen to your targets that need it, e.g.
#target_link_libraries(GOAST INTERFACE Eigen3::Eigen)

target_include_directories(GOAST INTERFACE /home/s24pjoha_hpc/openmesh/OpenMesh-11.0.0/include)

target_include_directories(GOAST INTERFACE /home/s24pjoha_hpc/SuiteSparse/include)

#target_link_libraries(GOAST INTERFACE /home/s24pjoha_hpc/SuiteSparse/lib)

#target_include_directories(GOAST INTERFACE /home/s24pjoha_hpc/SuiteSparse_install/include)

#target_include_directories(GOAST INTERFACE /home/s24pjoha_hpc/SuiteSparse_install/include/suitesparse)

#target_link_libraries(GOAST INTERFACE /home/s24pjoha_hpc/SuiteSparse_install/lib64)

# Find Python package
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
target_link_libraries(GOAST INTERFACE Python3::Python)
find_package(pybind11 REQUIRED)
target_link_libraries(GOAST INTERFACE pybind11::module)
# Including Python3 headers and NumPy headers
target_include_directories(GOAST INTERFACE 
    ${Python3_INCLUDE_DIRS} 
    ${Python3_NumPy_INCLUDE_DIRS}
)

set(CMAKE_PREFIX_PATH "/opt/software/easybuild-INTEL/software/OpenBLAS/0.3.20-GCC-11.3.0/lib/cmake/openblas" CACHE PATH "")

# Add the OpenBLAS library directory
link_directories(/opt/software/easybuild-INTEL/software/OpenBLAS/0.3.27-GCC-13.3.0/lib)

# If you want to be precise and portable:
find_library(OPENBLAS_LIB NAMES openblas
PATHS /opt/software/easybuild-INTEL/software/OpenBLAS/0.3.27-GCC-13.3.0/lib)

if (NOT OPENBLAS_LIB)
    message(FATAL_ERROR "OpenBLAS library not found!")
endif()

# Link against it
target_link_libraries(GOAST INTERFACE ${OPENBLAS_LIB})


if(MSVC)
    add_compile_options("/fp:strict" "/arch:AVX2" "/openmp")
endif()

###########################################
# Use MKL (via Eigen) in the examples etc
###########################################
option(GOAST_WITH_MKL "Use Intel's Math Kernel Library" OFF)
if (GOAST_WITH_MKL)
  # target_link_options(GOAST INTERFACE -m64 -Wl,--no-as-needed)
  # interface link flags not supported in CMake 3.5, hence
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed")

  target_compile_options(GOAST INTERFACE -m64)
  target_include_directories(GOAST INTERFACE $ENV{MKLROOT}/include)

  target_link_directories(GOAST INTERFACE $ENV{MKLROOT}/lib/intel64)
  target_link_libraries(GOAST INTERFACE mkl_intel_lp64 mkl_gnu_thread mkl_core gomp pthread m dl)

  target_compile_definitions(GOAST INTERFACE GOAST_WITH_MKL)
endif (GOAST_WITH_MKL)

################################################################################
# source directories
################################################################################
option(GOAST_BUILD_TESTS "Build GOAST tests" ON)
option(GOAST_BUILD_EXAMPLES "Build examples of GOAST" ON)
option(GOAST_BUILD_TOOLS "Build tools included with GOAST" ON)

if (GOAST_BUILD_TESTS)
  add_subdirectory(tests)
endif ()

if (GOAST_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif ()

if (GOAST_BUILD_TOOLS)
  add_subdirectory(tools)
endif ()


################################################################################
# install directories
################################################################################
include(GNUInstallDirs)


install(TARGETS GOAST GOAST_all
        EXPORT GOASTTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        )

install(EXPORT GOASTTargets
        FILE GOASTTargets.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GOAST
        )

include(CMakePackageConfigHelpers)
write_basic_package_version_file("GOASTConfigVersion.cmake"
                                 VERSION ${GOAST_VERSION}
                                 COMPATIBILITY SameMajorVersion)

configure_file(cmake/GOASTConfig.cmake.in GOASTConfig.cmake @ONLY)

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/GOASTConfig.cmake"
        "${CMAKE_CURRENT_BINARY_DIR}/GOASTConfigVersion.cmake"
        "cmake/Modules/FindSuiteSparse.cmake"
        "cmake/Modules/FindADOLC.cmake"
        "cmake/Modules/FindEigen3.cmake"
        "cmake/Modules/FindIpopt.cmake"
        "cmake/Modules/FindMOSEK.cmake"
        "cmake/Modules/FindOpenMesh.cmake"
        "cmake/Modules/Findtrlib.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/GOAST)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/goast DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
