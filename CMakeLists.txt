cmake_minimum_required(VERSION 3.15...3.22)

# check if called from python skbuild
if (NOT DEFINED SKBUILD)
  message(WARNING
          "You should call this through Python as this CMake file is meant to be
          executed using 'scikit-build'. If you are installing through the source
          yourself, run the script 'get_dependencies' included in this package
          before running the command 'pip install .'.")
endif()


# get version from cbuccaneer
file(STRINGS "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp"
    buccaneer_version_line REGEX "CCP4Program prog\\( \"cbuccaneer\",")
string(FIND ${buccaneer_version_line} "\"," buccaneer_version_pos)
string(SUBSTRING ${buccaneer_version_line} ${buccaneer_version_pos} 10 buccaneer_version)
string(REGEX REPLACE "[,\" ]" "" buccaneer_version_str ${buccaneer_version})
message(STATUS "cbuccaneer version ${buccaneer_version_str}")

# get version from python/version.hpp for python module version
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/python/version.hpp"
      bobkit_version_line REGEX "#define BOBKIT_VERSION \"")
string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
      bobkit_version_str ${bobkit_version_line})
project(BOBKIT LANGUAGES CXX VERSION ${bobkit_version_str})
message(STATUS "BOBKIT version ${PROJECT_VERSION}")
# options
# Use -fPIC even if statically compiled
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
option(BUILD_SHARED_LIBS "Build as shared library" OFF)
option(SHOW_COMPILE_TIME "Show compile time" OFF)
if(SHOW_COMPILE_TIME)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
endif()

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-DBUCCANEER_PROFILE=1")
# Common installation directories, CMake helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# find python
find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
if (EXISTS "${CMAKE_HOME_DIRECTORY}/pybind11")
  message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/pybind11 (internal copy).")
  add_subdirectory(pybind11)
else()
  find_package(pybind11 2.6 CONFIG)
  if(NOT pybind11_FOUND)
    #try once more, using pybind11-config
    execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
    string(STRIP "${pybind11_DIR}" pybind11_DIR)
    find_package(pybind11 2.6 CONFIG REQUIRED)
  endif()
  message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
endif()
# building with scikit-build-core (pip install)
if (DEFINED SKBUILD)
  # wheel.install-dir="/data" is set in pyproject.toml
  # set paths to /platlib and /scripts here
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

set(WRK_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
message(CMAKE_BINARY_DIR="${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
# install fftw3 first
include(ExternalProject)
message(STATUS "Trying to install fftw3 ... ")

set(lib_fftw_name
    ${CMAKE_STATIC_LIBRARY_PREFIX}fftw3f${CMAKE_STATIC_LIBRARY_SUFFIX})
set(fftw_src_dir ${WRK_DIR}/checkouts/fftw-3.3.4)
ExternalProject_add(project_fftw
                    SOURCE_DIR ${fftw_src_dir}
                    #URL https://fftw.org/pub/fftw/fftw-3.3.4.tar.gz
                    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
                    CONFIGURE_COMMAND ${fftw_src_dir}/configure
                    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install --enable-float 
                    "CFLAGS=-fPIC"
                    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_fftw_name})

add_library(fftw STATIC IMPORTED)
#set(lib_fftw_name
#    ${CMAKE_SHARED_LIBRARY_PREFIX}fftw3f${CMAKE_SHARED_LIBRARY_SUFFIX})
set_target_properties(fftw PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
add_dependencies(fftw project_fftw)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include)

set(FFTW3_INCDIR "${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include")
set(FFTW3_LIBDIR "${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib")
include_directories(${FFTW3_INCDIR})
link_directories(${FFTW3_LIBDIR})
message(STATUS "FFTW3 lib = ${FFTW3_LIBDIR}")
#find_library(FFTW3_DEP NAMES fftw3f
#            HINTS ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib)
#
#message(STATUC "FFTW3 found = ${FFTW3_DEP}")

#set(ccp4c_src_path "${CMAKE_SOURCE_DIR}/checkouts/libccp4")
#
#ExternalProject_Add(ccp4_ex
#    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/libccp4
#    SOURCE_DIR ${ccp4c_src_path}
#    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/libccp4/install
#    BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/libccp4/bin
#    BUILD_BYPRODUCTS ${CMAKE_SOURCE_DIR}/lib/libccp4.so
#    #DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2-src"
#    UPDATE_COMMAND ""
#    BUILD_IN_SOURCE false
#    PATCH_COMMAND ${CMAKE_COMMAND} -E copy
#      "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/libccp4/CMakeLists.txt" ${ccp4c_src_path}/CMakeLists.txt
#    #CONFIGURE_COMMAND ""
#    BUILD_COMMAND ${CMAKE_COMMAND} --build $<BINARY_DIR> #--prefix=${CMAKE_SOURCE_DIR}
#    INSTALL_COMMAND ${CMAKE_COMMAND} --install $<BINARY_DIR> --prefix=${CMAKE_CURRENT_BINARY_DIR}/libccp4
#)


# set subdirectories for dependencies
add_subdirectory(${WRK_DIR}/dependencies/gemmi)
add_subdirectory(${WRK_DIR}/dependencies/clipper)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/core)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/gemmi)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/cns)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/cif)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/phs)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/contrib)
#add_subdirectory(${WRK_DIR}/dependencies/fftw3)


include_directories(${CMAKE_CURRENT_SOURCE_DIR}
                    ${WRK_DIR}/checkouts/gemmi/include
                    ${WRK_DIR}/checkouts/clipper
                    
                    #${WRK_DIR}/checkouts/fftw3
                    )

link_directories(${WRK_DIR}/dependencies/clipper
                 #${WRK_DIR}/dependencies/fftw3
                 
                 ${WRK_DIR}/dependencies/include
                 ${WRK_DIR}/dependencies/gemmi
                 ${CMAKE_BINARY_DIR})

file(GLOB BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.cpp")
file(GLOB BUCCANEER_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.h")
# remove cbuccaneer.cpp and csequins.cpp, not compiling these
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/csequins.cpp")
add_library(buccaneer_lib STATIC ${BUCCANEER_SRC})
add_dependencies(buccaneer_lib project_fftw)
#set(LINK_LIBS "gemmi_cpp buccaneer_lib fftw3f clipper-cif clipper-core clipper-contrib clipper-minimol clipper-gemmi")
message(STATUS "clipper-core = " $<TARGET_FILE:clipper-core>)
message(STATUS "clipper-contrib = " $<TARGET_FILE:clipper-contrib>)
message(STATUS "clipper-gemmi = " $<TARGET_FILE:clipper-gemmi>)
message(STATUS "clipper-minimol =  " $<TARGET_FILE:clipper-minimol>)
message(STATUS "clipper-cif =  " $<TARGET_FILE:clipper-cif>)
message(STATUS "gemmi_cpp =  " $<TARGET_FILE:gemmi_cpp>)

find_package(Threads REQUIRED)
target_link_directories(buccaneer_lib PUBLIC
                      ${CMAKE_BINARY_DIR}/dependencies/clipper
                      ${CMAKE_BINARY_DIR}/dependencies/clipper/cif
                      ${CMAKE_BINARY_DIR}/dependencies/clipper/core
                      ${CMAKE_BINARY_DIR}/dependencies/clipper/contrib
                      ${CMAKE_BINARY_DIR}/dependencies/clipper/minimol
                      ${CMAKE_BINARY_DIR}/dependencies/clipper/gemmi
                      ${FFTW3_INCDIR}
                      #${CMAKE_BINARY_DIR}/dependencies/fftw3
                      ${CMAKE_BINARY_DIR}/dependencies/gemmi)

target_link_libraries(buccaneer_lib PUBLIC
                      clipper-cif
                      clipper-minimol
                      clipper-contrib
                      clipper-gemmi
                      clipper-core
                      fftw
                      #fftw
                      gemmi_cpp
                      Threads::Threads)

message(STATUS "buccaneer_lib = $<TARGET_FILE:buccaner_lib>")
file(GLOB PY_SRC "python/*.cpp")
# pybind11 module
pybind11_add_module(buc_py MODULE ${PY_SRC} ${BUCCANEER_SRC})# ${BUCCANEER_HEADERS})
#add_dependencies(buc_py project_fftw)
target_include_directories(buc_py PUBLIC
                           ${WRK_DIR}/checkouts/gemmi
                           ${WRK_DIR}/checkouts/clipper
                           ${WRK_DIR}/checkouts/clipper/clipper/cif
                           ${WRK_DIR}/checkouts/clipper/clipper/minimol
                           ${WRK_DIR}/checkouts/clipper/clipper/gemmi
                           ${WRK_DIR}/checkouts/clipper/clipper/contrib
                           ${WRK_DIR}/checkouts/clipper/clipper/core
                           #${FFTW3_INCDIR}
                           ${CMAKE_SOURCE_DIR}/include
                           #${WRK_DIR}/checkouts/fftw3
                           )
                  
target_link_directories(buc_py PUBLIC
                        ${CMAKE_BINARY_DIR}/dependencies/clipper/cif
                        ${CMAKE_BINARY_DIR}/dependencies/clipper/minimol
                        ${CMAKE_BINARY_DIR}/dependencies/clipper/gemmi
                        ${CMAKE_BINARY_DIR}/dependencies/clipper/contrib
                        ${CMAKE_BINARY_DIR}/dependencies/clipper/core
                        #${FFTW3_LIBDIR}
                        #${CMAKE_BINARY_DIR}/dependencies/fftw3
                        ${CMAKE_BINARY_DIR}/dependencies/gemmi)

target_link_libraries(buc_py PRIVATE buccaneer_lib
                      clipper-cif
                      clipper-contrib
                      clipper-minimol
                      clipper-gemmi
                      clipper-core
                      #fftw
                      #${FFTW3_LIBDIR}
                      gemmi_cpp
                      Threads::Threads)

#install
#dinstall(TARGETS buc_py LIBRARY DESTINATION bobkit)
# install
set_property(TARGET buc_py PROPERTY OUTPUT_NAME bobkit)
if (DEFINED PYTHON_INSTALL_DIR)
  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
endif()
file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
#install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")

#target_include_directories(buc_py PRIVATE "dependencies/include")
#target_include_directories(buc_py PRIVATE "include")
##find_library(MMDB2DEP NAMES mmdb2)
##find_library(CCP4DEP NAMES ccp4c)
#find_library(CLIPPERCOREDEP NAMES clipper-core
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#find_library(CLIPPERMMDBDEP NAMES clipper-mmdb
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#find_library(CLIPPERMINIMOLDEP NAMES clipper-minimol
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#find_library(CLIPPERCONTRIBDEP NAMES clipper-contrib
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#find_library(CLIPPERCCP4DEP NAMES clipper-ccp4
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
##message(STATUS "MMDB2  ${MMDB2DEP}")
#find_library(CLIPPERGEMMI NAMES clipper-gemmi
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#
##message(STATUS "CCP4C ${CCP4CDEP}")
#message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
##message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
#message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
#message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
#message(STATUS "CLIPPERGEMMI ${CLIPPERGEMMI}")
##message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")
#
#find_package(Threads REQUIRED)
#
#target_link_libraries(buc_py PRIVATE #${MMDB2DEP}
#                                     #${CCP4DEP}
#                                     ${CLIPPERCOREDEP}
#                                     #${CLIPPERMMDBDEP}
#                                     ${CLIPPERMINIMOLDEP}
#                                     ${CLIPPERCONTRIBDEP}
#                                     ${CLIPPERGEMMI}
#                                     #${CLIPPERCCP4DEP}
#                                     Threads::Threads)
#set_property(TARGET buc_py PROPERTY OUTPUT_NAME bobkit)
## install
#if (DEFINED PYTHON_INSTALL_DIR)
#  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
#  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
#endif()
#file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
##install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
#install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")