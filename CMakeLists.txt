cmake_minimum_required(VERSION 3.15)
project(buildkit LANGUAGES CXX) #VERSION ${buildkit_version_str})
# check if called from python skbuild
if (NOT DEFINED SKBUILD)
  message(WARNING "You should call this through Python so that all variables are set; pip install -v .")
endif()
include(ExternalProject)
include(FetchContent)
# still need mmdb2 since we need libccp4
# check for mmdb2 lib in /usr/local and then ccp4 lib dir
#find_library(MMDB2DEP NAMES mmdb2 PATH ${CLIB})
# get mmdb2 from ccp4 gitlab if mmdb2 not there
#set(mmdb2_dep_path ${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2)
#set(mmdb2_bin_path ${mmdb2_dep_path}/src/mmdb-build)
#FetchContent_Declare(mmdb2
#    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/mmdb2.git"
#    GIT_TAG main
#    OVERRIDE_FIND_PACKAGE
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND ""
#    )
#FetchContent_GetProperties(mmdb2)
#if(NOT mmdb2_POPULATED)
#  FetchContent_Populate(mmdb2)
#  add_subdirectory(${mmdb2_SOURCE_DIR}
#  ${mmdb2_BINARY_DIR}
#  EXCLUDE_FROM_ALL)
#endif()
#
#FetchContent_Declare(ccp4
#    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/libccp4.git"
#    GIT_TAG main
#    OVERRIDE_FIND_PACKAGE
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND ""
#    )
#add_dependencies(ccp4 mmdb2)
#FetchContent_MakeAvailable(mmdb2 ccp4)
#ExternalProject_Add(mmdb2
#    PREFIX ${mmdb2_dep_path}
#    SOURCE_DIR ${mmdb2_dep_path}/src/mmdb2
#    BINARY_DIR ${mmdb2_bin_path}
#    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/mmdb2.git"
#    GIT_TAG main
#    BUILD_IN_SOURCE false
#    PATCH_COMMAND ${CMAKE_COMMAND} -E copy
#      "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/mmdb2/CMakeLists.txt" ${mmdb2_dep_path}src/mmdb2/CMakeLists.txt
#    UPDATE_COMMAND ""
#    CONFIGURE_COMMAND "" # cmake 
#    BUILD_COMMAND "" # cmake 
#    INSTALL_COMMAND "" # don't need 
#)
#
#
## get version from cbuccaneer
#file(STRINGS "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp"
#    buccaneer_version_line REGEX "CCP4Program prog\\( \"cbuccaneer\",")
#string(FIND ${buccaneer_version_line} "\"," buccaneer_version_pos)
#string(SUBSTRING ${buccaneer_version_line} ${buccaneer_version_pos} 10 buccaneer_version)
#string(REGEX REPLACE "[,\" ]" "" buccaneer_version_str ${buccaneer_version})
#message(STATUS "cbuccaneer version ${buccaneer_version_str}")
#
## get version from python/version.hpp for python module version
#file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/python/version.hpp"
#      buildkit_version_line REGEX "#define BUILDKIT_VERSION \"")
#string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
#      buildkit_version_str ${buildkit_version_line})
#project(buildkit LANGUAGES CXX VERSION ${buildkit_version_str})
#message(STATUS "Buildkit version ${PROJECT_VERSION}")
#
## building with scikit-build-core (pip install)
if (DEFINED SKBUILD)
  # wheel.install-dir="/data" is set in pyproject.toml
  # set paths to /platlib and /scripts here
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()
#
## options
## Use -fPIC even if statically compiled
#set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
option(USE_PYTHON "Build Python bindings" OFF)
#option(BUILD_SHARED_LIBS "Build as shared library" OFF)
#option(SHOW_COMPILE_TIME "Show compile time" OFF)
#if(SHOW_COMPILE_TIME)
#  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
#endif()
#
## C++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
## Common installation directories, CMake helpers etc.
#include(GNUInstallDirs)
#include(CMakePackageConfigHelpers)
#
#include(FetchContent)
#
## find python
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
#

find_library(CLIPPERCOREDEP NAMES clipper-core
            HINTS ${CMAKE_SOURCE_DIR}/lib )
find_library(CLIPPERMMDBDEP NAMES clipper-mmdb
            HINTS ${CMAKE_SOURCE_DIR}/lib )
find_library(CLIPPERMINIMOLDEP NAMES clipper-minimol
            HINTS ${CMAKE_SOURCE_DIR}/lib )
find_library(CLIPPERCONTRIBDEP NAMES clipper-contrib
            HINTS ${CMAKE_SOURCE_DIR}/lib )
find_library(CLIPPERCCP4DEP NAMES clipper-ccp4
            HINTS ${CMAKE_SOURCE_DIR}/lib )
##message(STATUS "MMDB2  ${MMDB2DEP}")
##message(STATUS "CCP4C ${CCP4DEP}")
message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")
#
find_package(Threads REQUIRED)
#
## extra parts
##include(ExternalProject)
find_library(MMDB2LIB NAMES mmdb2
            HINTS ${CMAKE_SOURCE_DIR}/lib)
if(NOT MMDB2LIB)
set(mmdb2_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2")
set(mmdb2_src_path "${mmdb2_dep_path}/src/mmdb2")
set(mmdb2_bin_path "${mmdb2_dep_path}/src/mmdb2_ex-build/")
#
## mmdb2
#FetchContent_Declare(mmdb2
#  GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/mmdb2.git"
#  GIT_TAG main
#  OVERRIDE_FIND_PACKAGE)
#
#FetchContent_Declare(ccp4
#  GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/libccp4.git"
#  GIT_TAG main
#  OVERRIDE_FIND_PACKAGE)
#
#
ExternalProject_Add(mmdb2_ex
    PREFIX ${mmdb2_dep_path}
    SOURCE_DIR ${mmdb2_src_path}
    BINARY_DIR ${mmdb2_bin_path}
    BUILD_BYPRODUCTS ${CMAKE_SOURCE_DIR}/lib/libmmdb2.so
    #DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2-src"
    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/mmdb2.git"
    GIT_TAG main
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE false
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy
      "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/mmdb2/CMakeLists.txt" ${mmdb2_src_path}/CMakeLists.txt
    
    #CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} --build ${mmdb2_bin_path} #--prefix=${CMAKE_SOURCE_DIR}
    #INSTALL_COMMAND ""
    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${mmdb2_bin_path} --prefix=${CMAKE_SOURCE_DIR}
    #CONFIGURE_COMMAND ${mmdb2_dep_path}/src/mmdb2/configure CXXFLAGS="-O2 -std=c++11" CFLAGS="-O2" --prefix=${CMAKE_INSTALL_PREFIX} --exec_prefix=${CMAKE_INSTALL_PREFIX}/bin --enable-silent-rules --enable-shared --disable-static
    #BUILD_COMMAND make -j8
    #INSTALL_COMMAND make install #"${CMAKE_INSTALL_PREFIX}/lib/
)
#ExternalProject_Get_Property(mmdb2 install_dir)
add_library(mmdb2 SHARED IMPORTED)
set_property(TARGET mmdb2 PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/libmmdb2.so)
add_dependencies(mmdb2 mmdb2_ex)
set(MMDB2_LIB mmdb2)
else()
set(MMDB2_LIB ${MMDB2LIB})
endif()

find_library(CCP4LIB NAMES ccp4
            HINTS ${CMAKE_SOURCE_DIR}/lib)

if( NOT CCP4LIB)
## ccp4
set(ccp4c_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/ccp4")
set(ccp4c_src_path "${ccp4c_dep_path}/src/ccp4")
set(ccp4c_bin_path "${ccp4c_dep_path}/src/ccp4_ex-build/")
ExternalProject_Add(ccp4_ex
    PREFIX ${ccp4c_dep_path}
    SOURCE_DIR ${ccp4c_src_path}
    BINARY_DIR ${ccp4_bin_path}
    BUILD_BYPRODUCTS ${CMAKE_SOURCE_DIR}/lib/libccp4.so
    #DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2-src"
    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/libccp4.git"
    GIT_TAG main
    UPDATE_COMMAND ""
    BUILD_IN_SOURCE false
    PATCH_COMMAND ${CMAKE_COMMAND} -E copy
      "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/ccp4/CMakeLists.txt" ${ccp4c_src_path}/CMakeLists.txt
    #CONFIGURE_COMMAND ""
    BUILD_COMMAND ${CMAKE_COMMAND} --build ${ccp4c_bin_path} #--prefix=${CMAKE_SOURCE_DIR}
    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${ccp4c_bin_path} --prefix=${CMAKE_SOURCE_DIR}
    #CONFIGURE_COMMAND ${mmdb2_dep_path}/src/mmdb2/configure CXXFLAGS="-O2 -std=c++11" CFLAGS="-O2" --prefix=${CMAKE_INSTALL_PREFIX} --exec_prefix=${CMAKE_INSTALL_PREFIX}/bin --enable-silent-rules --enable-shared --disable-static
    #BUILD_COMMAND make -j8
    #INSTALL_COMMAND make install #"${CMAKE_INSTALL_PREFIX}/lib/
)
##add library mmdb2 interface
#file(GLOB MMDB2_SRC "${mmdb2_src_path}/mmdb2/*.cpp")
#file(GLOB MMDB2_HEADERS "${mmdb2_src_path}/mmdb2/*.h")
#add_library(mmdb2 PRIVATE ${MMDB2_SRC} ${MMDB2_HEADERS})
#add_library(mmdb2::mmdb2 ALIAS mmdb2)
#set_property(TARGET mmdb2 PROPERTY POSITION_INDEPENDENT_CODE ON)
#target_include_directories(mmdb2 PRIVATE
#        "$<BUILD_INTERFACE:${mmdb2_src_path}/mmdb2>"
#        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
#
##add ccp4 lib
#file(GLOB CCP4C_SRC "${ccp4c_src_path}/*.cpp")
#file(GLOB CCP4C_HEADERS "${ccp4c_src_path}/*.h")
#add_library(ccp4 PRIVATE ${CCP4C_SRC} ${CCP4C_HEADERS})
#add_library(ccp4::ccp4 ALIAS ccp4)
#set_property(TARGET ccp4 PROPERTY POSITION_INDEPENDENT ON)
#target_include_directories(ccp4 PRIVATE
#        "$<BUILD_INTERFACE:${ccp4c_src_path}/ccp4"
#        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}")
#
## clipper
#set(clipper_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/clipper")
#set(clipper_src_path "${ccp4c_dep_path}/src/clipper")
#set(clipper_bin_path "${ccp4c_dep_path}/src/clipper-build/")
#ExternalProject_Add(ccp4
#    PREFIX ${clipper_dep_path}
#    SOURCE_DIR ${clipper_src_path}
#    BINARY_DIR ${clipper_bin_path}
#    
#    #DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2-src"
#    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/clipper.git"
#    GIT_TAG main
#    UPDATE_COMMAND ""
#    BUILD_IN_SOURCE false
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND ""
#    INSTALL_COMMAND ""
#    #CONFIGURE_COMMAND ${mmdb2_dep_path}/src/mmdb2/configure CXXFLAGS="-O2 -std=c++11" CFLAGS="-O2" --prefix=${CMAKE_INSTALL_PREFIX} --exec_prefix=${CMAKE_INSTALL_PREFIX}/bin --enable-silent-rules --enable-shared --disable-static
#    #BUILD_COMMAND make -j8
#    #INSTALL_COMMAND make install #"${CMAKE_INSTALL_PREFIX}/lib/
#)
##add library clipper need to do it for each individually
#file(GLOB CLIPPER_SRC "${clipper_src_path}/clipper/*.cpp")
#file(GLOB CLIPPER_HEADERS "${clipper_src_path}/mmdb2/*.h")
#add_library(clipper PRIVATE ${MMDB2_SRC} ${MMDB2_HEADERS})
#add_library(clipper::clipper ALIAS clipper)
#set_property(TARGET mmdb2 PROPERTY POSITION_INDEPENDENT_CODE ON)
#target_include_directories(mmdb2 PRIVATE
#        "$<BUILD_INTERFACE:${mmdb2_src_path}/mmdb2>"
#        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
#
#
## buccaneer
#set(buc_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/clipperprogs")
#set(buc_src_path "${buc_dep_path}/src/clipper-progs/src/buccaneer")
#set(buc_bin_path "${buc_dep_path}/src/buc-build/")
#ExternalProject_Add(clipperprogs
#    PREFIX ${buc_dep_path}
#    SOURCE_DIR ${buc_src_path}
#    BINARY_DIR ${ccp4_bin_path}
#    
#    #DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_deps/mmdb2-src"
#    GIT_REPOSITORY "https://ccp4forge.rc-harwell.ac.uk/ccp4/clipper-progs.git"
#    GIT_TAG main
#    UPDATE_COMMAND ""
#    BUILD_IN_SOURCE false
#    CONFIGURE_COMMAND ""
#    BUILD_COMMAND ""
#    INSTALL_COMMAND ""
#    #CONFIGURE_COMMAND ${mmdb2_dep_path}/src/mmdb2/configure CXXFLAGS="-O2 -std=c++11" CFLAGS="-O2" --prefix=${CMAKE_INSTALL_PREFIX} --exec_prefix=${CMAKE_INSTALL_PREFIX}/bin --enable-silent-rules --enable-shared --disable-static
#    #BUILD_COMMAND make -j8
#    #INSTALL_COMMAND make install #"${CMAKE_INSTALL_PREFIX}/lib/
#)
##add cbuccaneer
#file(GLOB BUC_SRC "${buc_src_path}/*.cpp")
#file(GLOB BUC_HEADERS "${buc_src_path}/*.h")
## remove *-NN* and csequins.cpp, not compiling these
#list(REMOVE_ITEM BUC_SRC "${buc_src_path}/*-NN*")
#list(REMOVE_ITEM BUC_SRC "${buc_src_path}/k2c*")
#list(REMOVE_ITEM BUC_HEADERS "${buc_src_path}/*-NN")
#list(REMOVE_ITEM BUC_HEADERS "${buc_src_path}/k2c*")
#list(REMOVE_ITEM BUC_SRC "${buc_src_path}/csequins.cpp")
#
#add_executable(cbuccaneer PRIVATE ${BUC_SRC} ${BUC_HEADERS})
#add_executable(cbuccaneer::cbuccaneer ALIAS cbuccaneer)
#set_property(TARGET cbuccaneer PROPERTY POSITION_INDEPENDENT_CODE ON)
#target_include_directories(cbuccaneer PRIVATE mmdb2 ccp4 clipper)
#
#
## test cbuccaneer build
###add_executable(cbuccaneer ${BUCCANEER_SRC} "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
###add_executable(cbuccaneer::prog ALIAS cbuccaneer)
###target_link_libraries(cbuccaneer PRIVATE clipper)
#ExternalProject_Get_Property(ccp4 install_dir)
add_library(ccp4 SHARED IMPORTED)
set_property(TARGET ccp4 PROPERTY IMPORTED_LOCATION ${CMAKE_SOURCE_DIR}/lib/libccp4.so)
add_dependencies(ccp4 ccp4_ex)
set(CCP4_LIB ccp4)
else()
set(CCP4_LIB ${CCP4LIB})
endif()

#### buccanener
#include_directories(${CMAKE_SOURCE_DIR}/include/mmdb2 ${CMAKE_SOURCE_DIR}/include/ccp4)
#include_directories("/home/swh514/Projects/development/buildkit/buildkit/dependencies/clipper/trunk/clipper")

file(GLOB BUCCANEER_SRC "include/buccaneer/*.cpp")
## remove cbuccaneer.cpp and csequins.cpp, not compiling these
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/csequins.cpp")
#file(GLOB CLIP_SRC "/home/swh514/Projects/development/buildkit/buildkit/dependencies/include")
file(GLOB PY_SRC "python/*.cpp")
## pybind11 module
pybind11_add_module(buc_py MODULE ${PY_SRC} ${BUCCANEER_SRC})# ${CLIP_SRC})

#
##find_library(CCP4DEP NAMES CCP4C)
#
#  
#
#
target_include_directories(buc_py PRIVATE "dependencies/include")
target_include_directories(buc_py PRIVATE "include")
target_include_directories(buc_py PRIVATE "/home/swh514/Projects/development/buildkit/buildkit/dependencies/include")
#ExternalProject_Get_Property(ccp4 ${CMAKE_SOURCE_DIR}/lib)
#include_directories(${CMAKE_SOURCE_DIR}/include)
#add_dependencies(buc_py mmdb2 ccp4)
#add_dependencies(MMDB2LIB mmdb2)
#add_dependencies(CCP4LIB ccp4)
#add_dependencies(buc_py MMDB2LIB CCP4LIB)

target_link_libraries(buc_py PRIVATE ${MMDB2_LIB} #${MMDB2DEP}
                                     ${CCP4_LIB} #${CCP4DEP}
                                     ${CLIPPERCOREDEP}
                                     ${CLIPPERMMDBDEP}
                                     ${CLIPPERMINIMOLDEP}
                                     ${CLIPPERCONTRIBDEP}
                                     ${CLIPPERCCP4DEP}
                                     Threads::Threads)

set_property(TARGET buc_py PROPERTY OUTPUT_NAME buildkit)
#set_target_properties(buc_py PROPERTIES
#                      BUILD_WITH_INSTALL_RPATH TRUE
#                      INSTALL_RPATH_USE_LINK_PATH TRUE
#                      INSTALL_RPATH "\$ORIGIN/lib:${INSTALL_RPATH}")
## install
if(USE_PYTHON)

if (DEFINED PYTHON_INSTALL_DIR)
  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
endif()
file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
####install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")# RUNTIME_DEPENDENCIES ${CMAKE_SOURCE_DIR}/lib)
# ${MMDB2_LIB} ${CCP4_LIB} )
#install(RUNTIME_DEPENDCY_SET )
endif()