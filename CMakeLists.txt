cmake_minimum_required(VERSION 3.15...3.22)
# check if called from python skbuild
if (NOT DEFINED SKBUILD)
  message(WARNING
          "You should call this through Python so that all variables are set; pip install -v .")
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

# building with scikit-build-core (pip install)
if (DEFINED SKBUILD)
  # wheel.install-dir="/data" is set in pyproject.toml
  # set paths to /platlib and /scripts here
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()

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

file(GLOB BUCCANEER_SRC "include/buccaneer/*.cpp")
#file(GLOB BUCCANEER_HEADERS "include/buccaneer/*.h")
# remove cbuccaneer.cpp and csequins.cpp, not compiling these
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/csequins.cpp")
file(GLOB PY_SRC "python/*.cpp")
# pybind11 module
pybind11_add_module(buc_py MODULE ${PY_SRC} ${BUCCANEER_SRC})# ${BUCCANEER_HEADERS})
target_include_directories(buc_py PRIVATE "dependencies/include")
target_include_directories(buc_py PRIVATE "include")
#find_library(MMDB2DEP NAMES mmdb2)
#find_library(CCP4DEP NAMES ccp4c)
find_library(CLIPPERCOREDEP NAMES clipper-core
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
find_library(CLIPPERMMDBDEP NAMES clipper-mmdb
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
find_library(CLIPPERMINIMOLDEP NAMES clipper-minimol
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
find_library(CLIPPERCONTRIBDEP NAMES clipper-contrib
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
find_library(CLIPPERCCP4DEP NAMES clipper-ccp4
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#message(STATUS "MMDB2  ${MMDB2DEP}")
find_library(CLIPPERGEMMI NAMES clipper-gemmi
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})

#message(STATUS "CCP4C ${CCP4CDEP}")
message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
#message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
message(STATUS "CLIPPERGEMMI ${CLIPPERGEMMI}")
#message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")

find_package(Threads REQUIRED)

target_link_libraries(buc_py PRIVATE #${MMDB2DEP}
                                     #${CCP4DEP}
                                     ${CLIPPERCOREDEP}
                                     #${CLIPPERMMDBDEP}
                                     ${CLIPPERMINIMOLDEP}
                                     ${CLIPPERCONTRIBDEP}
                                     ${CLIPPERGEMMI}
                                     #${CLIPPERCCP4DEP}
                                     Threads::Threads)
set_property(TARGET buc_py PROPERTY OUTPUT_NAME bobkit)
# install
if (DEFINED PYTHON_INSTALL_DIR)
  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
endif()
file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
#install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")