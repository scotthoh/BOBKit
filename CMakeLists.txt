cmake_minimum_required(VERSION 3.15...3.22)

# check if called from python skbuild
if (NOT DEFINED SKBUILD)
  message(WARNING "\
You should call this through Python as this CMake file is meant to be \
executed using 'scikit-build'. If you are installing through the \
source yourself, run the script 'get_dependencies' included in this \
package before running the command 'pip install .'.")
endif()

if(DEFINED ENV{VIRTUAL_ENV})
  message(STATUS "VIRTUAL ENV PATH = $ENV{VIRTUAL_ENV}")
  set(VENV_INSTALL_DIR $ENV{VIRTUAL_ENV})
endif()

# get version from cbuccaneer
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp"
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
message(STATUS "BOBKit version ${PROJECT_VERSION}")

# options
# Use -fPIC even if statically compiled
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(SHOW_COMPILE_TIME "Show compile time" OFF)
if(SHOW_COMPILE_TIME)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
endif()

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "-DBUCCANEER_PROFILE=1") # for profiling in BUCCANEER
# Common installation directories, CMake helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# find python
# We don't need to set policy CMP0069 for this file, because the minimum
# version is 3.15 and for versions >= 3.9 it's automatically set to NEW.
# But a FindPython script in cmake 3.15 have "cmake_policy(VERSION 3.7)"
# so we need to set CMP0069 because of that.
cmake_policy(SET CMP0069 NEW)
cmake_policy(SET CMP0135 NEW)

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
message(STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
# install fftw3 first
include(ExternalProject)
message(STATUS "Trying to install fftw3 ... ")
#if (VENV_INSTALL_DIR)
#  set(FFTW_INSTALL_PATH ${VENV_INSTALL_DIR})
#else()
#  set(FFTW_INSTALL_PATH "${CMAKE_CURRENT_BINARY_DIR}/fftw/install")
#endif()
if(BUILD_SHARED_LIBS)
  set(lib_prefix ${CMAKE_SHARED_LIBRARY_PREFIX})
  set(lib_suffix ${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(fftwconf_args "--enable-shared")
else()
  set(lib_prefix ${CMAKE_STATIC_LIBRARY_PREFIX})
  set(lib_suffix ${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(fftwconf_args "")
endif()
# install fftw3
set(fftw_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/fftw")
set(fftw_src_path "${fftw_dep_path}/src/fftw")
set(fftw_bin_path "${fftw_dep_path}/src/fftw_ex-build/")
set(lib_fftw_name ${lib_prefix}fftw3f${lib_suffix})


ExternalProject_add(project_fftw
                    PREFIX ${fftw_dep_path}
                    SOURCE_DIR ${fftw_src_path}
                    BINARY_DIR ${fftw_bin_path}
                    INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_fftw_name}
                    URL https://fftw.org/pub/fftw/fftw-3.3.4.tar.gz
                    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
                    #BUILD_IN_SOURCE false
                    CONFIGURE_COMMAND <SOURCE_DIR>/configure 
                    --prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install --enable-float 
                    --disable-fortran ${fftwconf_args} "CFLAGS=-fPIC"
                    )
if (BUILD_SHARED_LIBS)
add_library(fftw SHARED IMPORTED)
else()
add_library(fftw STATIC IMPORTED)
endif()
set_target_properties(fftw PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
add_dependencies(fftw project_fftw)
set(FFTW3_INCDIR "${CMAKE_CURRENT_BINARY_DIR}/fftw/install/include")
set(FFTW3_LIBDIR "${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib")
include_directories(${FFTW3_INCDIR})
link_directories(${FFTW3_LIBDIR})
#if (EXISTS ${FFTW3_LIBDIR}/${lib_fftw_name})
#message(STATUS "FFTW3 lib = ${FFTW3_LIBDIR}/${lib_fftw_name}")
#else()
#message(FATAL_ERROR "FFTW3 lib not found!")
#endif()

# install gemmi
set(gemmi_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/gemmi")
set(gemmi_src_path "${gemmi_dep_path}/src/gemmi")
set(gemmi_bin_path "${gemmi_dep_path}/src/gemmi_ex-build/")

set(lib_gemmi_name ${lib_prefix}gemmi_cpp${lib_suffix})

ExternalProject_add(project_gemmi
                    PREFIX ${gemmi_dep_path}
                    SOURCE_DIR ${gemmi_src_path}
                    BINARY_DIR ${gemmi_bin_path}
                    GIT_REPOSITORY https://github.com/project-gemmi/gemmi.git
                    GIT_TAG v0.6.4
                    UPDATE_COMMAND ""
                    BUILD_IN_SOURCE false
                    PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/gemmi/CMakeLists.txt" 
                    "${gemmi_src_path}/CMakeLists.txt"
                    BUILD_COMMAND ${CMAKE_COMMAND} --build ${gemmi_bin_path}
                    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${gemmi_bin_path} --prefix=${VENV_INSTALL_DIR}
                    INSTALL_DIR ${VENV_INSTALL_DIR} #${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_gemmi_name}
                    CMAKE_ARGS -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS})

if(BUILD_SHARED_LIBS)
add_library(gemmi_cpp SHARED IMPORTED)
else()
add_library(gemmi_cpp STATIC IMPORTED)
endif()
set_target_properties(gemmi_cpp PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_gemmi_name})
add_dependencies(gemmi_cpp project_gemmi)
set(GEMMI_INCDIR "${gemmi_src_path}/include")
set(GEMMI_LIBDIR "${VENV_INSTALL_DIR}/lib")
include_directories(${GEMMI_INCDIR})
link_directories(${GEMMI_LIBDIR})


# install parts of clipper
set(clipper_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/clipper")
set(clipper_src_path "${clipper_dep_path}/src/clipper")
set(clipper_bin_path "${clipper_dep_path}/src/clipper_ex-build/")

set(lib_clippercore ${lib_prefix}clipper-core${lib_suffix})
set(lib_clippercontrib ${lib_prefix}clipper-contrib${lib_suffix})
set(lib_clipperminimol ${lib_prefix}clipper-minimol${lib_suffix})
set(lib_clippercif ${lib_prefix}clipper-cif${lib_suffix})
set(lib_clippergemmi ${lib_prefix}clipper-gemmi${lib_suffix})
    
ExternalProject_add(project_clipper
                    PREFIX ${clipper_dep_path}
                    SOURCE_DIR ${clipper_src_path}
                    BINARY_DIR ${clipper_bin_path}
                    GIT_REPOSITORY https://github.com/scotthoh/clipper-gemmi.git
                    GIT_TAG v3.0-dev
                    UPDATE_COMMAND ""
                    BUILD_IN_SOURCE false
                    PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/clipper/CMakeLists.txt" 
                    "${clipper_src_path}/CMakeLists.txt"
                    BUILD_COMMAND ${CMAKE_COMMAND} --build ${clipper_bin_path}
                    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${clipper_bin_path} --prefix=${VENV_INSTALL_DIR}
                    INSTALL_DIR ${VENV_INSTALL_DIR} #${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_clippercore}
                    <INSTALL_DIR>/lib/${lib_clippercif} <INSTALL_DIR>/lib/${lib_clipperminimol}
                    <INSTALL_DIR>/lib/${lib_clippercontrib} <INSTALL_DIR>/lib/${lib_clippergemmi}
                    CMAKE_ARGS -DFFTW_LIBPATH=${FFTW3_LIBDIR}/${lib_fftw_name} -DFFTW_INCPATH=${FFTW3_INCPATH} 
                    -DGEMMI_INCPATH=${GEMMI_INCDIR} -DGEMMI_LIBPATH=${GEMMI_LIBDIR}/${lib_gemmi_name} 
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS})
if(BUILD_SHARED_LIBS)
add_library(clipper-core SHARED IMPORTED)
add_library(clipper-contrib SHARED IMPORTED)
add_library(clipper-gemmi SHARED IMPORTED)
add_library(clipper-minimol SHARED IMPORTED)
add_library(clipper-cif SHARED IMPORTED)
else()
add_library(clipper-core STATIC IMPORTED)
add_library(clipper-contrib STATIC IMPORTED)
add_library(clipper-gemmi STATIC IMPORTED)
add_library(clipper-minimol STATIC IMPORTED)
add_library(clipper-cif STATIC IMPORTED)
endif()
add_dependencies(project_clipper project_fftw project_gemmi)
set_target_properties(clipper-core PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercore})
add_dependencies(clipper-core project_clipper project_fftw)
set_target_properties(clipper-contrib PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercontrib})
add_dependencies(clipper-contrib project_clipper project_fftw)
set_target_properties(clipper-gemmi PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippergemmi})
add_dependencies(clipper-gemmi project_clipper project_gemmi)
set_target_properties(clipper-minimol PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clipperminimol})
add_dependencies(clipper-minimol project_clipper project_gemmi)
set_target_properties(clipper-cif PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercif})
add_dependencies(clipper-cif project_clipper project_gemmi)

set(CLIPPER_INCDIR "${clipper_src_path}/clipper")

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
#add_subdirectory(${WRK_DIR}/dependencies/gemmi)
#add_subdirectory(${WRK_DIR}/dependencies/clipper)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/core)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/gemmi)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/cns)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/cif)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/phs)
#add_subdirectory(${WRK_DIR}/dependencies/clipper/contrib)
#add_subdirectory(${WRK_DIR}/dependencies/fftw3)


#include_directories(${CMAKE_CURRENT_SOURCE_DIR}
#                    #${WRK_DIR}/checkouts/gemmi/include
#                    ${WRK_DIR}/checkouts/clipper
#                    
#                    #${WRK_DIR}/checkouts/fftw3
#                    )

#link_directories(${WRK_DIR}/dependencies/clipper
#                 #${WRK_DIR}/dependencies/fftw3
#                 
#                 ${WRK_DIR}/dependencies/include
#                 #${WRK_DIR}/dependencies/gemmi
#                 ${CMAKE_BINARY_DIR})

file(GLOB BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.cpp")
file(GLOB BUCCANEER_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.h")
# remove cbuccaneer.cpp and csequins.cpp, not compiling these
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/csequins.cpp")
if(BUILD_SHARED_LIBS)
add_library(buccaneer_lib SHARED ${BUCCANEER_SRC})
else()
add_library(buccaneer_lib STATIC ${BUCCANEER_SRC})
endif()
#add_dependencies(buccaneer_lib project_fftw)
#set(LINK_LIBS "gemmi_cpp buccaneer_lib fftw3f clipper-cif clipper-core clipper-contrib clipper-minimol clipper-gemmi")
find_package(Threads REQUIRED)
target_link_directories(buccaneer_lib PUBLIC
                      ${GEMMI_LIBDIR}
                      #/clipper
                      #${VENV_INSTALL_DIR}/lib/clipper/cif
                      #${VENV_INSTALL_DIR}/lib/clipper/core
                      #${VENV_INSTALL_DIR}/lib/clipper/contrib
                      #${VENV_INSTALL_DIR}/lib/clipper/minimol
                      #${VENV_INSTALL_DIR}/lib/clipper/gemmi
                      ${FFTW3_INCDIR})
                      #${CMAKE_BINARY_DIR}/dependencies/fftw3
                      #${CMAKE_BINARY_DIR}/dependencies/gemmi)

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

#message(STATUS "buccaneer_lib = $<TARGET_FILE:buccaner_lib>")
file(GLOB PY_SRC "python/*.cpp")
# pybind11 module
pybind11_add_module(buc_py MODULE ${PY_SRC} ${BUCCANEER_SRC})# ${BUCCANEER_HEADERS})
#add_dependencies(buc_py project_fftw)
target_include_directories(buc_py PUBLIC
                           ${CLIPPER_INCDIR}
                           ${GEMMI_INCDIR}
                           #${WRK_DIR}/checkouts/gemmi
                           #${WRK_DIR}/checkouts/clipper
                           #${WRK_DIR}/checkouts/clipper/clipper/cif
                           #${WRK_DIR}/checkouts/clipper/clipper/minimol
                           #${WRK_DIR}/checkouts/clipper/clipper/gemmi
                           #${WRK_DIR}/checkouts/clipper/clipper/contrib
                           #${WRK_DIR}/checkouts/clipper/clipper/core
                           ${FFTW3_INCDIR}
                           ${CMAKE_SOURCE_DIR}/include
                           #${WRK_DIR}/checkouts/fftw3
                           )
                  
target_link_directories(buc_py PUBLIC
                        ${GEMMI_LIBDIR}
                        #${CMAKE_BINARY_DIR}/dependencies/clipper/cif
                        #${CMAKE_BINARY_DIR}/dependencies/clipper/minimol
                        #${CMAKE_BINARY_DIR}/dependencies/clipper/gemmi
                        #${CMAKE_BINARY_DIR}/dependencies/clipper/contrib
                        #${CMAKE_BINARY_DIR}/dependencies/clipper/core
                        ${FFTW3_LIBDIR})
                        #${CMAKE_BINARY_DIR}/dependencies/fftw3
                        #${CMAKE_BINARY_DIR}/dependencies/gemmi)

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
install(TARGETS buccaneer_lib LIBRARY DESTINATION lib )
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")

#this is for convenience to delete fftw
#set(INSTALL_SHAREDIR ${FFTW_INSTALL_PATH}/share/fftw3)
#install(CODE "string(REPLACE \";\" \"\\n\" MY_CMAKE_INSTALL_MANIFEST_CONTENT \"\$\{CMAKE_INSTALL_MANIFEST_FILES\}\")\n\
#            file(WRITE ${CMAKE_BINARY_DIR}/FFTW3_INSTALL_MANIFEST.txt \"\$\{MY_CMAKE_INSTALL_MANIFEST_CONTENT\}\\n${INSTALL_SHAREDIR}/FFTW3_INSTALL_MANIFEST.txt\")")
#install(FILES "${CMAKE_BINARY_DIR}/FFTW3_INSTALL_MANIFEST.txt" DESTINATION ${INSTALL_SHAREDIR})

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