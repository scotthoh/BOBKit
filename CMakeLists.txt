cmake_minimum_required(VERSION 3.20)

# check if called from python skbuild
if (NOT DEFINED SKBUILD)
  message(WARNING "\
You should call this through Python as this CMake file is meant to be \
executed using 'scikit-build'. If you are installing through the \
source yourself, run the script 'get_dependencies' included in this \
package before running the command 'pip install .'.")
endif()

#if(NOT APPLE)
# set(CMAKE_INSTALL_RPATH $ORIGIN)
#endif()
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
option(USE_PYTHON "Build Python bindings" ON)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
option(BUILD_OWN_FFTW3 "Build own version of fftw3 ignoring the one found on system" ON)
option(BUCCANEER_PROFILE "Set Buccaneer profiling on/off" ON)
#option(GENERATE_STUBS "Generate Python type stubs" OFF)
if(SHOW_COMPILE_TIME)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
endif()

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
if(BUCCANEER_PROFILE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBUCCANEER_PROFILE=1") # for profiling in BUCCANEER
endif()
# Common installation directories, CMake helpers
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
if(DEFINED ENV{VIRTUAL_ENV})
  message(STATUS "VIRTUAL ENV PATH = $ENV{VIRTUAL_ENV}")
  set(VENV_INSTALL_DIR $ENV{VIRTUAL_ENV})
else()
  set(VENV_INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif()

if (NOT BUILD_OWN_FFTW3)
  if(DEFINED VENV_INSTALL_DIR)
    set(TMP_CMAKE_LIB_PATH ${CMAKE_LIBRARY_PATH})
    set(CMAKE_LIBRARY_PATH ${VENV_INSTALL_DIR}/lib)
    message(STATUS ${CMAKE_LIBRARY_PATH})
    message(STATUS ${CMAKE_INSTALL_PREFIX})
    message(STATUS ${CMAKE_INCLUDE_PATH})
  endif()
  message(STATUS "Finding FFTW3 package.")
  #set(BUILD_OWN_FFTW3 OFF CACHE BOOL "" FORCE)
  find_package(PkgConfig REQUIRED)
  pkg_search_module(FFTW REQUIRED IMPORTED_TARGET GLOBAL fftw3f>=3.3.4)
  if(FFTW_FOUND)
    message(STATUS "Found FFTW3 - ${FFTW_LINK_LIBRARIES}")
    message(STATUS "FFTW3 include dir - ${FFTW_INCLUDE_DIRS}")
    set(CMAKE_LIBRARY_PATH ${TMP_CMAKE_LIB_PATH})
  elseif()
    set(BUILD_OWN_FFTW3 ON CACHE BOOL "" FORCE)
  endif()
endif()

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
if(POLICY CMP0069)
  cmake_policy(SET CMP0069 NEW)
endif()
if(POLICY CMP0135)
  cmake_policy(SET CMP0135 NEW)
endif()
if (NOT DEFINED CMAKE_INTERPROCEDURAL_OPTIMIZATION)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT is_ipo_supported OUTPUT ipo_error)
  if(is_ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_MINSIZEREL ON)
  else()
    message(STATUS "Interprocedural optimization not supported: ${error}")
  endif()
endif()

if (BUILD_SHARED_LIBS)
  if (APPLE)
    set(origin "@loader_path")
  else()
    set(origin "$ORIGIN")
  endif()
endif()

if(USE_PYTHON)
find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
if (EXISTS "${CMAKE_HOME_DIRECTORY}/nanobind") #pybind11")
  message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/nanobind (internal copy).")
  add_subdirectory(nanobind)
else()
  # Detect the installed nanobind package and import it into CMake
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nb_cmake_dir)
  list(APPEND CMAKE_PREFIX_PATH "${nb_cmake_dir}")
  find_package(nanobind 2.4.0 CONFIG REQUIRED)
  message(STATUS "Found nanobind ${nanobind_VERSION}: ${NB_DIR}")
  ## use pybind11-config (if available) to determine pybind11_DIR
  #execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
  #string(STRIP "${pybind11_DIR}" pybind11_DIR)
  #find_package(pybind11 2.6 CONFIG)
  ##if(NOT pybind11_FOUND)
  ##  #try once more, using pybind11-config
  ##  execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
  ##  string(STRIP "${pybind11_DIR}" pybind11_DIR)
  ##  find_package(pybind11 2.6 CONFIG REQUIRED)
  ##endif()
  #message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
endif()
# building with scikit-build-core (pip install)

if (DEFINED SKBUILD)
  # wheel.install-dir="/data" is set in pyproject.toml
  # set paths to /platlib and /scripts here
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
  if (DEFINED PYTHON_INSTALL_DIR)
    message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
    set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
  endif()
  file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
endif()
endif()

set(WRK_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "CMAKE_BINARY_DIR = ${CMAKE_BINARY_DIR}")
# install fftw3 first
include(ExternalProject)
if(BUILD_SHARED_LIBS)
  set(lib_prefix ${CMAKE_SHARED_LIBRARY_PREFIX})
  set(lib_suffix ${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(fftwconf_args "--enable-shared")
  if (USE_PYTHON)
  if (DEFINED SKBUILD)
    set(rpath_rel "../../../../${CMAKE_INSTALL_LIBDIR}")
  else()
    file(RELATIVE_PATH rpath_rel "${Python_SITEARCH}" "${VENV_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}")
  endif()
  endif()
else()
  message(STATUS "static lib prefix" ${CMAKE_STATIC_LIBRARY_PREFIX})
  set(lib_prefix ${CMAKE_STATIC_LIBRARY_PREFIX})
  set(lib_suffix ${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(fftwconf_args "")
endif()
set(lib_fftw_name ${lib_prefix}fftw3f${lib_suffix})
if(BUILD_OWN_FFTW3)
message(STATUS "Setting up to install FFTW3...")
# install fftw3
set(fftw_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/fftw")
set(fftw_src_path "${fftw_dep_path}/src/fftw")
set(fftw_bin_path "${fftw_dep_path}/src/fftw_ex-build")

ExternalProject_add(project_fftw
                    PREFIX ${fftw_dep_path}
                    SOURCE_DIR ${fftw_src_path}
                    BINARY_DIR ${fftw_bin_path}
                    INSTALL_DIR ${VENV_INSTALL_DIR}
                    #INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_fftw_name}
                    URL https://fftw.org/pub/fftw/fftw-3.3.4.tar.gz
                    #PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw
                    BUILD_IN_SOURCE false
                    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=${VENV_INSTALL_DIR} --enable-float 
                    --disable-fortran ${fftwconf_args} "CFLAGS=-fPIC"
                    #--prefix=${CMAKE_CURRENT_BINARY_DIR}/fftw/install --enable-float 
                    )
if (BUILD_SHARED_LIBS)
add_library(fftw3::fftw3 SHARED IMPORTED)
else()
add_library(fftw3::fftw3 STATIC IMPORTED)
endif()
#set_target_properties(fftw PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/fftw/install/lib/${lib_fftw_name})
set(FFTW3_INCDIR "${VENV_INSTALL_DIR}/include")
set(FFTW3_LIBDIR "${VENV_INSTALL_DIR}/lib")
message(STATUS "FFTW3 lib : " ${FFTW3_LIBDIR}/${lib_fftw_name})
set_target_properties(fftw3::fftw3 PROPERTIES
                      IMPORTED_LOCATION ${FFTW3_LIBDIR}/${lib_fftw_name}
                      #INTERFACE_INCLUDE_DIRECTORIES ${FFTW3_INCDIR}
)
add_dependencies(fftw3::fftw3 project_fftw)
#set(fftw3_extlib fftw3::fftw3)
else()
set(FFTW3_INCDIR "${FFTW_INCLUDE_DIRS}")
set(FFTW3_LIBDIR "${FFTW_LIBRARY_DIRS}")
#cmake_path(GET FFTW3_LIBRARY PARENT_PATH FFTW3_LIBDIR)
add_library(fftw3::fftw3 ALIAS PkgConfig::FFTW)
#if (BUILD_SHARED_LIBS)
#add_library(project_fftw SHARED IMPORTED)
#else()
#add_library(project_fftw STATIC IMPORTED)
#endif()
#add_library(fftw3::fftw3 ALIAS project_fftw)
##add_dependencies(fftw3::fftw3 project_fftw)
#set_target_properties(project_fftw PROPERTIES IMPORTED_LOCATION ${FFTW_LINK_LIBRARIES})
#                      #INTERFACE_INCLUDE_DIRECTORIES ${FFTW_INCLUDE_DIRS})
endif() # if(BUILD_OWN_FFTW3)

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
                    GIT_TAG v0.7.3
                    #GIT_TAG v0.6.7
                    UPDATE_COMMAND ""
                    BUILD_IN_SOURCE false
                    #PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/gemmi/CMakeLists.txt" 
                    #"${gemmi_src_path}/CMakeLists.txt"
                    BUILD_COMMAND ${CMAKE_COMMAND} --build ${gemmi_bin_path}
                    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${gemmi_bin_path} #--prefix=${VENV_INSTALL_DIR}
                    INSTALL_DIR ${VENV_INSTALL_DIR} #${CMAKE_CURRENT_BINARY_DIR}/fftw/install
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_gemmi_name}
                    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${VENV_INSTALL_DIR} -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DBUILD_GEMMI_PROGRAM=OFF -DUSE_PYTHON=ON -DSTANDALONE_PYTHON_MODULE=OFF)

if(BUILD_SHARED_LIBS)
add_library(gemmi::gemmi_cpp SHARED IMPORTED)
else()
add_library(gemmi::gemmi_cpp STATIC IMPORTED)
endif()
set(GEMMI_INCDIR "${VENV_INSTALL_DIR}/include")
set(GEMMI_LIBDIR "${VENV_INSTALL_DIR}/lib")
set_target_properties(gemmi::gemmi_cpp PROPERTIES
                      IMPORTED_LOCATION ${GEMMI_LIBDIR}/${lib_gemmi_name})
                      #INTERFACE_INCLUDE_DIRECTORIES ${GEMMI_INCDIR})
add_dependencies(gemmi::gemmi_cpp project_gemmi)

# install parts of clipper
set(clipper_dep_path "${CMAKE_CURRENT_BINARY_DIR}/_deps/clipper")
set(clipper_src_path "${clipper_dep_path}/src/clipper")
#set(clipper_src_path "/Users/swh514/Work/development/clipper/clipper-gemmi")
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
                    #GIT_REPOSITORY https://github.com/scotthoh/clipper-gemmi.git
                    GIT_REPOSITORY https://ccp4forge.rc-harwell.ac.uk/ccp4/clipper.git
                    GIT_TAG b1fce0e931043aaf6c466cae1354842d511a230d
                    #GIT_TAG origin/large_maps_support
                    UPDATE_COMMAND ""
                    BUILD_IN_SOURCE false
                    PATCH_COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/clipper/CMakeLists.txt" 
                    "${clipper_src_path}/CMakeLists.txt"
                    BUILD_COMMAND ${CMAKE_COMMAND} --build ${clipper_bin_path}
                    INSTALL_COMMAND ${CMAKE_COMMAND} --install ${clipper_bin_path} --prefix=${VENV_INSTALL_DIR}
                    INSTALL_DIR ${VENV_INSTALL_DIR}
                    BUILD_BYPRODUCTS <INSTALL_DIR>/lib/${lib_clippercore}
                    <INSTALL_DIR>/lib/${lib_clippercif} <INSTALL_DIR>/lib/${lib_clipperminimol}
                    <INSTALL_DIR>/lib/${lib_clippercontrib} <INSTALL_DIR>/lib/${lib_clippergemmi}
                    CMAKE_ARGS -DFFTW_LIBPATH=${FFTW3_LIBDIR}/${lib_fftw_name} -DFFTW_INCPATH=${FFTW3_INCDIR} 
                    -DGEMMI_INCPATH=${GEMMI_INCDIR} -DGEMMI_LIBPATH=${GEMMI_LIBDIR}/${lib_gemmi_name} 
                    -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DREL_RPATH=${rpath_rel})
if(BUILD_SHARED_LIBS)
add_library(clipper::core SHARED IMPORTED)
add_library(clipper::contrib SHARED IMPORTED)
add_library(clipper::gemmi SHARED IMPORTED)
add_library(clipper::minimol SHARED IMPORTED)
add_library(clipper::cif SHARED IMPORTED)
else()
add_library(clipper::core STATIC IMPORTED)
add_library(clipper::contrib STATIC IMPORTED)
add_library(clipper::gemmi STATIC IMPORTED)
add_library(clipper::minimol STATIC IMPORTED)
add_library(clipper::cif STATIC IMPORTED)
endif()
if (BUILD_OWN_FFTW3)
  add_dependencies(project_clipper project_gemmi project_fftw )
else()
  add_dependencies(project_clipper project_gemmi)
endif()
set(CLIPPER_INCDIR "${VENV_INSTALL_DIR}/include")
set_target_properties(clipper::core PROPERTIES
                      IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercore})
                      #INTERFACE_INCLUDE_DIRECTORIES ${CLIPPER_INCDIR})
add_dependencies(clipper::core project_clipper) # project_fftw)
set_target_properties(clipper::contrib PROPERTIES
                      IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercontrib})
                      #INTERFACE_INCLUDE_DIRECTORIES ${CLIPPER_INCDIR})
add_dependencies(clipper::contrib project_clipper) # project_fftw)
set_target_properties(clipper::gemmi
                      PROPERTIES IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippergemmi})
                      #INTERFACE_INCLUDE_DIRECTORIES ${CLIPPER_INCDIR})
add_dependencies(clipper::gemmi project_clipper project_gemmi)
set_target_properties(clipper::minimol PROPERTIES
                      IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clipperminimol})
                      #INTERFACE_INCLUDE_DIRECTORIES ${CLIPPER_INCDIR})
add_dependencies(clipper::minimol project_clipper project_gemmi)
set_target_properties(clipper::cif PROPERTIES
                      IMPORTED_LOCATION ${VENV_INSTALL_DIR}/lib/${lib_clippercif})
                      #'INTERFACE_INCLUDE_DIRECTORIES ${CLIPPER_INCDIR})
add_dependencies(clipper::cif project_clipper project_gemmi)

file(GLOB BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.cpp")
file(GLOB BUCCANEER_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/*.h")
# remove cbuccaneer.cpp and csequins.cpp, not compiling these
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/include/buccaneer/csequins.cpp")
if(BUILD_SHARED_LIBS)
  add_library(buccaneer_lib SHARED ${BUCCANEER_SRC} ${BUCCANEER_HEADERS})
  set_target_properties(buccaneer_lib PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
else()
add_library(buccaneer_lib STATIC ${BUCCANEER_SRC} ${BUCCANEER_HEADERS})
endif()
add_library(buccaneer::buccaneer_lib ALIAS buccaneer_lib)
set_target_properties(buccaneer_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)
add_dependencies(buccaneer_lib project_gemmi project_clipper)
find_package(Threads REQUIRED)
target_include_directories(buccaneer_lib PUBLIC ${VENV_INSTALL_DIR}/include ${CMAKE_SOURCE_DIR}/include)

target_link_directories(buccaneer_lib PRIVATE
                      ${VENV_INSTALL_DIR}/lib
                      ${clipper_bin_path}
                      ${gemmi_bin_path}
                      ${FFTW3_LIBDIR})
                      
target_link_libraries(buccaneer_lib PRIVATE
                      clipper::cif
                      clipper::minimol
                      clipper::contrib
                      clipper::gemmi
                      clipper::core
                      fftw3::fftw3
                      gemmi::gemmi_cpp)
if (Threads_FOUND)
  target_link_libraries(buccaneer_lib PRIVATE Threads::Threads)
endif()
message(STATUS "VENV_INSTALL_DIR is ${VENV_INSTALL_DIR}")
install(TARGETS buccaneer_lib LIBRARY DESTINATION "${VENV_INSTALL_DIR}/lib")
file(GLOB PY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/python/*.cpp")
list(REMOVE_ITEM PY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/python/c_map_io.cpp")
list(REMOVE_ITEM PY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/python/c_container.cpp")
#list(REMOVE_ITEM PY_SRC "${CMAKE_CURRENT_SOURCE_DIR}/python/c_gemmi.cpp")


# nanobind module
if(USE_PYTHON)
nanobind_add_module(buc_py ${PY_SRC})# $<TARGET_OBJECTS:buccaneer_lib>) #${BUCCANEER_SRC})# ${BUCCANEER_HEADERS})
add_dependencies(buc_py project_gemmi project_clipper buccaneer_lib)

target_include_directories(buc_py PUBLIC
                           ${CLIPPER_INCDIR}
                           ${GEMMI_INCDIR}
                           ${FFTW3_INCDIR}
                           ${CMAKE_SOURCE_DIR}/include
                          )

#target_link_directories(buc_py PRIVATE
#                        ${VENV_INSTALL_DIR}/lib
#                        ${FFTW3_LIBDIR}
#                        ${clipper_bin_path}
#                        ${gemmi_bin_path})
#
target_link_libraries(buc_py PRIVATE buccaneer::buccaneer_lib 
clipper::cif
clipper::contrib
clipper::minimol
clipper::gemmi
clipper::core
fftw3::fftw3
gemmi::gemmi_cpp)
#endif() # STANDALONE_PYTHON_MODULE

if (Threads_FOUND)
  target_link_libraries(buc_py PRIVATE Threads::Threads)
endif()
# install
set_property(TARGET buc_py PROPERTY OUTPUT_NAME bobkit_ext)
if (CMAKE_CONFIGURATION_TYPES)
  set(py_dir $<CONFIG>/py)
else()
  set(py_dir py)
endif()
set_property(TARGET buc_py PROPERTY LIBRARY_OUTPUT_DIRECTORY ${py_dir}/bobkit)
add_custom_command(TARGET buc_py PRE_BUILD
  COMMENT "Copying bobkit/__init__.py"
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/python/bobkit/__init__.py" ${py_dir}/bobkit/__init__.py
  BYPRODUCTS ${py_dir}/bobkit/__init__.py
)
add_custom_command(TARGET buc_py PRE_BUILD
  COMMENT "Copying bobkit/buccaneer.py"
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/python/bobkit/buccaneer.py" ${py_dir}/bobkit/buccaneer.py
  BYPRODUCTS ${py_dir}/bobkit/buccaneer.py
)
add_custom_command(TARGET buc_py PRE_BUILD
  COMMENT "Copying bobkit/util.py"
  COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_CURRENT_SOURCE_DIR}/python/bobkit/util.py" ${py_dir}/bobkit/util.py
  BYPRODUCTS ${py_dir}/bobkit/util.py
)
add_custom_command(TARGET buc_py PRE_BUILD
  COMMENT "Copying bobkit/_scripts"
  COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/python/bobkit/_scripts" ${py_dir}/bobkit/_scripts
  BYPRODUCTS ${py_dir}/bobkit/_scripts
)


#execute_process(
#  COMMAND ${Python_EXECUTABLE} -c "import scipy"
#  RESULT_VARIABLE scipy_import_result
#)
#if(NOT scipy_import_result EQUAL 0)
#  message(FATAL_ERROR "Python cannot import scipy during build")
#endif()


if (GENERATE_STUBS AND CMAKE_CROSSCOMPILING)
  message(WARNING "Stubs cannot be generated while cross-compiling - skipping.")
elseif(GENERATE_STUBS)
  # RECURSIVE didn't work keeps getting error when building/installing
  function(add_nanostub target marker underscore)
    if (underscore)
      set(MODNAME bobkit._${target})
    else()
      set(MODNAME bobkit.${target})
    endif()
    if (marker)
      nanobind_add_stub(
        bobkit_${target}_stub
        MODULE ${MODNAME}
        PYTHON_PATH ${py_dir}
        DEPENDS buc_py #project_gemmi
        OUTPUT ${py_dir}/bobkit/${target}.pyi
        )
    else()
      nanobind_add_stub(
        bobkit_stub
        MODULE bobkit
        PYTHON_PATH ${py_dir}
        MARKER_FILE ${py_dir}/bobkit/py.typed
        DEPENDS buc_py project_clipper project_gemmi
        OUTPUT ${py_dir}/bobkit/__init__.pyi
      )
    endif()
  endfunction()
  add_nanostub(bobkit false false)
  add_nanostub(clipper true false) 
  add_nanostub(buccaneer true true)
  add_nanostub(protein_db true false)
  add_nanostub(util true true)
  add_nanostub(data true false)
  add_nanostub(mm true false)

  #nanobind_add_stub(
  #  bobkit_stub
  #  MODULE bobkit
  #  PYTHON_PATH ${py_dir}
  #  DEPENDS buc_py project_gemmi
  #  MARKER_FILE ${py_dir}/bobkit/py.typed
  #  OUTPUT ${py_dir}/bobkit/__init__.pyi
  #  )
  #nanobind_add_stub(
  #  bobkit_clipper_stub
  #  MODULE bobkit.clipper
  #  PYTHON_PATH ${py_dir}
  #  DEPENDS buc_py project_gemmi
  #  OUTPUT ${py_dir}/bobkit/clipper.pyi
  #  )
  #nanobind_add_stub(
  #  bobkit_buccaneer_stub
  #  MODULE bobkit.buccaneer
  #  PYTHON_PATH ${py_dir}
  #  DEPENDS buc_py project_gemmi
  #  OUTPUT ${py_dir}/bobkit/buccaneer.pyi
  #  )
  #set(sep ":")
  #if (WIN32)
  #  set(setp ";")
  #endif()
  #add_custom_command(TARGET buc_py POST_BUILD COMMAND ${CMAKE_COMMAND} -E env "PYTHONPATH=$<TARGET_FILE_DIR:buc_py>${sep}$ENV{PYTHONPATH}"
  #${Python_EXECUTABLE} -m pybind11_stubgen bobkit
  #--enum-class-locations=Ignore)
endif()
# here 21Jan


if (BUILD_SHARED_LIBS)
  set_target_properties(buc_py PROPERTIES INSTALL_RPATH "${origin}/${rpath_rel}")
endif()
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/__init__.py"
              "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/buccaneer.py"
              "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/util.py"
              DESTINATION "${Python_SITEARCH}/bobkit" COMPONENT py)
#install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/_buccaneer.py" DESTINATION "${Python_SITEARCH}/bobkit" COMPONENT py)
#install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/_util.py" DESTINATION "${Python_SITEARCH}/bobkit" COMPONENT py)
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/_scripts" DESTINATION "${Python_SITEARCH}/bobkit"
        FILES_MATCHING PATTERN "*.py" )
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}/bobkit" COMPONENT py)
if (GENERATE_STUBS AND NOT CMAKE_CROSSCOMPILIING)
  install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/__init__.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/clipper.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/buccaneer.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/protein_db.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/util.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/data.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/mm.pyi"
                "${CMAKE_CURRENT_BINARY_DIR}/${py_dir}/bobkit/py.typed"
                DESTINATION "${Python_SITEARCH}/bobkit" COMPONENT py)
endif()

endif()
#this is for convenience to delete fftw
#set(INSTALL_SHAREDIR ${FFTW_INSTALL_PATH}/share/fftw3)
#install(CODE "string(REPLACE \";\" \"\\n\" MY_CMAKE_INSTALL_MANIFEST_CONTENT \"\$\{CMAKE_INSTALL_MANIFEST_FILES\}\")\n\
#            file(WRITE ${CMAKE_BINARY_DIR}/FFTW3_INSTALL_MANIFEST.txt \"\$\{MY_CMAKE_INSTALL_MANIFEST_CONTENT\}\\n${INSTALL_SHAREDIR}/FFTW3_INSTALL_MANIFEST.txt\")")
#install(FILES "${CMAKE_BINARY_DIR}/FFTW3_INSTALL_MANIFEST.txt" DESTINATION ${INSTALL_SHAREDIR})
