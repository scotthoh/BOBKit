cmake_minimum_required(VERSION 3.15...3.22)
set(gemmi_src ${WRK_DIR}/checkouts/gemmi)
# get version string from version.hpp
file(STRINGS "${gemmi_src}/include/gemmi/version.hpp"
     gemmi_version_line REGEX "#define GEMMI_VERSION \"")
string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
       gemmi_version_str ${gemmi_version_line})

project(gemmi LANGUAGES C CXX VERSION ${gemmi_version_str})
message(STATUS "Gemmi version ${PROJECT_VERSION}")

include(GNUInstallDirs)  # for CMAKE_INSTALL_LIBDIR, etc

if (DEFINED SKBUILD)  # building with scikit-build-core (pip install)
  # we set wheel.install-dir="/data" in pyproject.toml,
  # so here we need to only set paths to /platlib and /scripts
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
  #set(CMAKE_INSTALL_BINDIR "${SKBUILD_SCRIPTS_DIR}")
endif()

option(BUILD_SHARED_LIBS "Build shared library" OFF)
#option(BUILD_GEMMI_PROGRAM "Build gemmi command-line program" OFF)
option(INSTALL_DEV_FILES "Install C++ headers, cmake configs, static lib" ON)
option(USE_FORTRAN "Build Fortran bindings" OFF)
option(USE_PYTHON "Build Python bindings" ON)
option(EXTRA_WARNINGS "Set extra warning flags" OFF)
option(USE_WMAIN "(Windows only) take Unicode arguments in gemmi program" ON)
set(GEMMI_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/gemmi" CACHE STRING
    "Install path for gemmi CMake files")
option(STRIP_BINARY "Strip symbols from program" OFF)

# uncomment to show compilation times for each compilation unit
#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")

if (DEFINED ENV{FC} OR CMAKE_Fortran_COMPILER)
  set(USE_FORTRAN ON CACHE BOOL "Build Fortran bindings" FORCE)
endif()

if (NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
  set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)
endif()
if (NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
  set(CMAKE_CXX_VISIBILITY_PRESET hidden)
endif()

# We don't need to set policy CMP0069 for this file, because the minimum
# version is 3.15 and for versions >= 3.9 it's automatically set to NEW.
# But a FindPython script in cmake 3.15 have "cmake_policy(VERSION 3.7)"
# so we need to set CMP0069 because of that.
set(CMAKE_POLICY_DEFAULT_CMP0069 NEW)

# Unless it's explicitly specifed, we want LTO for non-debug builds only.
# pybind11Common.cmake has its own logic to set LTO options.
# Having different options in pybind11_add_module than for gemmi_cpp
# may cause problems -- we may need to handle it at some point.
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

if (INFO)
  set(GEMMI_VERSION_INFO ${INFO} CACHE STRING "Extra text for gemmi -V" FORCE)
endif()

if (USE_FORTRAN)
   message(STATUS "Currently, Fortran bindings are not built by cmake.")
   message(STATUS "They can be built by running make in fortran directory.")
  enable_language(Fortran)
else()
  #message(STATUS "Skipping Fortran bindings. Add -D USE_FORTRAN=1 to build them.")
endif()

if (DEFINED ENV{EXTRA_WFLAGS})
  set(EXTRA_WARNINGS ON CACHE BOOL "Set extra warning flags" FORCE)
endif()

set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(CheckIncludeFile)  # for check_include_file

if (DEFINED ENV{CXXFLAGS})
  set(USING_ENV_CXXFLAGS ON CACHE BOOL "" FORCE)
endif()

# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

find_package(ZLIB)
if (ZLIB_FOUND)
  include_directories("${ZLIB_INCLUDE_DIR}")
else()
  message(STATUS "The build will use zlib code from third_party/zlib.")
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/third_party/zlib")
endif()

if (NOT DEFINED SKBUILD AND CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(benchmark 1.3 QUIET)
endif()
if (benchmark_FOUND)
  message(STATUS "Found benchmark: ${benchmark_DIR}")
else (NOT benchmark_FOUND)
  message(STATUS "Benchmarks not configured.")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 -Wredundant-decls -Wfloat-conversion -Wdisabled-optimization -Wshadow $ENV{EXTRA_WFLAGS}")
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" cxx_flags_config)
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS} ${${cxx_flags_config}}")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Qvec-report:1")
endif()

if (USE_FORTRAN)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=f2003 -fimplicit-none -Wall -Wextra -pedantic")
    message(STATUS "Fortran flags set to: ${CMAKE_Fortran_FLAGS}")
    set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -fbounds-check")
  endif()
endif()


if (ZLIB_FOUND)
  macro(support_gz exe)
    target_link_libraries(${exe} PRIVATE ZLIB::ZLIB)
  endmacro()
else()
  add_library(ungz OBJECT
      "third_party/zlib/adler32.c"
      "third_party/zlib/crc32.c"
      "third_party/zlib/gzlib.c"
      "third_party/zlib/gzread.c"
      "third_party/zlib/inflate.c"
      "third_party/zlib/inftrees.c"
      "third_party/zlib/inffast.c"
      "third_party/zlib/zutil.c")
  check_include_file(unistd.h has_unistd_h)
  target_compile_definitions(ungz PRIVATE NO_GZCOMPRESS=1)
  if (has_unistd_h)
    target_compile_definitions(ungz PRIVATE Z_HAVE_UNISTD_H=1)
  endif()
  if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    # /wd4996 - the POSIX name for this item is deprecated.
    # /wd4267 - conversion from 'size_t' to 'unsigned int', possible loss of data
    target_compile_options(ungz PRIVATE /wd4996 /wd4267)
  endif()
  set_property(TARGET ungz PROPERTY POSITION_INDEPENDENT_CODE ON)
  macro(support_gz exe)
    target_sources(${exe} PUBLIC $<TARGET_OBJECTS:ungz>)
  endmacro()
endif()

if (WIN32 AND USE_WMAIN)
  add_definitions(-D_UNICODE=1)
endif()

#add_library(cgemmi STATIC fortran/grid.cpp fortran/symmetry.cpp)
#
#if (USE_FORTRAN)
#  add_library(fgemmi STATIC fortran/gemmi.f90)
#  target_link_libraries(fgemmi PRIVATE cgemmi)
#endif()
#include_directories("${gemmi_src}/include"
#        "${gemmi_src}/third_party")

file(GLOB gemmi_CPPSRC "${gemmi_src}/src/*.cpp")
file(GLOB gemmi_HEADERS "${gemmi_src}/include/gemmi/*.hpp")
file(GLOB gemmi_third_party-headers_HEADERS "${gemmi_src}/include/gemmi/third_party/*.h")
file(GLOB gemmi_py_SRC "${gemmi_src}/python/*.cpp")

add_library(gemmi_headers INTERFACE)
add_library(gemmi::headers ALIAS gemmi_headers)
target_include_directories(gemmi_headers INTERFACE
    "$<BUILD_INTERFACE:${gemmi_src}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_compile_features(gemmi_headers INTERFACE cxx_std_11)
set_target_properties(gemmi_headers PROPERTIES EXPORT_NAME headers)

add_library(gemmi_cpp ${gemmi_CPPSRC})
add_library(gemmi::gemmi_cpp ALIAS gemmi_cpp)
set_property(TARGET gemmi_cpp PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_property(TARGET gemmi_cpp PROPERTY CXX_VISIBILITY_PRESET hidden)
target_compile_definitions(gemmi_cpp PRIVATE GEMMI_BUILD)
target_include_directories(gemmi_cpp PRIVATE "${gemmi_src}/third_party")

if (BUILD_SHARED_LIBS)
  target_compile_definitions(gemmi_cpp PUBLIC GEMMI_SHARED)
endif()
target_link_libraries(gemmi_cpp PUBLIC gemmi_headers)
support_gz(gemmi_cpp)

### Python bindings ###
if (USE_PYTHON)
  message(STATUS "The python module will be built.")
  # CMake >=3.18 has subcomponent Development.Module, scikit-build-core also has it
  if (${CMAKE_VERSION} VERSION_LESS 3.18 AND NOT SKBUILD)
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development)
  else()
    find_package(Python ${PYTHON_VERSION} REQUIRED COMPONENTS Interpreter Development.Module)
  endif()
  if (EXISTS "${CMAKE_HOME_DIRECTORY}/pybind11")
    message(STATUS "Using ${CMAKE_HOME_DIRECTORY}/pybind11 (internal copy).")
    add_subdirectory(pybind11)
  else()
    # use pybind11-config (if available) to determine pybind11_DIR
    execute_process(COMMAND pybind11-config --cmakedir OUTPUT_VARIABLE pybind11_DIR)
    string(STRIP "${pybind11_DIR}" pybind11_DIR)
    find_package(pybind11 2.6 CONFIG REQUIRED)
    message(STATUS "Found pybind11 ${pybind11_VERSION}: ${pybind11_INCLUDE_DIRS}")
  endif()
  pybind11_add_module(gemmi_py ${gemmi_py_SRC} $<TARGET_OBJECTS:gemmi_cpp>)
  target_link_libraries(gemmi_py PRIVATE gemmi_headers)
  set_property(TARGET gemmi_py PROPERTY OUTPUT_NAME gemmi)
  if(CMAKE_CXX_FLAGS MATCHES "-Wshadow")
    target_compile_options(gemmi_py PRIVATE "-Wno-shadow")
  endif()
  support_gz(gemmi_py)
else()
  message(STATUS "Skipping Python module. Add -D USE_PYTHON=1 to build it.")
endif()

# avoid GCC warning: the ABI of passing structure with 'complex float' member
# has changed in GCC 4.4
set_property(SOURCE ${gemmi_src}/python/recgrid.cpp ${gemmi_src}/python/hkl.cpp ${gemmi_src}/python/mtz.cpp ${gemmi_src}/prog/mtz.cpp
             PROPERTY COMPILE_FLAGS $<$<CXX_COMPILER_ID:GNU>:-Wno-psabi>)

include(CMakePackageConfigHelpers)
configure_package_config_file(${gemmi_src}/tools/gemmi-config.cmake.in gemmi-config.cmake
                              INSTALL_DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
write_basic_package_version_file(gemmi-config-version.cmake
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

if (INSTALL_DEV_FILES)
  # In CMake 3.23 we'll be able to use FILE_SET to install the headers
  install(DIRECTORY ${gemmi_src}/include/gemmi DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
  install(TARGETS gemmi_headers gemmi_cpp EXPORT GemmiTargets)
  install(EXPORT GemmiTargets FILE gemmi-targets.cmake NAMESPACE gemmi::
          DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
  install(FILES
          "${CMAKE_CURRENT_BINARY_DIR}/gemmi-config.cmake"
          "${CMAKE_CURRENT_BINARY_DIR}/gemmi-config-version.cmake"
          DESTINATION "${GEMMI_INSTALL_CMAKEDIR}")
elseif (BUILD_SHARED_LIBS)
  install(TARGETS gemmi_cpp)
endif()

if (USE_PYTHON)
  if (DEFINED PYTHON_INSTALL_DIR)
    message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
    set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
  endif()
  file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
  # Using Python_SITEARCH/SITELIB is not good, because they are absolute
  # and don't respect CMAKE_INSTALL_PREFIX.
  # https://discourse.cmake.org/t/findpython3-how-to-specify-local-installation-directory-for-python-module/3580/5
  install(TARGETS gemmi_py DESTINATION "${Python_SITEARCH}")
endif()
