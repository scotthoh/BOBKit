cmake_minimum_required(VERSION 3.15...3.22)
# check if called from python skbuild
#if (NOT DEFINED SKBUILD)
#  message(WARNING
#          "You should call this through Python so that all variables are set; pip install -v .")
#endif()

#file(STRINGS "${CMAKE_SOURCE_DIR}/src/cbuccaneer.cpp"
#    buccaneer_version_line REGEX "CCP4Program prog\\( \"cbuccaneer\",")
#string(FIND ${buccaneer_version_line} "\"," buccaneer_version_pos)
#string(SUBSTRING ${buccaneer_version_line} ${buccaneer_version_pos} 10 buccaneer_version)
#string(REGEX REPLACE "[,\" ]" "" buccaneer_version_str ${buccaneer_version})
#message(STATUS "cbuccaneer version ${buccaneer_version_str}")

#get version from python/version.hpp for python module version
file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/python/version.hpp"
      buildkit_version_line REGEX "#define BUILDKIT_VERSION \"")
string(REGEX REPLACE ".+\"([0-9]+\.[0-9]+\.[0-9]+)(-dev)?\"" "\\1"
      buildkit_version_str ${buildkit_version_line})
project(buildkit LANGUAGES CXX VERSION ${buildkit_version_str})
message(STATUS "Buildkit version ${PROJECT_VERSION}")

# building with scikit-build-core (pip install)
if (DEFINED SKBUILD)
  # wheel.install-dir="/data" is set in pyproject.toml
  # set paths to /platlib and /scripts here
  set(PYTHON_INSTALL_DIR "${SKBUILD_PLATLIB_DIR}")
endif()


# options
# Use -fPIC even if statically compiled
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
#set(BUILDKIT_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/buildkit" CACHE STRING
#  "Install path for buildkit CMake files")
set(THREADS_PREFER_PTHREAD_FLAG ON)
#option(INSTALL_DEV_FILES "Install C++ headesr, cmake configs, static lib" ON)
option(BUILD_SHARED_LIBS "Build as shared library" OFF)
#option(USE_PYTHON "Install python bindings" ON)
option(SHOW_COMPILE_TIME "Show compile time" OFF)
if(SHOW_COMPILE_TIME)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
endif()

# C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
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
# remove cbuccaneer.cpp and csequins.cpp
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/cbuccaneer.cpp")
list(REMOVE_ITEM BUCCANEER_SRC "${CMAKE_SOURCE_DIR}/include/buccaneer/csequins.cpp")
file(GLOB PY_SRC "python/*.cpp")

pybind11_add_module(buc_py MODULE ${PY_SRC} ${BUCCANEER_SRC})
target_include_directories(buc_py PRIVATE "dependencies/include")
target_include_directories(buc_py PRIVATE "include")
find_library(MMDB2DEP NAMES mmdb2)
find_library(CCP4DEP NAMES ccp4c)
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
message(STATUS "MMDB2  ${MMDB2DEP}")
message(STATUS "CCP4C ${CCP4CDEP}")
message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")

find_package(Threads REQUIRED)

target_link_libraries(buc_py PRIVATE ${MMDB2DEP}
                                     ${CCP4DEP}
                                     ${CLIPPERCOREDEP}
                                     ${CLIPPERMMDBDEP}
                                     ${CLIPPERMINIMOLDEP}
                                     ${CLIPPERCONTRIBDEP}
                                     ${CLIPPERCCP4DEP}
                                     Threads::Threads)
set_property(TARGET buc_py PROPERTY OUTPUT_NAME buildkit)
# install
if (DEFINED PYTHON_INSTALL_DIR)
  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
endif()
file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
#install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
install(TARGETS buc_py DESTINATION "${Python_SITEARCH}")


## Buccaneer library
##set(BUCCANEER_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/buccaneer" CACHE STRING 
##"Install path for Buccaneer library CMake files")
#
##find_package(Python3 3.8 REQUIRED COMPONENTS Interpreter Development)
#find_library(MMDB2DEP NAMES mmdb2
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
#find_library(CCP4CDEP NAMES ccp4c
#            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
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
#find_package(Threads REQUIRED)
##find_package(pybind11 REQUIRED)
#
## clipper-cif is needed only for cmakerefence, it will added explicitly
##list(REMOVE_ITEM CCP4_LIBRARIES ${CLIPPER-CIF_LIBRARY})
#message(STATUS "MMDB2  ${MMDB2DEP}")
#message(STATUS "CCP4C ${CCP4CDEP}")
#message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
#message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
#message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
#message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
#message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")
#
##message(STATUS "Python executable ${PYTHON_EXECUTABLE}")
##message(STATUS "PYTHONLIB ${PYTHON_LIBRARIES}")
#
#
#
## get cbuccaner version
#set (BUCCANEER_SRC_DIR ${CMAKE_SOURCE_DIR}/cpp/include/buccaneer)
##set (BUCCANEER_SRC_DIR ${CMAKE_SOURCE_DIR}/cpp/src)
#
#file(STRINGS "${BUCCANEER_SRC_DIR}/cbuccaneer.cpp"
#        buccaneer_version_line REGEX "CCP4Program prog\\( \"cbuccaneer\",")
#string(FIND ${buccaneer_version_line} "\"," buccaneer_version_pos)
#string(SUBSTRING ${buccaneer_version_line} ${buccaneer_version_pos} 10 buccaneer_version)
#string(REGEX REPLACE "[,\" ]" "" buccaneer_version_str ${buccaneer_version})
#message(STATUS "cbuccaneer version ${buccaneer_version_str}")
#file (GLOB SRC_FILES "${BUCCANEER_SRC_DIR}/*.cpp")
#list(REMOVE_ITEM SRC_FILES "${BUCCANEER_SRC_DIR}/cbuccaneer.cpp")
#file (GLOB HEADER_FILES "${BUCCANEER_SRC_DIR}/*.h")
## buc headers
#add_library(buccaneer_headers INTERFACE)
#add_library(buccaneer::headers ALIAS buccaneer_headers)
#
#target_include_directories(buccaneer_headers INTERFACE
#            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/include>"
#            "$<INSTALL_INTERFACE:cpp/include>")
#target_include_directories(buccaneer_headers INTERFACE
#            "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/dependencies/include>"
#            "$<INSTALL_INTERFACE:dependencies/include>")
#target_link_libraries(buccaneer_headers INTERFACE ${MMDB2DEP}
#                                        ${CCP4CDEP}
#                                        ${CLIPPERCOREDEP}
#                                        ${CLIPPERCCP4DEP}
#                                        ${CLIPPERMMDBDEP}
#                                        ${CLIPPERMINIMOLDEP}
#                                        ${CLIPPERCONTRIBDEP}
#                                        )
#target_link_libraries(buccaneer_headers INTERFACE Threads::Threads)
#target_compile_features(buccaneer_headers INTERFACE cxx_std_11)
#set_target_properties(buccaneer_headers PROPERTIES EXPORT_NAME headers)
###end buc headers
#
## for buccaneer library
##add_library(buccaneer_lib ${SRC_FILES})# ${HEADER_FILES})
##add_library(Buccaneer::buccaneer_lib ALIAS buccaneer_lib)
##set_property(TARGET buccaneer_lib PROPERTY POSITION_INDEPENDENT_CODE ON)
###set_target_properties(buccaneer_lib PROPERTIES PUBLIC_HEADER ${HEADER_FILES})
##
##
###target_compile_definitions(buccaneer_lib PUBLIC BUCCANEER_SHARED)
###target_link_libraries(buccaneer_lib PUBLIC buccaneer_headers)
###target_link_libraries(buccaneer_lib PRIVATE Threads::Threads)
##target_include_directories(buccaneer_lib PUBLIC
##                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cpp/src ${CMAKE_CURRENT_SOURCE_DIR}/cpp/include/buccaneer>"
##                           "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
##target_link_libraries(buccaneer_lib PRIVATE Threads::Threads)
#
## test cell
##add_executable(test_cell "${CMAKE_SOURCE_DIR}/cpp/test/test_cell.cpp")
#
#
##install(
##  TARGETS buccaneer_lib #buccaneer_headers #test_cell
##  EXPORT BuccaneerLibTargets
##  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
##  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
##  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
##  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
##)
##
### install exported targets
##install(
##  EXPORT BuccaneerLibTargets
##  FILE BuccaneerLibTargets.cmake
##  NAMESPACE Buccaneer::
##  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/BuccaneerLib
##)
#
## create a cmake package
#
#
## Prepare the Config.cmake.in content
#set(PACKAGE_INIT_MACRO "@PACKAGE_INIT@")
#set(CONFIG_CMAKE_IN "\
#@PACKAGE_INIT_MACRO@\n\
#include(\"\${CMAKE_CURRENT_LIST_DIR}/buildkit-targets.cmake\")\n\
#check_required_components(buildkit)\n\
#if(NOT buildkit_FIND_QUIETLY)\n\
#  message(STATUS \"Found buildkit version \${\${CMAKE_FIND_PACKAGE_NAME}_VERSION}\")\n"
#)
#
## Create Config.cmake.in
#file(CONFIGURE
#    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/buildkit-config.cmake.in
#    CONTENT ${CONFIG_CMAKE_IN}
#    @ONLY)
#
## Create BuccaneerLibConfig.cmake
#configure_package_config_file(
#    ${CMAKE_CURRENT_BINARY_DIR}/buildkit-config.cmake.in
#    ${CMAKE_CURRENT_BINARY_DIR}/buildkit-config.cmake
#    INSTALL_DESTINATION ${BUILDKIT_INSTALL_CMAKEDIR})
#
## Create BuccaneerLibConfigVersion.cmake
#write_basic_package_version_file(
#    "${CMAKE_CURRENT_BINARY_DIR}/buildkit-config-version.cmake"
#    VERSION "${PROJECT_VERSION}"
#    COMPATIBILITY AnyNewerVersion
#)
#
## Install CMake package files
## install
#if(INSTALL_DEV_FILES)
#  install(DIRECTORY cpp/include/buccaneer/ DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/buildkit")
#  install(TARGETS buccaneer_headers EXPORT BuildkitTargets)
#  install(EXPORT BuildkitTargets FILE buildkit-targets.cmake NAMESPACE buildkit::
#          DESTINATION "${BUILDKIT_INSTALL_CMAKEDIR}")
#  install(FILES
#    "${CMAKE_CURRENT_BINARY_DIR}/buildkit-config.cmake"
#    "${CMAKE_CURRENT_BINARY_DIR}/buildkit-config-version.cmake"
#    DESTINATION "${BUILDKIT_INSTALL_CMAKEDIR}"
#)
#endif()
#
#
## Python bindings
#if (USE_PYTHON)
#  
#  endif()
#
#file (GLOB PYTHON_FILES "python/*.cpp")
#file (GLOB PYTHON_FILES_H "python/*.h")
#
#pybind11_add_module(buildkit ${PYTHON_FILES})# ${PYTHON_FILES_H})
#
#
#target_link_libraries(buildkit PRIVATE ${MMDB2DEP}
#                                        ${CCP4CDEP}
#                                        ${CLIPPERCOREDEP}
#                                        ${CLIPPERCCP4DEP}
#                                        ${CLIPPERMMDBDEP}
#                                        ${CLIPPERMINIMOLDEP}
#                                        ${CLIPPERCONTRIBDEP}
#                                        )
#target_link_libraries(buildkit PRIVATE buccaneer_headers)
#target_link_libraries(buildkit PRIVATE Threads::Threads)
#
## The name of the shared library must match the module name
##set_target_properties(buildkit PROPERTIES
##    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/buccaneer_lib
##    OUTPUT_NAME "buildkit")
#set_property(TARGET buildkit PROPERTY OUTPUT_NAME buildkit)
#
## install
#if (DEFINED PYTHON_INSTALL_DIR)
#  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
#  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
#endif()
#file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
##install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
#install(TARGETS buildkit DESTINATION "${Python_SITEARCH}")
#endif()
#
### Directories with files to be linked
##include_directories(${CMAKE_SOURCE_DIR}/cpp/include/buccaneer
##                    ${CMAKE_SOURCE_DIR}/cpp/src
##                    ${CMAKE_SOURCE_DIR}/python
##                    ${CMAKE_SOURCE_DIR}/dependencies/include)
##
##
##
##                            
##target_link_libraries(buccaneer_lib ${MMDB2DEP}
##                                    ${CCP4CDEP}
##                                    ${CLIPPERCOREDEP}
##                                    ${CLIPPERMMDBDEP}
##                                    ${CLIPPERMINIMOLDEP}
##                                    ${CLIPPERCONTRIBDEP}
##                                    ${CLIPPERCCP4DEP}
##                                    ${CLIPEPRCIFDEP}
##                                    ${PYTHON_LIBRARY}
##)
##target_include_directories(buccaneer_lib PRIVATE ${BUCCANEER_HEADER_DIR})
##target_link_libraries(buccaneer_lib Threads::Threads)
###file (GLOB SOURCE_FILES "cpp/src/*.cpp")
###file (GLOB HEADER_FILES "cpp/include/buccaneer/*.h")
##file (GLOB PYTHON_FILES "python/*.cpp")
##file (GLOB PYTHON_FILES_H "python/*.h")
##### need to add the dependencies similar to privateer
##### install dependencies locally
##
##
##
##
##target_link_libraries(buildkit PRIVATE 
##                      buccaneer_lib
##                      ${MMDB2DEP}
##                      ${CCP4CDEP}
##                      ${CLIPPERCOREDEP}
##                      ${CLIPPERMMDBDEP}
##                      ${CLIPPERMINIMOLDEP}
##                      ${CLIPPERCONTRIBDEP}
##                      ${CLIPPERCCP4DEP}
##                      ${CLIPEPRCIFDEP}
##                      ${PYTHON_LIBRARY})
##
###include(CTest)
###enable_testing()
##
###add_library(buildkit buildkit.cpp)
##
###set(CPACK_PROJECT_NAME ${PROJECT_NAME})
###set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
###include(CPack)
##
##install(DIRECTORY ${BUCCANEER_HEADER_DIR} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
##install(TARGETS buccaneer_lib EXPORT BucTargets)
##install(EXPORT BucTargets FILE buccaneer-targets.cmake NAMESPACE buccaneer:: DESTINATION "${BUCCANEER_INSTALL_CMAKEDIR}")
##
##install ( TARGETS buccaneer_lib LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
##if (DEFINED PYTHON_INSTALL_DIR)
##  message(STATUS "Install directory for Python module: ${PYTHON_INSTALL_DIR}")
##  set(Python_SITEARCH "${PYTHON_INSTALL_DIR}")
##  file(TO_CMAKE_PATH "${Python_SITEARCH}" Python_SITEARCH)
##  #install(TARGETS buccaneer_cpp LIBRARY DESTINATION "${Python_SITEARCH}")
##  install(TARGETS buildkit DESTINATION "${Python_SITEARCH}")
##endif()
###install ( TARGETS buildkit DESTINATION ${PROJECT_SOURCE_DIR})
##