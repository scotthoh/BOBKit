cmake_minimum_required(VERSION 3.15...3.22)
project(buildkit VERSION 0.1.2)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CLIBENV $ENV{CLIB})
message(STATUS "CLIBENV ${CLIBENV}")
# copied from cmakelist clipper-progs
#find_package(CCP4 REQUIRED mmdb2 ccp4c clipper-core clipper-ccp4 clipper-mmdb2
#                            clipper-contrib clipper-minimol)
find_library(MMDB2DEP NAMES mmdb2
            PATHS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
find_library(CCP4CDEP NAMES ccp4c
            HINTS ${CMAKE_SOURCE_DIR}/dependencies/lib ${CLIBENV})
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
find_package(Threads REQUIRED)
find_package(pybind11 REQUIRED)
# clipper-cif is needed only for cmakerefence, it will added explicitly
#list(REMOVE_ITEM CCP4_LIBRARIES ${CLIPPER-CIF_LIBRARY})
message(STATUS "MMDB2  ${MMDB2DEP}")
message(STATUS "CCP4C ${CCP4CDEP}")
message(STATUS "CLIPPERCORE ${CLIPPERCOREDEP}")
message(STATUS "CLIPPERMMDB ${CLIPPERMMDBDEP}")
message(STATUS "CLIPPERMINI ${CLIPPERMINIMOLDEP}")
message(STATUS "CLIPPERCONT ${CLIPPERCONTRIBDEP}")
message(STATUS "CLIPPERCCP4 ${CLIPPERCCP4DEP}")
message(STATUS "CLIPPERCIF ${CLIPEPRCIFDEP}")
message(STATUS "PYTHONLIB ${PYTHON_LIBRARY}")

option(BUILD_SHARED_LIBS "Build as shared library" ON)
set(THREADS_PREFER_PTHREAD_FLAG ON)

# Directories with files to be linked
include_directories(${CMAKE_SOURCE_DIR}/cpp/include/buccaneer
                    ${CMAKE_SOURCE_DIR}/cpp/src
                    ${CMAKE_SOURCE_DIR}/python
                    ${CMAKE_SOURCE_DIR}/dependencies/include)

set (BUCCANEER_HEADER_DIR ${CMAKE_SOURCE_DIR}/cpp/include/buccaneer)
set (BUCCANEER_SRC_DIR ${CMAKE_SOURCE_DIR}/cpp/src)

# get cbuccaner version
file(STRINGS "${BUCCANEER_SRC_DIR}/cbuccaneer.cpp"
        buccaneer_version_line REGEX "CCP4Program prog\\( \"cbuccaneer\",")
string(FIND ${buccaneer_version_line} "\"," buccaneer_version_pos)
string(SUBSTRING ${buccaneer_version_line} ${buccaneer_version_pos} 10 buccaneer_version)
string(REGEX REPLACE "[,\" ]" "" buccaneer_version_str ${buccaneer_version})
message(STATUS "cbuccaneer version ${buccaneer_version_str}")

# for buccaneer library
add_library(buccaneer_lib   ${BUCCANEER_SRC_DIR}/buccaneer-build.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-prune.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-lib.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-correct.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-sequence.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-link.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-filter.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-util.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-ncsbuild.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-find.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-prep.cpp
                            ${BUCCANEER_SRC_DIR}/simplex-lib.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-grow.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-prot.cpp
                            ${BUCCANEER_SRC_DIR}/simulate-lib.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-join.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-known.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-merge.cpp
                            ${BUCCANEER_SRC_DIR}/buccaneer-tidy.cpp)
                            
target_link_libraries(buccaneer_lib ${MMDB2DEP}
                                    ${CCP4CDEP}
                                    ${CLIPPERCOREDEP}
                                    ${CLIPPERMMDBDEP}
                                    ${CLIPPERMINIMOLDEP}
                                    ${CLIPPERCONTRIBDEP}
                                    ${CLIPPERCCP4DEP}
                                    ${CLIPEPRCIFDEP}
                                    ${PYTHON_LIBRARY}
)
target_include_directories(buccaneer_lib PUBLIC $BUCCANEER_HEADER_DIR)
target_link_libraries(buccaneer_lib ${MMDB2DEP}
                                    ${CCP4CDEP}
                                    ${CLIPPERCOREDEP}
                                    ${CLIPPERMMDBDEP}
                                    ${CLIPPERMINIMOLDEP}
                                    ${CLIPPERCONTRIBDEP}
                                    ${CLIPPERCCP4DEP}
                                    ${CLIPEPRCIFDEP}
                                    ${PYTHON_LIBRARY}
                                    Threads::Threads)
#file (GLOB SOURCE_FILES "cpp/src/*.cpp")
#file (GLOB HEADER_FILES "cpp/include/buccaneer/*.h")
file (GLOB PYTHON_FILES "python/*.cpp")
file (GLOB PYTHON_FILES_H "python/*.h")
### need to add the dependencies similar to privateer
### install dependencies locally


pybind11_add_module(buildkit ${PYTHON_FILES} ${PYTHON_FILES_H})

target_link_libraries(buildkit PRIVATE 
                      buccaneer_lib
                      ${MMDB2DEP}
                      ${CCP4CDEP}
                      ${CLIPPERCOREDEP}
                      ${CLIPPERMMDBDEP}
                      ${CLIPPERMINIMOLDEP}
                      ${CLIPPERCONTRIBDEP}
                      ${CLIPPERCCP4DEP}
                      ${CLIPEPRCIFDEP}
                      ${PYTHON_LIBRARY})



#include(CTest)
#enable_testing()

#add_library(buildkit buildkit.cpp)

#set(CPACK_PROJECT_NAME ${PROJECT_NAME})
#set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
#include(CPack)
install ( TARGETS buccaneer_lib LIBRARY DESTINATION ${PROJECT_SOURCE_DIR})
install ( TARGETS buildkit DESTINATION ${PROJECT_SOURCE_DIR})
