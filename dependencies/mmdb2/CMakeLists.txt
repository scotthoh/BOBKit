cmake_minimum_required(VERSION 3.15)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/mmdb2/mmdb_defs.h"
    mmdb_vers_line REGEX "")
#option(INSTALL_DEV_FILES "Install C++ headers, cmake configs, static lib" OFF)
# c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_RPATH_USE_ORIGIN ON)
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
if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 -Wredundant-decls -Wfloat-conversion -Wdisabled-optimization -Wshadow $ENV{EXTRA_WFLAGS}")
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" cxx_flags_config)
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS} ${${cxx_flags_config}}")
endif()

# automatically parse the version number
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/mmdb2/mmdb_defs.h" mmdb2_version_file)
string(REGEX MATCH "MAJOR_VERSION[ = ]+([0-9]+)" mmdb2_major_match "${mmdb2_version_file}")
set(MMDB2_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "MINOR_VERSION[ = ]+([0-9]+)" mmdb2_minor_match "${mmdb2_version_file}")
set(MMDB2_MINOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "MICRO_VERSION[ = ]+([0-9]+)" mmdb2_micro_match "${mmdb2_version_file}")
set(MMDB2_MICRO_VERSION "${CMAKE_MATCH_1}")
set(MMDB2_VERSION_NUMBER ${MMDB2_MAJOR_VERSION}.${MMDB2_MINOR_VERSION}.${MMDB2_MICRO_VERSION})

project(mmdb2 LANGUAGES CXX VERSION ${MMDB2_VERSION_NUMBER})
include(GNUInstallDirs)
message(STATUS "mmdb2 version ${MMDB2_VERSION_NUMBER}")
option(SHOW_COMPILE_TIME "Show compile time" ON)
if(SHOW_COMPILE_TIME)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "\"${CMAKE_COMMAND}\" -E time")
endif()
set(MMDB2_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/mmdb2" CACHE STRING
    "Install path for mmdb2 CMake files")
# c++ standard
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_BUILD_TYPE Release)
# Common installation directories, CMake helpers
option(BUILD_SHARED_LIBS "Build shared library" ON)
option(INSTALL_DEV_FILES "Install headers, cmake configs" ON)
#set(MMDB_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/mmdb2" CACHE STRING "Install path for mmdb2 CMake files")

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU" AND EXTRA_WARNINGS)
  set(CMAKE_CXX_FLAGS
      "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -Wformat=2 -Wredundant-decls -Wfloat-conversion -Wshadow $ENV{EXTRA_WFLAGS}")
  string(TOUPPER "CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}" cxx_flags_config)
  message(STATUS "C++ flags set to: ${CMAKE_CXX_FLAGS} ${${cxx_flags_config}}")
endif()

#add_library(mmdb2 INTERFACE)
message(STATUS "CURRENT CMAKE SOURCE DIR ${CMAKE_CURRENT_SOURCE_DIR}")
file(GLOB mmdb_src "${CMAKE_CURRENT_SOURCE_DIR}/mmdb2/*.cpp")
#message(STATUS ${mmdb_src})
file(GLOB mmdb_headers "${CMAKE_CURRENT_SOURCE_DIR}/mmdb2/*.h")
add_library(mmdb2 ${mmdb_src} ${mmdb_headers})
add_library(mmdb2::mmdb2 ALIAS mmdb2)
#target_include_directories(mmdb2 PRIVATE #PUBLIC #INTERFACE
#    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mmdb2>"
#    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>")
target_include_directories(mmdb2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mmdb2/*.cpp)
#target_compile_features(mmdb2 INTERFACE cxx_std_11)
#set_target_properties(mmdb2 PROPERTIES EXPORT_NAME MMDB2)
target_compile_options(mmdb2 PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(mmdb2 PROPERTIES PUBLIC_HEADER "${mmdb_headers}")
if (INSTALL_DEV_FILES)
  install(TARGETS mmdb2 #EXPORT mmdb2-export
          #ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
          #LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib #${CMAKE_INSTALL_LIBDIR}
          LIBRARY DESTINATION lib #${CMAKE_INSTALL_LIBDIR}
          #RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
          PUBLIC_HEADER DESTINATION include/mmdb2) #${CMAKE_INSTALL_INCLUDEDIR}/mmdb2)
          #PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include/mmdb2) #${CMAKE_INSTALL_INCLUDEDIR}/mmdb2)
  #install(EXPORT mmdb2-export FILE mmdb2-targets.cmake NAMESPACE mmdb2::
  #        DESTINATION "${MMDB2_INSTALL_CMAKEDIR}")
  #install(FILES "${CMAKE_CURRENT_BINARY_DIR}/mmdb2-config.cmake"
  #              "${CMAKE_CURRENT_BINARY_DIR}/mmdb2-config-version.cmake"
  #              DESTINATION "${MMDB2_INSTALL_CMAKEDIR}")
endif()


#include(CMakePackageConfigHelpers)
#configure_package_config_file()
#file(GLOB MMDB2_SOURCE "mmdb2/*.cpp")
#file(GLOB MMDB2_HEADER "mmdb2/*.h")

#target_include_directories()
