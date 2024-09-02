cmake_minimum_required(VERSION 3.15...3.22)
project(clipper CXX)
include(GNUInstallDirs)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/source.cmake)
set(THREADS_PREFER_PTHREAD_FLAG ON)
# Set default build mode (based on CMake FAQ)
if (NOT CMAKE_BUILD_TYPE AND NOT USING_ENV_CXXFLAGS)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
find_package(Threads)
#find_package(gemmi 0.6.4 CONFIG REQUIRED)
#if(gemmi_FOUND)
#  message(STATUS gemmi::gemmi_cpp)
#endif()
if(NOT APPLE)
 set(CMAKE_INSTALL_RPATH $ORIGIN)
endif()
if(FFTW_LIBPATH)
  message(STATUS "fftw libpath = " ${FFTW_LIBPATH})
endif()
option(BUILD_SHARED_LIBS "Build shared library" ON)
if(GEMMI_LIBPATH)
message(STATUS "gemmi libpath = " ${GEMMI_LIBPATH})
if(BUILD_SHARED_LIBS)
  add_library(fftw::fftw SHARED IMPORTED)
  add_library(gemmi::gemmi_cpp SHARED IMPORTED)
else()
  add_library(gemmi::gemmi_cpp STATIC IMPORTED)
  add_library(fftw::fftw STATIC IMPORTED)
endif()
set_target_properties(gemmi::gemmi_cpp PROPERTIES IMPORTED_LOCATION ${GEMMI_LIBPATH})
set_target_properties(fftw::fftw PROPERTIES IMPORTED_LOCATION ${FFTW_LIBPATH})
endif()
# clipper-core
add_library(clipper-core ${clipper-core_sources} ${clipper-core_headers})
add_library(clipper::core ALIAS clipper-core)
#add_dependencies(clipper-core project_fftw)
target_link_libraries(clipper-core PRIVATE fftw::fftw) #${FFTW_LIBPATH})

if (Threads_FOUND)
target_link_libraries(clipper-core PRIVATE Threads::Threads)
endif()
target_include_directories(clipper-core PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR}
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/core
                          #${FFTW3_INCDIR}
                          ${FFTW_INCPATH}
                          ${CMAKE_BINARY_DIR})
target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER
                     "${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper.h")
                      #"${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper/clipper.h")
set_target_properties(clipper-core PROPERTIES PRIVATE_HEADER "${clipper-core_headers}")
set_property(TARGET clipper-core PROPERTY POSITION_INDEPENDENT_CODE ON)
#set_target_properties(clipper-core PROPERTIES SOVERSIO)
if(MSVC)
  set_target_properties(clipper-core PROPERTIES DEFINE_SYMBOL "CLIPPER_DLL_EXPORT")
endif()
set_target_properties(clipper-core PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS clipper-core
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/core
)

# clipper-contrib
add_library(clipper-contrib ${clipper-contrib_sources} ${clipper-contrib_headers})
add_library(clipper::contrib ALIAS clipper-contrib)
#add_dependencies(clipper-contrib project_fftw)
target_link_libraries(clipper-contrib PRIVATE clipper::core fftw::fftw)# ${FFTW_LIBPATH})#${FFTW3_LIBDIR}/${lib_fftw_name})

if(Threads_FOUND)
  target_link_libraries(clipper-contrib PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-contrib PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR} #/clipper
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/core #/clipper/clipper/core
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/contrib #/clipper/clipper/contrib
                          ${FFTW_INCPATH}
                          #${FFTW3_INCDIR}
                          ${CMAKE_BINARY_DIR} )
target_compile_options(clipper-contrib PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-contrib PROPERTIES PUBLIC_HEADER
                     "${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper-contrib.h")
                      #"${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper/clipper-contrib.h")
set_target_properties(clipper-contrib PROPERTIES PRIVATE_HEADER "${clipper-contrib_headers}")
set_property(TARGET clipper-contrib PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(clipper-contrib PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS clipper-contrib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/contrib
)

add_library(clipper-gemmi ${clipper-gemmi_sources} ${clipper-gemmi_headers})
add_library(clipper::gemmi ALIAS clipper-gemmi)
#add_dependencies(clipper-gemmi project_gemmi)

target_link_libraries(clipper-gemmi PRIVATE clipper::core fftw::fftw)# $<TARGET_FILE:gemmi_cpp>) #${GEMMI_LIBPATH})
if(Threads_FOUND)
target_link_libraries(clipper-gemmi PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-gemmi PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR}                   #/clipper
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/core      #/clipper/core
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/gemmi  #/clipper/gemmi
                          ${GEMMI_INCPATH}
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-gemmi PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper-gemmi.h") #clipper/clipper-gemmi.h")
set_target_properties(clipper-gemmi PROPERTIES PRIVATE_HEADER "${clipper-gemmi_headers}")
set_property(TARGET clipper-gemmi PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(clipper-gemmi PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS clipper-gemmi
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/gemmi
)

add_library(clipper-minimol ${clipper-minimol_headers} ${clipper-minimol_sources})
add_library(clipper::minimol ALIAS clipper-minimol)
#add_dependencies(clipper-minimol project_gemmi)

target_link_libraries(clipper-minimol PRIVATE
                      clipper::core
                      clipper::gemmi
                      gemmi::gemmi_cpp #${GEMMI_LIBPATH}
                      #gemmi::gemmi_cpp
)

if(Threads_FOUND)
target_link_libraries(clipper-minimol PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-minimol PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR}                           #/clipper
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/core              #clipper/core
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/gemmi             #clipper/gemmi
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/minimol           #clipper/minimol
                          ${GEMMI_INCPATH}
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-minimol PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper-minimol.h") #clipper/clipper-minimol.h")
set_target_properties(clipper-minimol PROPERTIES PRIVATE_HEADER "${clipper-minimol_headers}")
set_property(TARGET clipper-minimol PROPERTY POSITION_INDEPENDENT_CODE ON)
if(MSVC)
  set_target_properties(clipper-minimol PROPERTIES DEFINE_SYMBOL "CLIPPER_MINIMOL_DLL_EXPORT")
endif()
set_target_properties(clipper-minimol PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS clipper-minimol
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/minimol
)

# clipper-cif
add_library(clipper-cif ${clipper-cif_headers} ${clipper-cif_sources})
add_library(clipper::cif ALIAS clipper-cif)
#add_dependencies(clipper-cif project_gemmi)

target_link_libraries(clipper-cif PRIVATE
                      clipper::core
                      gemmi::gemmi_cpp
                      #"${GEMMI_LIBPATH}"
                      #gemmi::gemmi_cpp
)

if(Threads_FOUND)
target_link_libraries(clipper-cif PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-cif PRIVATE
                          ${CMAKE_CURRENT_SOURCE_DIR}                           #/clipper
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/core              #clipper/core
                          ${CMAKE_CURRENT_SOURCE_DIR}/clipper/cif               #clipper/cif
                          ${GEMMI_INCPATH}
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-cif PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-cif PROPERTIES PUBLIC_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/clipper/clipper-cif.h")#clipper/clipper-cif.h")
set_target_properties(clipper-cif PROPERTIES PRIVATE_HEADER "${clipper-cif_headers}")
set_property(TARGET clipper-cif PROPERTY POSITION_INDEPENDENT_CODE ON)
set_target_properties(clipper-cif PROPERTIES INSTALL_RPATH "$ORIGIN")
install(TARGETS clipper-cif
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/cif
)