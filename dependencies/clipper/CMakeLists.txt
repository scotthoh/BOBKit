cmake_minimum_required(VERSION 3.15)
project(clipper CXX)
include(GNUInstallDirs)
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/source.cmake)
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads)
#option(BUILD_SHARED_LIBS "build shared libs (true)" TRUE)

# clipper-core
add_library(clipper-core ${clipper-core_sources} ${clipper-core_headers})
add_library(clipper::core ALIAS clipper-core)
add_dependencies(clipper-core project_fftw)
target_link_libraries(clipper-core PRIVATE ${FFTW3_LIBDIR}/${lib_fftw_name})

if (Threads_FOUND)
target_link_libraries(clipper-core PRIVATE Threads::Threads)
endif()
target_include_directories(clipper-core PRIVATE
                          ${WRK_DIR}/checkouts/clipper
                          ${WRK_DIR}/checkouts/clipper/clipper/core
                          ${FFTW3_INCDIR}
                          ${CMAKE_BINARY_DIR})
target_compile_options(clipper-core PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-core PROPERTIES PUBLIC_HEADER
                      "${WRK_DIR}/checkouts/clipper/clipper/clipper.h")
set_target_properties(clipper-core PROPERTIES PRIVATE_HEADER "${clipper-core_headers}")
#set_target_properties(clipper-core PROPERTIES SOVERSIO)
if(MSVC)
  set_target_properties(clipper-core PROPERTIES DEFINE_SYMBOL "CLIPPER_DLL_EXPORT")
endif()

install(TARGETS clipper-core
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/core
)

# clipper-contrib
add_library(clipper-contrib ${clipper-contrib_sources} ${clipper-contrib_headers})
add_library(clipper::contrib ALIAS clipper-contrib)
add_dependencies(clipper-contrib project_fftw)
target_link_libraries(clipper-contrib PRIVATE clipper::core ${FFTW3_LIBDIR}/${lib_fftw_name})

if(Threads_FOUND)
  target_link_libraries(clipper-contrib PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-contrib PRIVATE
                          ${WRK_DIR}/checkouts/clipper
                          ${WRK_DIR}/checkouts/clipper/clipper/core
                          ${WRK_DIR}/checkouts/clipper/clipper/contrib
                          ${FFTW3_INCDIR}
                          ${CMAKE_BINARY_DIR} )
target_compile_options(clipper-contrib PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-contrib PROPERTIES PUBLIC_HEADER
                      "${WRK_DIR}/checkouts/clipper/clipper/clipper-contrib.h")
set_target_properties(clipper-contrib PROPERTIES PRIVATE_HEADER "${clipper-contrib_headers}")

install(TARGETS clipper-contrib
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/contrib
)

add_library(clipper-gemmi ${clipper-gemmi_sources} ${clipper-gemmi_headers})
add_library(clipper::gemmi ALIAS clipper-gemmi)

target_link_libraries(clipper-gemmi PRIVATE clipper::core gemmi_cpp)
if(Threads_FOUND)
target_link_libraries(clipper-gemmi PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-gemmi PRIVATE
                          ${WRK_DIR}/checkouts/clipper
                          ${WRK_DIR}/checkouts/clipper/clipper/core
                          ${WRK_DIR}/checkouts/clipper/clipper/gemmi
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-gemmi PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-gemmi PROPERTIES PUBLIC_HEADER "${WRK_DIR}/checkouts/clipper/clipper/clipper-gemmi.h")
set_target_properties(clipper-gemmi PROPERTIES PRIVATE_HEADER "${clipper-gemmi_headers}")

install(TARGETS clipper-gemmi
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/gemmi
)

add_library(clipper-minimol ${clipper-minimol_headers} ${clipper-minimol_sources})
add_library(clipper::minimol ALIAS clipper-minimol)

target_link_libraries(clipper-minimol PRIVATE
                      clipper::core
                      clipper::gemmi
                      gemmi_cpp
)

if(Threads_FOUND)
target_link_libraries(clipper-minimol PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-minimol PRIVATE
                          ${WRK_DIR}/checkouts/clipper
                          ${WRK_DIR}/checkouts/clipper/clipper/core
                          ${WRK_DIR}/checkouts/clipper/clipper/gemmi
                          ${WRK_DIR}/checkouts/clipper/clipper/minimol
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-minimol PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-minimol PROPERTIES PUBLIC_HEADER "${WRK_DIR}/checkouts/clipper/clipper/clipper-minimol.h")
set_target_properties(clipper-minimol PROPERTIES PRIVATE_HEADER "${clipper-minimol_headers}")
if(MSVC)
  set_target_properties(clipper-minimol PROPERTIES DEFINE_SYMBOL "CLIPPER_MINIMOL_DLL_EXPORT")
endif()

install(TARGETS clipper-minimol
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/minimol
)

# clipper-cif
add_library(clipper-cif ${clipper-cif_headers} ${clipper-cif_sources})
add_library(clipper::cif ALIAS clipper-cif)

target_link_libraries(clipper-cif PRIVATE
                      clipper::core
                      gemmi_cpp
)

if(Threads_FOUND)
target_link_libraries(clipper-cif PRIVATE Threads::Threads)
endif()

target_include_directories(clipper-cif PRIVATE
                          ${WRK_DIR}/checkouts/clipper
                          ${WRK_DIR}/checkouts/clipper/clipper/core
                          ${WRK_DIR}/checkouts/clipper/clipper/cif
                          ${CMAKE_BINARY_DIR}
)
target_compile_options(clipper-cif PUBLIC "-DFFTW_ENABLE_FLOAT")
set_target_properties(clipper-cif PROPERTIES PUBLIC_HEADER "${WRK_DIR}/checkouts/clipper/clipper/clipper-cif.h")
set_target_properties(clipper-cif PROPERTIES PRIVATE_HEADER "${clipper-cif_headers}")

install(TARGETS clipper-cif
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include/clipper
        PRIVATE_HEADER DESTINATION include/clipper/cif
)
message(STATUS "clipper-core = " $<TARGET_FILE:clipper-core>)
message(STATUS "clipper-contrib = " $<TARGET_FILE:clipper-contrib>)
message(STATUS "clipper-gemmi = " $<TARGET_FILE:clipper-gemmi>)
message(STATUS "clipper-minimol =  " $<TARGET_FILE:clipper-minimol>)
message(STATUS "clipper-cif =  " $<TARGET_FILE:clipper-cif>)
