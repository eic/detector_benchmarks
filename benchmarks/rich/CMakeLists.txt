cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

get_filename_component(det_name ${CMAKE_CURRENT_LIST_DIR} NAME)
project(${CMAKE_PROJECT_NAME}_${det_name} LANGUAGES CXX)

cmake_policy(SET CMP0074 NEW) # use `<PackageName>_ROOT` variables
include(GNUInstallDirs)

# ROOT
find_package(ROOT 6 REQUIRED COMPONENTS Core RIO Hist)
include(${ROOT_USE_FILE})

# data model
find_package(podio REQUIRED)
find_package(EDM4HEP REQUIRED)
find_package(EDM4EIC REQUIRED)

# logging
find_package(spdlog REQUIRED)
add_compile_definitions(SPDLOG_FMT_EXTERNAL)

# sources
set(algo_exe_source ${CMAKE_CURRENT_SOURCE_DIR}/src/benchmark.cc)
file(GLOB algo_headers CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB algo_sources CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
list(REMOVE_ITEM algo_sources ${algo_exe_source})

# library
set(algo_lib ${PROJECT_NAME})
add_library(${algo_lib} SHARED ${algo_sources})
target_include_directories(${algo_lib} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(${algo_lib} PRIVATE -Wall -Wno-misleading-indentation)

# executable
set(algo_exe benchmark_${det_name}_reconstruction)
add_executable(${algo_exe} ${algo_exe_source})

# linking
target_link_libraries(${algo_lib}
  PUBLIC
  ROOT::Core
  ROOT::Hist
  podio::podio
  podio::podioRootIO
  EDM4EIC::edm4eic
  EDM4HEP::edm4hep
  spdlog::spdlog
  )
target_link_libraries(${algo_exe}
  PRIVATE
  ${algo_lib}
  )

# installation
install(FILES ${algo_headers} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}/${det_name})
install(TARGETS ${algo_exe} ${algo_lib})
