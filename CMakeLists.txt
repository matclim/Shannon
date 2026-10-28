cmake_minimum_required(VERSION 3.21)
project(shannon VERSION 0.5.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Ensure plugin modules retain DT_NEEDED entries for all linked libraries,
# even if no symbols are directly referenced at link time (e.g. ROOT
# dictionary-only dependencies like SHiPDataModel's collection proxies).
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--no-as-needed")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} -Wl,--no-as-needed")

# Dependencies
find_package(SHiPDataModel REQUIRED)
find_package(SHiPRandom 0.4 REQUIRED)

if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW) # use BoostConfig.cmake, not the removed FindBoost
endif()
find_package(spdlog REQUIRED)
find_package(phlex 0.4.1 REQUIRED)
find_package(
    ROOT
    REQUIRED
    COMPONENTS Core RIO Hist MathCore ROOTNTuple
)

# The config file with versions etc
configure_file(ShannonConfig.h.in ShannonConfig.h)

# Sim file reader plugin
add_library(read_sim_file MODULE src/read_sim_file.cpp)
target_link_libraries(
    read_sim_file
    PRIVATE phlex::module SHiP::DataModel ROOT::RIO ROOT::ROOTNTuple
)

# Digitisation plugin
add_library(digitise_hits MODULE src/digitise_hits.cpp)
target_include_directories(
    digitise_hits
    PUBLIC "${PROJECT_BINARY_DIR}"
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(
    digitise_hits
    PRIVATE phlex::module SHiP::DataModel SHiP::Random
)

# Digitised output module plugin
add_library(digitised_output_module MODULE src/digitised_output_module.cpp)
target_link_libraries(
    digitised_output_module
    PRIVATE
        phlex::module
        SHiP::DataModel
        ROOT::Hist
        ROOT::RIO
        ROOT::ROOTHistUtil
        ROOT::ROOTNTuple
        spdlog::spdlog
)

# Smoke-test input generator (build-only, not installed)
# 2018 muon-flux test beam: drift tubes and RPC strips from aegir's
# workflows/muflux2018_replay.jsonnet (src/muflux/muflux_digitiser.hpp)
add_library(muflux_digitise MODULE src/muflux/muflux_digitise.cpp)
target_include_directories(muflux_digitise PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(
    muflux_digitise
    PRIVATE phlex::module SHiP::DataModel SHiP::Random spdlog::spdlog
)

add_library(muflux_digits_output MODULE src/muflux/muflux_digits_output.cpp)
target_include_directories(muflux_digits_output PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(
    muflux_digits_output
    PRIVATE phlex::module ROOT::ROOTNTuple spdlog::spdlog
)

enable_testing()
add_executable(muflux_digitiser_test tests/muflux_digitiser_test.cpp)
target_include_directories(muflux_digitiser_test PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")
add_test(NAME muflux_digitiser COMMAND muflux_digitiser_test)

add_executable(make_test_input tools/make_test_input.cpp)
target_include_directories(
    make_test_input
    PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src"
)
target_link_libraries(
    make_test_input
    PRIVATE SHiP::DataModel ROOT::RIO ROOT::ROOTNTuple SHiP::Random
)

# Install
include(GNUInstallDirs)

set(PLUGIN_TARGETS
    read_sim_file
    digitise_hits
    digitised_output_module
    muflux_digitise
    muflux_digits_output
)

install(TARGETS ${PLUGIN_TARGETS} LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY workflows/ DESTINATION share/shannon/workflows)
