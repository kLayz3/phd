cmake_minimum_required(VERSION 3.16)
project(foot_simulations LANGUAGES CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

# Standalone telescope transport; deliberately excludes the reaction sources.
add_executable(foot_transport
    transport/foot_transport.cc src/TransportConfig.cc src/Telescope.cc
    ../includes/util/JSONParser.cxx ../includes/util/FilesystemUtils.cxx)
target_include_directories(foot_transport PRIVATE include ../includes)
target_compile_features(foot_transport PRIVATE cxx_std_17)
target_compile_options(foot_transport PRIVATE -Wall -Wextra)
target_compile_definitions(foot_transport PRIVATE FOOT_PARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../params")
target_link_libraries(foot_transport PRIVATE ${Geant4_LIBRARIES})

# Reuse the existing Nucleus implementation without loading libStructures/Go4.
find_package(ROOT REQUIRED COMPONENTS Core Tree Hist Graf Gpad Matrix)
target_sources(foot_transport PRIVATE
    ../includes/util/MPhysics.cxx ../includes/util/Geometry.cxx
    ../includes/util/FromChars.cxx ../includes/util/TextUtils.cxx)
target_include_directories(foot_transport PRIVATE ../includes/monad/indicators/include)
target_link_libraries(foot_transport PRIVATE ROOT::Core)

# ROOT-only postprocessing: Geant4 need not be rerun for trial resolutions.
add_executable(analyse_transport transport/analyse_transport.cc
    ../includes/util/JSONParser.cxx
    ../includes/util/CLI.cxx ../includes/util/TextUtils.cxx
    ../includes/util/FromChars.cxx ../includes/util/FilesystemUtils.cxx)
target_include_directories(analyse_transport PRIVATE
	../includes ../includes/monad/indicators/include)
target_compile_features(analyse_transport PRIVATE cxx_std_17)
target_compile_options(analyse_transport PRIVATE -Wall -Wextra)
target_compile_definitions(analyse_transport PRIVATE
    FOOT_PARAM_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../params"
    FOOT_TRACK_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/../quick/hit/cfg/foot_track_analysis/config.cfg")
target_link_libraries(analyse_transport PRIVATE ROOT::Core ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Matrix)

set_target_properties(foot_transport analyse_transport PROPERTIES
	RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/transport")

add_subdirectory(reactions)
