cmake_minimum_required(VERSION 3.16)
project(c9be9_fragmentation LANGUAGES CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

file(GLOB SOURCES CONFIGURE_DEPENDS src/*.cc)

add_executable(c9be9 main.cc ${SOURCES})
target_include_directories(c9be9 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(c9be9 PRIVATE ${Geant4_LIBRARIES})
target_compile_features(c9be9 PRIVATE cxx_std_17)
target_compile_options(c9be9 PRIVATE -Wall -Wextra -Wpedantic)

configure_file(run.mac run.mac COPYONLY)

# Standalone telescope transport; deliberately excludes the reaction sources.
add_executable(foot_transport
    transport/main.cc transport/Config.cc
    ../includes/util/JSONParser.cxx ../includes/util/FilesystemUtils.cxx)
target_include_directories(foot_transport PRIVATE transport ../includes)
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
add_executable(foot_transport_analysis transport/analyse_transport.cc
    ../includes/util/CLI.cxx ../includes/util/TextUtils.cxx
    ../includes/util/FromChars.cxx ../includes/util/FilesystemUtils.cxx)
target_include_directories(foot_transport_analysis PRIVATE ../includes)
target_compile_features(foot_transport_analysis PRIVATE cxx_std_17)
target_compile_options(foot_transport_analysis PRIVATE -Wall -Wextra)
target_compile_definitions(foot_transport_analysis PRIVATE MND_EXPORT_OSTREAM_G_NAMESPACE
    FOOT_TRACK_CONFIG="${CMAKE_CURRENT_SOURCE_DIR}/../quick/hit/cfg/foot_track_analysis/config.cfg")
target_link_libraries(foot_transport_analysis PRIVATE ROOT::Core ROOT::Tree ROOT::Hist ROOT::Graf ROOT::Gpad ROOT::Matrix)
