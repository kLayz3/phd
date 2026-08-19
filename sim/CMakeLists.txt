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
