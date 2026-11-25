find_package(ROOT REQUIRED COMPONENTS Physics Gui)
find_package(Python3 REQUIRED COMPONENTS Development)
set(util "${CMAKE_CURRENT_SOURCE_DIR}/../../includes/util")
add_executable(reactions main.cc Decay.cc Analysis.cc
	../src/Telescope.cc ../src/TransportConfig.cc
	${util}/PrettyHisto.cxx ${util}/MPhysics.cxx ${util}/Geometry.cxx ${util}/CLI.cxx
	${util}/MacroHelpers.cxx ${util}/JSONParser.cxx ${util}/FilesystemUtils.cxx ${util}/TextUtils.cxx ${util}/FromChars.cxx)
target_compile_features(reactions PRIVATE cxx_std_17)
target_compile_options(reactions PRIVATE -Wall -Wextra)
target_include_directories(reactions PRIVATE ../include ../../includes
	../../includes/monad/indicators/include ../../includes/pybind11/include)
target_compile_definitions(reactions PRIVATE
	REACTION_LIST="${CMAKE_CURRENT_SOURCE_DIR}/reaction_list.json")
target_link_libraries(reactions PRIVATE ${Geant4_LIBRARIES} ROOT::Physics ROOT::Gui
	ROOT::Hist ROOT::Gpad Python3::Python)

add_custom_command(TARGET reactions POST_BUILD
	COMMAND "${CMAKE_COMMAND}" -E create_symlink
		"$<TARGET_FILE:reactions>" "${CMAKE_CURRENT_SOURCE_DIR}/reactions"
	BYPRODUCTS "${CMAKE_CURRENT_SOURCE_DIR}/reactions"
	VERBATIM)
