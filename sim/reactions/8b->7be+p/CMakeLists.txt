add_executable(reaction_8b_to_7be_p
	main.cc
	src/ActionInitialization.cc
	src/C9FragmentationProcess.cc
	src/EventAction.cc
	src/FragmentationPhysics.cc
	src/PhysicsList.cc
	src/RunAction.cc
	src/SteppingAction.cc)
target_include_directories(reaction_8b_to_7be_p PRIVATE
	include "${PROJECT_SOURCE_DIR}/../includes")
target_link_libraries(reaction_8b_to_7be_p PRIVATE sim_setup)
target_compile_options(reaction_8b_to_7be_p PRIVATE -Wall -Wextra -Wpedantic)
set_target_properties(reaction_8b_to_7be_p PROPERTIES OUTPUT_NAME run)
configure_file(run.mac run.mac COPYONLY)
target_compile_definitions(reaction_8b_to_7be_p PRIVATE
	REACTION_DEFAULT_MACRO="${CMAKE_CURRENT_BINARY_DIR}/run.mac")
