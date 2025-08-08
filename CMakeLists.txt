if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for IntelLLVM support of IntelDPCPP package on Windows
    cmake_minimum_required(VERSION 3.23)
else()
    # CMake 3.20.5 is the minimum recommended for IntelLLVM on Linux
    cmake_minimum_required(VERSION 3.20.5)
endif()

project(MCsquare 
	LANGUAGES C 
	VERSION 1.1)
	
include(cmake/MCSquareGitState.cmake)

#mark_as_advanced(CMAKE_MSVC_RUNTIME_LIBRARY)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/CMakeModules")

option(ENABLE_FAST_MATH "Optimize with fast floating point math" ON)
if(ENABLE_FAST_MATH)
	if(MSVC)
		add_compile_options("/fp:fast")
	else()
		add_compile_options("-ffast-math")
	endif()
else()
	if(MSVC)
		add_compile_options("/fp:precise")
	endif()
endif()

add_executable(MCsquare "")
#set_target_properties(MCsquare PROPERTIES MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT MCsquare)

#Add Source Files
target_sources(MCsquare
	PRIVATE
		src/main.c
		src/struct.c
		src/compute_4D.c
		src/compute_beam_model.c
		src/compute_DVH.c
		src/compute_EM_interaction.c
		src/compute_Fippel.c
		src/compute_geometry.c
		src/compute_hadron.c
		src/compute_math.c
		src/compute_nuclear_interaction.c
		src/compute_random.c
		src/compute_range_shifter.c
		src/compute_scenario_simulation.c
		src/compute_scoring.c
		src/compute_semi_infinite_slab.c
		src/compute_simulation.c
		src/compute_simulation_beamlet.c
		src/compute_treatment_uncertainties.c
		src/data_4D.c
		src/data_ascii.c
		src/data_beam_model.c
		src/data_config.c
		src/data_ct.c
		src/data_ct_mhd.c
		src/data_ct_penelope.c
		src/data_materials.c
		src/data_mhd.c
		src/data_nuclear.c
		src/data_sparse.c
		src/data_Stop_Pow.c
		src/data_contours.c
		src/File_process.c
		src/include/struct.h
		src/include/compute_4D.h
		src/include/compute_beam_model.h
		src/include/compute_DVH.h
		src/include/compute_EM_interaction.h
		src/include/compute_Fippel.h
		src/include/compute_geometry.h
		src/include/compute_hadron.h
		src/include/compute_math.h
		src/include/compute_nuclear_interaction.h
		src/include/compute_random.h
		src/include/compute_range_shifter.h
		src/include/compute_scenario_simulation.h
		src/include/compute_scoring.h
		src/include/compute_semi_infinite_slab.h
		src/include/compute_simulation.h
		src/include/compute_simulation_beamlet.h
		src/include/compute_treatment_uncertainties.h
		src/include/data_4D.h
		src/include/data_ascii.h
		src/include/data_beam_model.h
		src/include/data_config.h
		src/include/data_ct.h
		src/include/data_ct_mhd.h
		src/include/data_ct_penelope.h
		src/include/data_materials.h
		src/include/data_mhd.h
		src/include/data_nuclear.h
		src/include/data_sparse.h
		src/include/data_Stop_Pow.h
		src/include/data_contours.h
		src/include/define.h
		src/include/File_process.h
)

set_property(TARGET MCsquare PROPERTY C_STANDARD 99)

find_package(IntelDPCPP REQUIRED)

#try to find MKL. Path has usually
find_package(MKL)
if (MKL_FOUND)
	option(USE_MKL_LIB "Use the Intel MKL Libraries" ON)
	
	if(USE_MKL_LIB AND CMAKE_C_COMPILER_ID STREQUAL "Intel")
		if(MSVC)
			target_compile_options(MCsquare PUBLIC "/Qmkl")
		else()
			target_compile_options(MCsquare PUBLIC "-Qmkl")
		endif()
		target_link_libraries(MCsquare PUBLIC "${MKL_LIBRARIES_PARALLEL}") 
	elseif(USE_MKL_LIB)
		include_directories("${MKL_INCLUDE}")
		target_link_libraries(MCsquare PUBLIC "${MKL_LIBRARIES_SEQUENTIAL}") #The thread library relies on intels OpenMP
	endif()
endif(MKL_FOUND)

#Find OpenMP package and determine version... not working well at the moment (Version not read)
find_package(OpenMP)
if (OpenMP_C_FOUND)
	#We can only use simd pragmas with visual studio 2019 (v142) and the openmp:experimental switch
	if (MSVC AND MSVC_VERSION LESS "1920" AND NOT (CMAKE_C_COMPILER_ID STREQUAL "Intel"))
		message("${MSVC_TOOLSET_VERSION}")
		message("No OpenMP simd pragmas available in OpenMP Version ${OpenMP_VERSION}!")
	else()
		if(MSVC AND NOT (CMAKE_C_COMPILER_ID STREQUAL "Intel"))
			set(OpenMP_C_FLAGS "/openmp:experimental" CACHE STRING "C Compiler Flags for OpenMP Parallelization" FORCE)
			message(STATUS "Adding MSVC-2019 specific /openmp:experimental flag for enabling OpenMP SIMD pragmas")
		elseif(MSVC AND (CMAKE_C_COMPILER_ID STREQUAL "Intel"))
			set(OpenMP_C_FLAGS "/Qopenmp" CACHE STRING "C Compiler Flags for OpenMP Parallelization" FORCE)
		else()
		endif()
	endif()
	target_compile_options(MCsquare PUBLIC ${OpenMP_C_FLAGS})
	target_link_libraries(MCsquare PUBLIC OpenMP::OpenMP_C)	
endif(OpenMP_C_FOUND)

#Explicitly link math library on unix systems
if(UNIX AND NOT (CMAKE_C_COMPILER_ID STREQUAL "Intel"))
	target_link_libraries(MCsquare PUBLIC m)
endif()

# I guess this could be extended to manage the definitions much more smoothly
#include(CheckSymbolExists)
#check_symbol_exists("rand_r" "<stdlib.h>" HAVE_RAND_R)
#check_symbol_exists("strtok_r" "<string.h>" HAVE_STRTOK_R)
#if(NOT HAVE_STRTOK_R)
#	check_symbol_exists("strtok_s" "<string.h>" HAVE_STRTOK_S)
#endif()

#Options to set PRECISION
option(MCSQUARE_USE_DOUBLE_PRECISION_DATA "Work with double precision data in MCsquare" OFF)
option(MCSQUARE_USE_DOUBLE_PRECISION_SCORING "Score double precision values in the computation" OFF)
option(MCSQUARE_USE_DOUBLE_PRECISION_COMPUTE "Use double precision throught MC computations" OFF)

if(MCSQUARE_USE_DOUBLE_PRECISION_DATA)
	set(VAR_DATA_PRECISION 2)
else()
	set(VAR_DATA_PRECISION 1)
endif()

if(MCSQUARE_USE_DOUBLE_PRECISION_SCORING)
	set(VAR_SCORING_PRECISION 2)
else()
	set(VAR_SCORING_PRECISION 1)
endif()

if(MCSQUARE_USE_DOUBLE_COMPUTE_SCORING)
	set(VAR_COMPUTE_PRECISION 2)
else()
	set(VAR_COMPUTE_PRECISION 1)
endif()

include(CheckIncludeFile)
include(CheckIncludeFiles)
check_include_file("mathimf.h" HAVE_MATHIMF_H)
check_include_file("math.h" HAVE_MATH_H)
check_include_file("BaseTsd.h" HAVE_BASETSD_H)

check_include_files("stdio.h;stdlib.h;time.h;string.h;float.h;stdint.h" HAVE_REQUIRED_HEADERS)

if(NOT HAVE_REQUIRED_HEADERS)
  message(FATAL_ERROR "A required header file could not be found for your compiler!")
endif()

include(CheckSymbolExists)
check_symbol_exists(strtok_r "stdio.h;string.h;" HAVE_STRTOK_R)
check_symbol_exists(strtok_s "stdio.h;string.h;" HAVE_STRTOK_S)

include(CheckCSourceCompiles)
check_c_source_compiles("void *aa_test(void *x) { return __assume_aligned(x, 16);}\nint main(void) { return 0; }" HAVE_ASSUME_ALIGNED)
check_c_source_compiles("void *aa_test(void *x) { return __builtin_assume_aligned(x, 16);}\nint main(void) { return 0; }" HAVE_BUILTIN_ASSUME_ALIGNED)

check_c_source_compiles("float aa_test() { __attribute__((aligned(64))) float x[16]; return x[0];} \nint main(void) { return 0; }" HAVE_ATTRIBUTE_ALIGNED)
check_c_source_compiles("float aa_test() { __declspec(align(64)) float x[16]; return x[0];} \nint main(void) { return 0; }" HAVE_DECLSPEC_ALIGNED)

include(CheckTypeSize)
check_type_size("ssize_t" SIZEOF_SSIZE_T)


if  (WIN32)
  set(RMDIR_CMD "rd /s /q  %s")
else()
  set(RMDIR_CMD "rm -r  %s")
endif()

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/define.h.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/include/define_cmake.h)
target_include_directories(MCsquare PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include/)
target_sources(MCsquare PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/include/define_cmake.h)

#target_compile_definitions(MCsquare PRIVATE $<$<BOOL:${MCSQUARE_USE_DOUBLE_PRECISION_DATA}>:USE_DOUBLE_PRECISION_DATA> $<$<BOOL:${MCSQUARE_USE_DOUBLE_PRECISION_SCORING}>:USE_DOUBLE_PRECISION_SCORING> $<$<BOOL:${MCSQUARE_USE_DOUBLE_PRECISION_COMPUTE}>:USE_DOUBLE_PRECISION_COMPUTE>)
target_compile_definitions(MCsquare PRIVATE CMAKE_BUILD)

install(TARGETS MCsquare DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY Scanners DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY Materials DESTINATION ${CMAKE_INSTALL_PREFIX})
install(DIRECTORY BDL DESTINATION ${CMAKE_INSTALL_PREFIX})

#include(CPack)




# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
    # message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
