if (CMAKE_HOST_WIN32)
    # need at least CMake 3.23 for IntelLLVM support of IntelDPCPP package on Windows
    cmake_minimum_required(VERSION 3.23)
else()
    # CMake 3.20.5 is the minimum recommended for IntelLLVM on Linux
    cmake_minimum_required(VERSION 3.20.5)
endif()

cmake_policy(SET CMP0091 NEW)
if(MSVC AND MSVC_VERSION LESS "1920")
	message(FATAL_ERROR "MCsquare requires Visual Studio Version 2019 with openmp:experimental option to compile with msvc")
endif()

project(MCSquare C)
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>") # CACHE STRING "Runtime Library Setting")


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

# This is needed to avoid that gcc uses a builtin fmin that wrongly handles NaNs
if (CMAKE_C_COMPILER_ID MATCHES "GNU")
	add_compile_options("-fno-builtin-fmin" "-fno-builtin-fmax")
endif() 

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

if(MCSQUARE_USE_DOUBLE_PRECISION_COMPUTE)
	set(VAR_COMPUTE_PRECISION 2)
else()
	set(VAR_COMPUTE_PRECISION 1)
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
        src/include/lib/win_dirent.h
)

set_property(TARGET MCsquare PROPERTY C_STANDARD 99)

set(MKL_LINK "static" CACHE STRING "MKL Link Type")
if(NOT MINGW)
	find_package(MKL CONFIG)
endif()

if (MKL_FOUND)
	option(USE_MKL_LIB "Use the Intel MKL Libraries" ON)
endif(MKL_FOUND)

if (MKL_FOUND AND USE_MKL_LIB)
    target_compile_options(MCsquare PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
    target_include_directories(MCsquare PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(MCsquare PUBLIC $<LINK_ONLY:MKL::MKL>)
endif(MKL_FOUND AND USE_MKL_LIB)

#Find OpenMP package and determine version... not working well at the moment (Version not read)
find_package(OpenMP)
if (OpenMP_C_FOUND)
    # Enable appropriate OpenMP flags based on compiler
    if(MSVC)
        message(STATUS "Configuring OpenMP with simd vectorization for MSVC compiler id ${CMAKE_C_COMPILER_ID} (version ${MSVC_VERSION}) with toolset ${MSVC_TOOLSET_VERSION}")

        if(CMAKE_C_COMPILER_ID MATCHES "Intel")
            set(OpenMP_C_FLAGS "/Qopenmp" CACHE STRING "C Compiler Flags for OpenMP Parallelization" FORCE)
        elseif(MSVC_VERSION GREATER_EQUAL 1920)
            set(OpenMP_C_FLAGS "/openmp:experimental" CACHE STRING "C Compiler Flags for OpenMP Parallelization" FORCE)
            message(STATUS "Adding MSVC>=2019 specific /openmp:experimental flag for enabling OpenMP SIMD pragmas")
        else()
            message("No OpenMP simd pragmas available in OpenMP Version ${OpenMP_VERSION} for compiler ${CMAKE_C_COMPILER_ID}")
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

include(CheckIncludeFile)
include(CheckIncludeFiles)
check_include_file("mathimf.h" HAVE_MATHIMF_H)
check_include_file("math.h" HAVE_MATH_H)
check_include_file("BaseTsd.h" HAVE_BASETSD_H)
check_include_file("sys/stat.h" HAVE_SYS_STAT_H)
check_include_file("direct.h" HAVE_DIRECT_H)
check_include_file("dirent.h" HAVE_DIRENT_H)

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

if(HAVE_DIRECT_H)
  SET(MKDIR_INCLUDES "
  #include <sys/types.h>
  #include <sys/stat.h>
  #include <direct.h>
  ")
else()
    SET(MKDIR_INCLUDES "
    #include <sys/types.h>
    #include <sys/stat.h>
    ")
endif()

check_c_source_compiles("${MKDIR_INCLUDES} int main() { return mkdir(\"test\", 0755); }" HAVE_MKDIR_TWO_ARGS)
check_c_source_compiles("${MKDIR_INCLUDES} int main() { return mkdir(\"test\"); }" HAVE_MKDIR_ONE_ARG)
check_c_source_compiles("${MKDIR_INCLUDES} int main() { return _mkdir(\"test\", 0755); }" HAVE__MKDIR_TWO_ARGS)
check_c_source_compiles("${MKDIR_INCLUDES} int main() { return _mkdir(\"test\"); }" HAVE__MKDIR_ONE_ARG)

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


install(TARGETS MCsquare RUNTIME DESTINATION ".")
install(DIRECTORY Scanners DESTINATION ".")
install(DIRECTORY Materials DESTINATION ".")
install(DIRECTORY BDL DESTINATION ".")
if (WIN32)
  # This generator expression expands to full paths of runtime DLLs.
  install(FILES $<TARGET_RUNTIME_DLLS:MCsquare> DESTINATION "." COMPONENT Runtime)
  if (OMP_DLL_DIR)
    install(FILES ${OMP_DLL_DIR}/libiomp5md.dll DESTINATION "." COMPONENT Runtime)
  endif()
endif()

# Enable CTest, but keep it easy to turn off with -DBUILD_TESTING=OFF
include(CTest)           # defines BUILD_TESTING and calls enable_testing()

if (BUILD_TESTING)
  # Run: <path-to-myapp> -v
  add_test(NAME MCsquare.version
           COMMAND $<TARGET_FILE:MCsquare> -v
           WORKING_DIRECTORY $<TARGET_FILE_DIR:MCsquare>)

  # Build a regex that matches your version string
  # (escape dots so 1.2.3 is treated literally)
  string(REGEX REPLACE "\\." "\\\\." _VER_RX "${PROJECT_VERSION}")

  set_tests_properties(MCsquare.version PROPERTIES
    PASS_REGULAR_EXPRESSION "(${PROJECT_NAME}[[:space:]]+)?${_VER_RX}"
    TIMEOUT 5)

  # If MKL was used set the properties
  if(MKL_ENV)
    set_tests_properties(MCsquare.version PROPERTIES ENVIRONMENT "${MKL_ENV}")
  endif()
endif()

include(CPack)
set(CPACK_GENERATOR "ZIP")

# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
    # message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()
