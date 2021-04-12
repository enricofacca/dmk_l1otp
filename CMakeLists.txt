#
# Main CMakeLists.txt file to build library "dmk_l1otp"
#

cmake_minimum_required( VERSION 3.10 )

message("**")
message("** PROJECT: dmk_l1opt")

# set compiler type
if ( NOT DEFINED My_Fortran_Compiler)
  SET( CMAKE_Fortran_COMPILER /usr/bin/gfortran )
else ()
  SET( CMAKE_Fortran_COMPILER ${My_Fortran_Compiler})
endif()
MESSAGE( "** Compiling with ${CMAKE_Fortran_COMPILER}")


# set default 
if ( NOT DEFINED BUILD_TYPE)
  set (BUILD_TYPE "RELEASE")
endif()
set (CMAKE_Fortran_FLAGS_BASE "-fPIC -ffree-line-length-none")
set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_BASE} -O5")
set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_BASE} -g -C -Wall -fcheck=all -O")
set (CMAKE_Fortran_FLAGS_RelWithIEEE "${CMAKE_Fortran_FLAGS_REALESE} -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")

# set build type (DEBUG/RELEASE)
if ( ${BUILD_TYPE} MATCHES "DEBUG" )
  message("** Compiling in debug mode")
  set ( CMAKE_BUILD_TYPE Debug )
elseif ( ${BUILD_TYPE} MATCHES "RELEASE"   )
  message("** Compiling in release mode")
  set ( CMAKE_BUILD_TYPE RELEASE )
elseif ( ${BUILD_TYPE} MATCHES "IEEE"   )
  message("** Compiling in release with IEEE mode")
  set ( CMAKE_BUILD_TYPE RelWithIEEE )
else ()
  set ( CMAKE_BUILD_TYPE user )
  set ( CMAKE_Fortran_FLAGS_user  "${CMAKE_Fortran_FLAGS_BASE} ${BUILD_TYPE}" )
  message("** Compiling in user mode:  ${CMAKE_Fortran_FLAGS_user}")
endif()

# set the project name
project( dmk_l1opt LANGUAGES Fortran VERSION 1.0 )


#####################
# chiama e costruisce la "dmk_solver"
if( NOT TARGET "dmk_solver" )
  # Look in the parent directory
  FILE( GLOB local ${CMAKE_CURRENT_SOURCE_DIR}/../dmk_solver )
  if (EXISTS ${local})
    MESSAGE( "** Repository p1galerkin found in parent directory" )
    set ( DMK_SOLVER_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../dmk_solver )
    set ( DMK_SOLVER_BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}../dmk_solver/build)
    ADD_SUBDIRECTORY( ${DMK_SOLVER_SRC_DIR} ${DMK_SOLVER_BIN_DIR})
  else()
     # print error message
    MESSAGE( SEND_ERROR "++++++++++++++++++++++++++++++++++++++++++++++")
    MESSAGE( SEND_ERROR  "Repository p1galerkin NOT FOUND" )
    MESSAGE( SEND_ERROR  "Check README for instructions" )
    MESSAGE( FATAL_ERROR "++++++++++++++++++++++++++++++++++++++++++++++")
  endif()
endif()
#####################

# include paths containg .mod files 
# and append the new ones 
include_directories( ${MODS} )                         
LIST( APPEND MODS ${MODS_BINARY_DIR} )
set(MODS ${MODS} PARENT_SCOPE )


# Use only specific one target
#ADD_DEPENDENCIES(optimalnetworkdesing f2py_interface_dmk)
add_executable (dmk_p1p0_from_file${EXECUTABLE_EXTENSION} srcs/main.f90)
target_link_libraries( dmk_p1p0_from_file${EXECUTABLE_EXTENSION}
  dmk p1galerkin geometry linalg agmg arpack blas lapack)

# LOOK FOR BLAS AND LIBRARY
# SETTING VARIABLES:
# BLAS_LIB, BLAS_DIR
#set(CMAKE_FIND_LIBRARY_SUFFIXES .a .so)
if ( NOT DEFINED BLAS_LIB OR NOT DEFINED BLAS_DIR ) 
  find_package(BLAS)
  message(STATUS "BLAS LIBRARIES: ${BLAS_LIBRARIES}")
  list(GET BLAS_LIBRARIES 0 BLAS_LIB)
endif()
get_filename_component(BLAS_DIR "${BLAS_LIB}" DIRECTORY)
message("BLAS LIB: ${BLAS_LIB} dir ${BLAS_DIR}")

# LOOK FOR LAPACK LIBRARY
# SETTING VARIABLES:
# LAPACK_LIB, LAPACK_DIR
if ( NOT DEFINED LAPACK_LIB)
  message(STATUS "LOOKING FOR LAPACK LIBRARIES")
  find_package(LAPACK REQUIRED)
  message(STATUS "LAPACK LIBRARIES: ${LAPACK_LIBRARIES}")
  list(GET LAPACK_LIBRARIES 0 LAPACK_LIB)
endif()
get_filename_component(LAPACK_DIR "${LAPACK_LIB}" DIRECTORY)
message(STATUS "LAPACK LIB: ${LAPACK_LIB} ")
message("LAPACK LIB: ${LAPACK_LIB}")

add_custom_target(local ALL DEPENDS f2py_interface_dmk timedata2vtk.out)




message("**")
