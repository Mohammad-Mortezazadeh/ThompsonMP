####################################################################
# COMMON FLAGS
####################################################################
set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -ffp-contract=off ")

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -march=native -funroll-all-loops -finline-functions")

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -O0 -g -fcheck=bounds -ffpe-trap=invalid,zero,overflow,underflow" )

####################################################################
# FLAGS FOR GPU
####################################################################

set( Fortran_GPU_FLAGS        "" )
