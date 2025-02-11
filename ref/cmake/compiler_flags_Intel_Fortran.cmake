####################################################################
# COMMON FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -fp-model precise " )

####################################################################
# RELEASE FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -xHost " )

####################################################################
# DEBUG FLAGS
####################################################################

set( CMAKE_Fortran_FLAGS_DEBUG   "-O0 -debug -nolib-inline -fno-inline-functions -assume protect_parens,minus0 -prec-div -prec-sqrt -check bounds -check uninit -fp-stack-check -init=snan,array -warn unused" )

####################################################################
# FLAGS FOR GPUS
####################################################################

set( Fortran_GPU_FLAGS        "" )
