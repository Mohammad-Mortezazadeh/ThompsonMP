if( NOT CMAKE_BUILD_TYPE MATCHES "Debug" )
  add_definitions( -DNDEBUG )
endif( )
#add_definitions( -DSINGLE_PREC )

# Set flags depending on which compiler is used
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  include("compiler_flags_GNU_Fortran")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  include("compiler_flags_Intel_Fortran")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  include("compiler_flags_PGI_Fortran")
endif()

# Set flags depending on which compiler iis used
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
  include("compiler_flags_GNU_C")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  include("compiler_flags_Intel_C")
elseif (CMAKE_C_COMPILER_ID STREQUAL "PGI")
  include("compiler_flags_PGI_C")
endif()
