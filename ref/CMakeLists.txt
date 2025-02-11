################################################################################
# thompson
################################################################################

cmake_minimum_required( VERSION 3.16 )

option( ENABLE_GPU "Build to run on GPU" OFF )

if ( ENABLE_GPU )
  project( thompson LANGUAGES C Fortran CUDA )
else()
  project( thompson LANGUAGES C Fortran )
endif()

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )


enable_testing()

include(${PROJECT_NAME}_compiler_flags)

find_package( OpenMP COMPONENTS C Fortran )
find_package( MPI COMPONENTS C Fortran )

if ( ENABLE_GPU )
  set( OpenACC_ACCEL_TARGET tesla )
  find_package( OpenACC REQUIRED )
  find_package( CUDAToolkit REQUIRED )
  add_compile_definitions(ENABLE_GPU)
  include(FindCUDA/select_compute_arch)
  CUDA_DETECT_INSTALLED_GPUS(INSTALLED_GPU_CCS_1)
  string(STRIP "${INSTALLED_GPU_CCS_1}" INSTALLED_GPU_CCS_2)
  string(REPLACE " " ";" INSTALLED_GPU_CCS_3 "${INSTALLED_GPU_CCS_2}")
  list(GET INSTALLED_GPU_CCS_3 0 INSTALLED_GPU_CCS_4)
  string(REPLACE "." "" CUDA_ARCH_LIST "${INSTALLED_GPU_CCS_4}")
  message( STATUS "CUDA_ARCH_LIST: ${CUDA_ARCH_LIST}" )
  SET(CMAKE_CUDA_ARCHITECTURES ${CUDA_ARCH_LIST})
  add_compile_options("-Minfo=accel")
endif()

add_subdirectory(src)
add_subdirectory(test)

