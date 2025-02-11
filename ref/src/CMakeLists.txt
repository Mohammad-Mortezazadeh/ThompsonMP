list( APPEND thompson_kernel_common_files
machine.F90
mt19937.F90
physcons.F90
mp_thompson.F90
mp_thompson_pre.F90
mp_thompson_post.F90
module_mp_radar.F90
module_mp_thompson_make_number_concentrations.F90
module_mp_thompson.F90
thompson_utils.F90
thompson_test.F90
)

add_executable( thompson_kernel_cpu ${thompson_kernel_common_files} )

if(ENABLE_GPU)
  add_executable( thompson_kernel_gpu ${thompson_kernel_common_files} )
  target_compile_options(thompson_kernel_gpu PUBLIC ${OpenACC_Fortran_OPTIONS}:cc${CUDA_ARCH_LIST},cuda${CUDAToolkit_VERSION_MAJOR}.${CUDAToolkit_VERSION_MINOR} -Minfo=accel)
  target_link_options(thompson_kernel_gpu PUBLIC ${OpenACC_Fortran_OPTIONS}:cc${CUDA_ARCH_LIST})
endif()

if(OpenMP_FOUND)
    target_link_libraries(thompson_kernel_cpu PUBLIC OpenMP::OpenMP_Fortran)
    if(ENABLE_GPU)
      target_link_libraries(thompson_kernel_gpu PUBLIC OpenMP::OpenMP_Fortran)
    endif()
endif()

if(MPI_FOUND)
    target_link_libraries(thompson_kernel_cpu PUBLIC MPI::MPI_Fortran)
    if(ENABLE_GPU)
      target_link_libraries(thompson_kernel_gpu PUBLIC MPI::MPI_Fortran)
    endif()
endif()

install(TARGETS thompson_kernel_cpu
        DESTINATION ${PROJECT_SOURCE_DIR}/exe)
if(ENABLE_GPU)
 install(TARGETS thompson_kernel_gpu
	 DESTINATION ${PROJECT_SOURCE_DIR}/exe)
endif()
