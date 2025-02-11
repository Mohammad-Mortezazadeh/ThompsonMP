list( APPEND tool_files
  compare.sh
)

foreach(_file IN LISTS tool_files)
  execute_process( COMMAND ${CMAKE_COMMAND} -E copy
     ${CMAKE_CURRENT_SOURCE_DIR}/${_file}
     ${CMAKE_BINARY_DIR}/bin/thompson_${_file} )
endforeach()
