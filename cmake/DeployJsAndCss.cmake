include(FindPackageHandleStandardArgs)

if(NOT UglifyJS_EXECUTABLE)
    find_program(UglifyJS_EXECUTABLE uglifyjs)
endif()

find_package_handle_standard_args(UglifyJS DEFAULT_MSG UglifyJS_EXECUTABLE)


if(NOT UglifyCSS_EXECUTABLE)
    find_program(UglifyCSS_EXECUTABLE uglifycss)
endif()

find_package_handle_standard_args(UglifyCSS DEFAULT_MSG UglifyCSS_EXECUTABLE)


macro( deploy_js_resource input output )
message( "Will deploy JS ${input} to ${output}")
  if( UglifyJS_EXECUTABLE )
    message( "Will uglify ${input} to ${output} using ${UglifyJS_EXECUTABLE}")
    add_custom_command(
        OUTPUT ${output}
        COMMAND ${UglifyJS_EXECUTABLE} -c -o \"${output}\" \"${input}\" MAIN_DEPENDENCY ${input} )
  else( UglifyJS_EXECUTABLE )
    message( "Will COPY ${input} to ${output}")
    add_custom_command( OUTPUT ${output} COMMAND ${CMAKE_COMMAND} -E copy ${input} ${output} MAIN_DEPENDENCY ${input} )
  endif( UglifyJS_EXECUTABLE )
endmacro( deploy_js_resource )

macro( deploy_css_resource input output )
   message( "Will deploy CSS ${input} to ${output}")
  if( UglifyCSS_EXECUTABLE )
    message( "Will uglify ${input} to ${output} using ${UglifyCSS_EXECUTABLE}")
    add_custom_command(
        OUTPUT ${output}
        COMMAND ${UglifyCSS_EXECUTABLE} --output \"${output}\" \"${input}\" MAIN_DEPENDENCY ${input} )
  else( UglifyCSS_EXECUTABLE )
    message( "Will COPY ${input} to ${output}")
    add_custom_command( OUTPUT ${output} COMMAND ${CMAKE_COMMAND} -E copy ${input} ${output} MAIN_DEPENDENCY ${input} )
  endif( UglifyCSS_EXECUTABLE )
endmacro( deploy_css_resource )