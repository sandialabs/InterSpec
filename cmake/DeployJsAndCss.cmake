include(FindPackageHandleStandardArgs)

# Skip minification for local dev builds (local-server / OSX app) and for the Xcode generator, where
#  node typically isn't on the build tool's sanitized PATH - terser/uglifyjs are `#!/usr/bin/env node`
#  scripts, so `env: node: No such file or directory` breaks the build.  Deployment/web builds still
#  minify so shipped resources stay small.  When skipped, the macros below fall back to a plain copy.
if( BUILD_AS_LOCAL_SERVER OR BUILD_AS_OSX_APP OR (CMAKE_GENERATOR STREQUAL "Xcode") )
    message( "Skipping JS/CSS minification (local/OSX/Xcode build): resources will be copied verbatim")
else()
    if(NOT UglifyJS_EXECUTABLE)
        # NAMES is required to search for more than one program.  In find_program's short signature,
        #  `find_program(<VAR> name1 [path1 ...])`, only the first token is a name and the rest are
        #  taken as search *paths* - so without NAMES this looked for `terser` alone and treated
        #  `uglifyjs` as a directory, silently shipping un-minified JS wherever only uglify-js is
        #  installed (which is what target/electron/build_linux_app_from_docker.sh installs).
        find_program(UglifyJS_EXECUTABLE NAMES terser uglifyjs)
    endif()

    find_package_handle_standard_args(UglifyJS DEFAULT_MSG UglifyJS_EXECUTABLE)


    if(NOT UglifyCSS_EXECUTABLE)
        find_program(UglifyCSS_EXECUTABLE uglifycss)
    endif()

    find_package_handle_standard_args(UglifyCSS DEFAULT_MSG UglifyCSS_EXECUTABLE)
endif()


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