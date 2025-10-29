cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project( InterSpecApp VERSION 1.0.13 )

if(POLICY CMP0114)
  # Add @rpath to the install_name of shared libraries and frameworks that are part of a macOS application bundle.
  cmake_policy(SET CMP0114 NEW)
endif()

set(BUILD_AS_OSX_APP ON CACHE INTERNAL "Build a OSX native app" )
set(USE_OSX_NATIVE_MENU ON CACHE INTERNAL "Use macOS native menu-system" )
set(USE_SPECRUM_FILE_QUERY_WIDGET ON CACHE INTERNAL "Enable spectrum file query widget" )
set(TRY_TO_STATIC_LINK ON CACHE INTERNAL "Try to statically link" )
set(CMAKE_OSX_DEPLOYMENT_TARGET "10.15" CACHE INTERNAL "Set OSX deployment target to 10.15" )
set(USE_REL_ACT_TOOL ON CACHE INTERNAL "Enable Rel. Act. Tools" )
set(USE_LEAFLET_MAP ON CACHE INTERNAL "Use leaflet map tool" )
option( USE_BATCH_TOOLS "Enables batch processing." ON )

add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../../ InterSpecLib )


SET(MACOSX_BUNDLE_INFO_STRING   "InterSpec - Sandia National Lab, Will Johnson")
SET(MACOSX_BUNDLE_ICON_FILE InterSpec )
SET(MACOSX_BUNDLE_LONG_VERSION_STRING "" )
SET(MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )
SET(MACOSX_BUNDLE_SHORT_VERSION_STRING "InterSpec" )
SET(MACOSX_BUNDLE_BUNDLE_VERSION "0" )
SET(MACOSX_BUNDLE_COPYRIGHT "Sandia National Labs, Will Johnson" )
SET(MACOSX_BUNDLE_GUI_IDENTIFIER "gov.sandia.macOS.InterSpec" )
#set MACOS_BUNDLE_VERSION_NUMBER for Info.plist.template to populate corresponding value in Xcode
set(MACOS_BUNDLE_VERSION_NUMBER 52)


SET( ${PRODUCT_NAME} "InterSpec" )
SET( MACOSX_BUNDLE_BUNDLE_NAME "InterSpec" )

FIND_LIBRARY(WEBKIT_LIBRARY WebKit)
MARK_AS_ADVANCED( WEBKIT_LIBRARY )

set(app_sources main.m AppDelegate.mm AppDelegate.h InterSpec.icns en.lproj/MainMenu.xib )


set( CMAKE_INSTALL_RPATH "@executable_path/../Frameworks" )
set( CMAKE_BUILD_RPATH "@executable_path/../Frameworks" )

# Xcode doesnt like having two bundles for some reason, so we will add the
# QuickLook project as an external project and add it to the final application
# at InterSpec.app/Contents/Library/QuickLook
include(ExternalProject)
ExternalProject_Add(SpecFilePreview
    BINARY_DIR   "${CMAKE_BINARY_DIR}/SpecFilePreview"
    URL ${CMAKE_CURRENT_SOURCE_DIR}/../macOsQuickLook/SpecFilePreview/SpecFilePreview/
    CMAKE_GENERATOR   "Unix Makefiles"
    CMAKE_CACHE_ARGS -DCMAKE_PREFIX_PATH:PATH=${CMAKE_PREFIX_PATH}
                     #-DWt_INCLUDE_DIR:PATH=${Wt_INCLUDE_DIR}
                     #-DBOOST_ROOT:PATH=${Boost_INCLUDE_DIR}/..
                     -DInterSpec_BASE_DIR:PATH=${CMAKE_CURRENT_SOURCE_DIR}/../../
                     #-DBOOST_INCLUDEDIR:PATH=${Boost_INCLUDE_DIR}
                     -DCMAKE_BUILD_TYPE:STRING=RELEASE
                     -DSpecUtils_ENABLE_D3_CHART:BOOL=OFF
                     INSTALL_COMMAND   ""
)

# Add dummy command to let cmake know SpecFilePreview.qlgenerator isnt available
# at configure time so it wont complain
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator"
    COMMAND ""
)


#In order to copy all the needed resources to the "Archive" builds of Xcode,
#  we need to specify each individual file, and later set their MACOSX_PACKAGE_LOCATION properties...
file(GLOB_RECURSE DATA_RES_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/../../data/*")
file(GLOB_RECURSE WT_RES_SOURCES "${WT_RESOURCES_DIRECTORY}/*")
file( GLOB_RECURSE INTERSPEC_RES_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources/*"
)

list( APPEND INTERSPEC_RES_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../../NOTICE.html )


file( GLOB_RECURSE EXAMPLE_RES_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../../example_spectra/*"
)


SET( EN_LOCALIZATION_FILES  "${CMAKE_CURRENT_SOURCE_DIR}/en.lproj/InfoPlist.strings"   
    "${CMAKE_CURRENT_SOURCE_DIR}/en.lproj/Credits.rtf"
)

find_program(IBTOOL ibtool HINTS "/usr/bin" "${OSX_DEVELOPER_ROOT}/usr/bin")
if(${IBTOOL} STREQUAL "IBTOOL-NOTFOUND")
    message(SEND_ERROR "could fine ibtool to compile the .xib files")
endif()

add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
    COMMAND
        ${IBTOOL} --errors --warnings --notices --output-format
            human-readable-text --compile
            "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
            "${CMAKE_CURRENT_SOURCE_DIR}/en.lproj/MainMenu.xib"
    COMMENT "Compiling MainMenu.xib"
)

ADD_EXECUTABLE( InterSpecApp MACOSX_BUNDLE ${app_sources}
    ${DATA_RES_SOURCES} ${WT_RES_SOURCES} ${INTERSPEC_RES_SOURCES} ${EXAMPLE_RES_SOURCES} ${EN_LOCALIZATION_FILES}
    "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
    "${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator"
)
    
SET_TARGET_PROPERTIES( InterSpecApp PROPERTIES MACOSX_BUNDLE TRUE )
#set_target_properties( InterSpecApp PROPERTIES SUFFIX "")
set_target_properties(InterSpecApp PROPERTIES OUTPUT_NAME "InterSpec")

target_link_libraries(InterSpecApp PUBLIC InterSpecLib ${WEBKIT_LIBRARY} )

SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_SOURCE_DIR}/InterSpec.icns 
    PROPERTIES MACOSX_PACKAGE_LOCATION Resources
)

add_dependencies(InterSpecApp SpecFilePreview)


if( USE_BATCH_TOOLS )
  add_subdirectory( ${CMAKE_CURRENT_SOURCE_DIR}/../batch InterSpec_batch )
  add_dependencies( InterSpecApp InterSpec_batch )
  
  # It seems like $<TARGET_BUNDLE_DIR:InterSpecLib> is a symlink to the actual build location, so we
  #  need to copy the symlinks path, but also, preserve symlinks inside the bundle.  So we'll use
  #  Unix cp to do this, after getting where the symlink points to
  add_custom_command(TARGET InterSpecApp POST_BUILD
    COMMAND mkdir -p $<TARGET_BUNDLE_DIR:InterSpecApp>/Contents/Frameworks && cp -a `readlink -f $<TARGET_BUNDLE_DIR:InterSpecLib>` $<TARGET_BUNDLE_DIR:InterSpecApp>/Contents/Frameworks/InterSpec.framework
  )

  # We want to look in InterSpec.app/Contents/Frameworks for InterSpecLib to link to
  set_target_properties( InterSpecApp PROPERTIES INSTALL_RPATH "@executable_path/../Frameworks" )
  set_target_properties( InterSpecApp PROPERTIES BUILD_WITH_INSTALL_RPATH 1 INSTALL_NAME_DIR "@rpath" )

  set_target_properties( InterSpec_batch PROPERTIES INSTALL_RPATH "@executable_path/../Frameworks" )
  set_target_properties( InterSpec_batch PROPERTIES BUILD_WITH_INSTALL_RPATH 1 INSTALL_NAME_DIR "@rpath" )

  # Copy the batch executable into the bundle - not really sure where is best, so we'll put it in Resources right now
  add_custom_command(TARGET InterSpecApp POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy
    $<TARGET_FILE:InterSpec_batch>
    $<TARGET_FILE_DIR:InterSpecApp>/
  )
endif( USE_BATCH_TOOLS )

#Set resources MACOSX_PACKAGE_LOCATION property so it is where we want in the app bundle
#ToDo: combine 5 loops below into a single loop...
# ToDo: Try to use: install(DIRECTORY InterSpec_resources DESTINATION $<TARGET_FILE_DIR:InterSpecApp>/../Resources/ )
# ToDo: Try using (something like):
#                 add_custom_command(TARGET ${PROJECT_NAME} PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory
#                                    "${CMAKE_SOURCE_DIR}/example_spectra" "$<TARGET_FILE_DIR:InterSpecApp>/Contents/Resources/")

# set( RES_SOURCES_DIRS DATA_RES_SOURCES WT_RES_SOURCES INTERSPEC_RES_SOURCES EXAMPLE_RES_SOURCES EN_LOCALIZATION_FILES )
# set( RES_DEST_DIRS data resources InterSpec_resources example_spectra en.lproj )

set_property( SOURCE ${CMAKE_CURRENT_BINARY_DIR}/SpecFilePreview/SpecFilePreview.qlgenerator
    PROPERTY MACOSX_PACKAGE_LOCATION "Library/QuickLook/"
)

foreach(RES_FILE ${DATA_RES_SOURCES})
    #Get the relative path from the data-folder to the particular file
    file( RELATIVE_PATH
        RES_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/../../data"
            ${RES_FILE}
    )
        
    get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
    #Set it's location inside the app package (under Resources)
    set_property( SOURCE ${RES_FILE}
        PROPERTY MACOSX_PACKAGE_LOCATION "Resources/data/${RES_PATH}"
    )
    source_group("Resources/data/${RES_PATH}" FILES "${RES_FILE}")
endforeach(RES_FILE)


foreach(RES_FILE ${WT_RES_SOURCES})
    file(RELATIVE_PATH RES_PATH "${WT_RESOURCES_DIRECTORY}" ${RES_FILE})
    get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
    set_property( SOURCE ${RES_FILE}
            PROPERTY MACOSX_PACKAGE_LOCATION "Resources/resources/${RES_PATH}"
    )
    source_group("Resources/resources/${RES_PATH}" FILES "${RES_FILE}")
endforeach(RES_FILE)


foreach(RES_FILE ${INTERSPEC_RES_SOURCES})
    file( RELATIVE_PATH RES_PATH
            "${CMAKE_CURRENT_SOURCE_DIR}/../../InterSpec_resources"
            ${RES_FILE}
    )
    get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
    set_property( SOURCE ${RES_FILE}
            PROPERTY MACOSX_PACKAGE_LOCATION "Resources/InterSpec_resources/${RES_PATH}"
    )
    source_group("Resources/InterSpec_resources/${RES_PATH}" FILES "${RES_FILE}")
endforeach(RES_FILE)


foreach(RES_FILE ${EXAMPLE_RES_SOURCES})
    file( RELATIVE_PATH 
        RES_PATH
        "${CMAKE_CURRENT_SOURCE_DIR}/../../example_spectra"
        ${RES_FILE}
    )

    get_filename_component(RES_PATH ${RES_PATH} DIRECTORY)
    set_property( SOURCE ${RES_FILE}
        PROPERTY MACOSX_PACKAGE_LOCATION "Resources/example_spectra/${RES_PATH}"
    )
    source_group("Resources/example_spectra/${RES_PATH}" FILES "${RES_FILE}")
endforeach(RES_FILE)


foreach(RES_FILE ${EN_LOCALIZATION_FILES})
    set_property( SOURCE ${RES_FILE}
        PROPERTY MACOSX_PACKAGE_LOCATION "Resources/en.lproj/"
    )
endforeach(RES_FILE)


set_property( SOURCE "${CMAKE_CURRENT_BINARY_DIR}/MainMenu.nib"
    PROPERTY MACOSX_PACKAGE_LOCATION "Resources/en.lproj/"
)


set_target_properties( InterSpecApp
    PROPERTIES MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/Info.plist.template"
)

set_target_properties( InterSpecApp
        PROPERTIES XCODE_ATTRIBUTE_CODE_SIGN_ENTITLEMENTS InterSpec.entitlements
)

# With out this next line, Xcode will (silently) fail to build Archive builds
set_target_properties( InterSpecApp
    PROPERTIES XCODE_ATTRIBUTE_INSTALL_PATH "/Applications"
)

#Figure out the development team ID, and set this info for Xcode so it will default to signing the app
#  This is probably a little brittle on other peoples machined
execute_process( COMMAND /usr/bin/security find-identity -v -p codesigning 
    OUTPUT_VARIABLE dev_team_ids
)

STRING(REGEX REPLACE ";" "\\\\;" dev_team_ids "${dev_team_ids}")
STRING(REGEX REPLACE "\n" ";" dev_team_ids "${dev_team_ids}")
LIST(FILTER dev_team_ids INCLUDE REGEX "Developer ID Application:.+" )
STRING(REGEX MATCH "\\(.+\\)" mac_team_id ${dev_team_ids} )
STRING(REGEX REPLACE "\\(|\\)" "" mac_team_id "${mac_team_id}")
message("Using Mac Developer Team ID='${mac_team_id}'")

set_target_properties( InterSpecApp PROPERTIES XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "Mac Developer" )
set_target_properties( InterSpecApp PROPERTIES XCODE_ATTRIBUTE_DEVELOPMENT_TEAM "${mac_team_id}" )
set_target_properties( InterSpecApp PROPERTIES XCODE_ATTRIBUTE_ENABLE_HARDENED_RUNTIME "YES" )

set_target_properties( InterSpecLib PROPERTIES XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "Mac Developer" )
set_target_properties( InterSpecLib PROPERTIES XCODE_ATTRIBUTE_DEVELOPMENT_TEAM "${mac_team_id}" )
set_target_properties( InterSpecLib PROPERTIES XCODE_ATTRIBUTE_ENABLE_HARDENED_RUNTIME "YES" )

if( USE_BATCH_TOOLS )
  set_target_properties( InterSpec_batch PROPERTIES XCODE_ATTRIBUTE_CODE_SIGN_IDENTITY "Mac Developer" )
  set_target_properties( InterSpec_batch PROPERTIES XCODE_ATTRIBUTE_DEVELOPMENT_TEAM "${mac_team_id}" )
  set_target_properties( InterSpec_batch PROPERTIES XCODE_ATTRIBUTE_ENABLE_HARDENED_RUNTIME "YES" )
endif( USE_BATCH_TOOLS )

