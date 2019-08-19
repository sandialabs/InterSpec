#adapted from https://github.com/Kitware/CMake/tree/master/Tests/CFBundleTest
cmake_minimum_required(VERSION 3.1.0)

#project( SpecFilePreview )


set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

set( SOURCES
    main.c
    GenerateThumbnailForURL.m
    GeneratePreviewForURL.m
    	
    QLSpectrumChart.h
    QLSpectrumChart.cpp
    SpecPreviewCommon.h
    SpecPreviewCommon.cpp
      
    QLPeakDef.h
    QLPeakDef.cpp
    QLSpecMeas.h
    QLSpecMeas.cpp

    QLSpectrumDataModel.h
    QLSpectrumDataModel.cpp
    
    Info.plist
)

add_library( SpecFilePreview MODULE ${SOURCES} )


set(sysroot)
if(CMAKE_OSX_SYSROOT)
  set(sysroot -isysroot ${CMAKE_OSX_SYSROOT})
endif()



# note that for some reason, the makefile and xcode generators use a different
# property to indicate where the Info.plist file is :-/ For that reason, we
# specify it twice so it will work both places
set_target_properties( SpecFilePreview PROPERTIES
    BUNDLE 1
    BUNDLE_EXTENSION qlgenerator
    XCODE_ATTRIBUTE_MACH_O_TYPE mh_bundle
    XCODE_ATTRIBUTE_INFOPLIST_FILE ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist
    MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist
    )
    #LINK_FLAGS "-Wl,-exported_symbols_list,\"${CMAKE_CURRENT_SOURCE_DIR}/ExportList_plugin.txt\"")

#export(TARGETS SpecFilePreview FILE ${CMAKE_CURRENT_BINARY_DIR}/exp.cmake)

target_link_libraries( SpecFilePreview PUBLIC 
  SpecUtils
  ${HPDF_LIBRARY}
  ${PNG_LIBRARY}
  ${LIBZ_LIBRARY}
  ${Wt_LIBRARIES}
  ${Boost_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_REGEX_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_RANDOM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${Boost_DATE_TIME_LIBRARY}
  /System/Library/Frameworks/QuickLook.framework/QuickLook
  /System/Library/Frameworks/CoreServices.framework/CoreServices
)

target_include_directories( SpecFilePreview PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../../../external_libs )
message( "SpecFilePreview: ${CMAKE_CURRENT_SOURCE_DIR}/../../../../external_libs" )


#install(TARGETS SpecFilePreview DESTINATION InterSpec.app/Contents/Library/QuickLook/ )
