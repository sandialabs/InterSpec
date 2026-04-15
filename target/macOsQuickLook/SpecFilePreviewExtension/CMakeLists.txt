cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project( SpecFilePreviewExtensions )

# The InterSpec source root directory
set( InterSpec_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../.." )

# Set to 1 to use the old libharu/WPdfImage rendering path,
# or 0 to use the native Core Graphics CGPaintDevice path.
option( MACOS_QUICK_LOOK_USE_PDF "Use libharu PDF rendering for Quick Look (vs native CG)" OFF )

# --- Quick Look Preview Extension ---
# App extensions must be executables (MH_EXECUTE), not bundles (MH_BUNDLE).
add_executable( SpecFilePreview
  PreviewViewController.mm
  SpecPreviewCommon.mm
  SpecPreviewCommon.h
  CGPaintDevice.mm
  CGPaintDevice.h
)

set_target_properties( SpecFilePreview PROPERTIES
  MACOSX_BUNDLE TRUE
  MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/PreviewInfo.plist"
  XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "gov.sandia.macOS.InterSpec.SpecFilePreview"
  XCODE_ATTRIBUTE_INFOPLIST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/PreviewInfo.plist"
  XCODE_ATTRIBUTE_WRAPPER_EXTENSION "appex"
  XCODE_PRODUCT_TYPE "com.apple.product-type.app-extension"
  # The .appex is at Contents/PlugIns/X.appex/Contents/MacOS/
  # The framework is at Contents/Frameworks/
  INSTALL_RPATH "@executable_path/../../../../Frameworks"
  BUILD_WITH_INSTALL_RPATH TRUE
)

target_compile_definitions( SpecFilePreview PRIVATE
  MACOS_QUICK_LOOK_USE_PDF=$<BOOL:${MACOS_QUICK_LOOK_USE_PDF}>
)

target_link_libraries( SpecFilePreview PRIVATE
  InterSpecLib
  "-framework QuickLook"
  "-framework Cocoa"
  "-framework Quartz"
  "-framework CoreText"
  "-framework CoreGraphics"
  "-framework ImageIO"
  "-e _NSExtensionMain"
)


# --- Quick Look Thumbnail Extension ---
add_executable( SpecFileThumbnail
  ThumbnailProvider.mm
  SpecPreviewCommon.mm
  SpecPreviewCommon.h
  CGPaintDevice.mm
  CGPaintDevice.h
)

set_target_properties( SpecFileThumbnail PROPERTIES
  MACOSX_BUNDLE TRUE
  MACOSX_BUNDLE_INFO_PLIST "${CMAKE_CURRENT_SOURCE_DIR}/ThumbnailInfo.plist"
  XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "gov.sandia.macOS.InterSpec.SpecFileThumbnail"
  XCODE_ATTRIBUTE_INFOPLIST_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ThumbnailInfo.plist"
  XCODE_ATTRIBUTE_WRAPPER_EXTENSION "appex"
  XCODE_PRODUCT_TYPE "com.apple.product-type.app-extension"
  INSTALL_RPATH "@executable_path/../../../../Frameworks"
  BUILD_WITH_INSTALL_RPATH TRUE
)

target_compile_definitions( SpecFileThumbnail PRIVATE
  MACOS_QUICK_LOOK_USE_PDF=$<BOOL:${MACOS_QUICK_LOOK_USE_PDF}>
)

target_link_libraries( SpecFileThumbnail PRIVATE
  InterSpecLib
  "-framework QuickLook"
  "-framework QuickLookThumbnailing"
  "-framework Cocoa"
  "-framework CoreText"
  "-framework CoreGraphics"
  "-framework ImageIO"
  "-e _NSExtensionMain"
)

# Both targets need the InterSpec headers
target_include_directories( SpecFilePreview PRIVATE "${InterSpec_BASE_DIR}" )
target_include_directories( SpecFileThumbnail PRIVATE "${InterSpec_BASE_DIR}" )

# Copy logo into thumbnail extension's Resources for watermark.
# CMake generator expressions resolve to .app but Xcode builds .appex, so compute path directly.
set( LOGO_PNG "${InterSpec_BASE_DIR}/InterSpec_resources/images/InterSpec128.png" )
set( THUMBNAIL_APPEX_RESOURCES "${CMAKE_CURRENT_BINARY_DIR}/$<CONFIG>/SpecFileThumbnail.appex/Contents/Resources" )
add_custom_command( TARGET SpecFileThumbnail POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E make_directory "${THUMBNAIL_APPEX_RESOURCES}"
  COMMAND ${CMAKE_COMMAND} -E copy_if_different "${LOGO_PNG}" "${THUMBNAIL_APPEX_RESOURCES}/InterSpec128.png"
)

# App extensions must be sandboxed
set_target_properties( SpecFilePreview PROPERTIES
  XCODE_ATTRIBUTE_CODE_SIGN_ENTITLEMENTS "${CMAKE_CURRENT_SOURCE_DIR}/SpecFilePreview.entitlements"
)
set_target_properties( SpecFileThumbnail PROPERTIES
  XCODE_ATTRIBUTE_CODE_SIGN_ENTITLEMENTS "${CMAKE_CURRENT_SOURCE_DIR}/SpecFilePreview.entitlements"
)
