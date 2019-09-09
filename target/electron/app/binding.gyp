
{
    "targets": [{
        "target_name": "InterSpecAddOn",
        "cflags!": [ "-fno-exceptions" ],
        "cflags_cc!": [ "-fno-exceptions" ],
        'xcode_settings': {
          'GCC_ENABLE_CPP_EXCEPTIONS': 'YES',
          'CLANG_CXX_LIBRARY': 'libc++',
          'MACOSX_DEPLOYMENT_TARGET': '10.10',
        },
        'msvs_settings': {
          'VCCLCompilerTool': { 'ExceptionHandling': 1 },
        },
        "sources": [
            "InterSpecAddOn.cpp"
        ],
        'conditions': [
          ['OS=="mac"', {
            'cflags+': ['-fvisibility=hidden'],
            'xcode_settings': {
              'GCC_SYMBOLS_PRIVATE_EXTERN': 'YES'
            }
          }],
        ],
        'include_dirs': [
            "<!@(node -p \"require('node-addon-api').include\")",
            "../../..",
            "../../../build_electron"
        ],
        'libraries': [
          "-Wl,-rpath,@loader_path",
          "-Wl,-rpath,/Users/wcjohns/rad_ana/InterSpec/build_electron/Release/",
          "../../../../build_electron/Release/libLibInterSpec.a"
          ],
        'dependencies': [
            "<!(node -p \"require('node-addon-api').gyp\")"
        ],
        'defines': [ 'NAPI_DISABLE_CPP_EXCEPTIONS' ]
    }]
}