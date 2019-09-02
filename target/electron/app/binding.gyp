
{
    "targets": [{
        "target_name": "InterSpecAddOn",
        "cflags!": [ "-fno-exceptions" ],
        "cflags_cc!": [ "-fno-exceptions" ],
        "sources": [
            "InterSpecAddOn.cpp"
        ],
        'include_dirs': [
            "<!@(node -p \"require('node-addon-api').include\")",
            "../../..",
            "../../../build_electron"
        ],
        'libraries': [
          "-Wl,-rpath,@loader_path",
          "-Wl,-rpath,/Users/wcjohns/rad_ana/InterSpec/build_electron/Release/",
          "../../../../build_electron/Release/libLibInterSpec.dylib"
          ],
        'dependencies': [
            "<!(node -p \"require('node-addon-api').gyp\")"
        ],
        'defines': [ 'NAPI_DISABLE_CPP_EXCEPTIONS' ]
    }]
}