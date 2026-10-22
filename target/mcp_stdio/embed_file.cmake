# Writes INPUT_FILE into OUTPUT_FILE as a C++ byte array named SYMBOL, with SYMBOL_size its length.
file( READ "${INPUT_FILE}" contents HEX )
string( LENGTH "${contents}" num_hex )
math( EXPR num_bytes "${num_hex} / 2" )
string( REGEX REPLACE "([0-9a-f][0-9a-f])" "0x\\1," bytes "${contents}" )
string( REGEX REPLACE "(0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,0x..,)" "\\1\n" bytes "${bytes}" )
file( WRITE "${OUTPUT_FILE}"
  "// Generated from ${INPUT_FILE}; do not edit.\n#include <cstddef>\n\n"
  "extern const unsigned char ${SYMBOL}[];\nextern const size_t ${SYMBOL}_size;\n\n"
  "const unsigned char ${SYMBOL}[] = {\n${bytes}\n};\nconst size_t ${SYMBOL}_size = ${num_bytes};\n" )
