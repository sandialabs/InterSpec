cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

project( SandiaDecayExample )

add_executable( sandia_decay_example sandia_decay_example.cpp )
target_link_libraries( sandia_decay_example SandiaDecay )


