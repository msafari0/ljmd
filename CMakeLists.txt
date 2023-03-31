# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C) 

# COMPILER
set(CMAKE_C_COMPILER gcc)

# add subdirectories
add_subdirectory(Obj-serial)
add_subdirectory(examples)

# set target for serial version
add_custom_target(serial
    COMMAND make
)

# set target for cleaning
add_custom_target(clean_all
    COMMAND ${CMAKE_MAKE_PROGRAM} -C Obj-serial clean_obj
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples clean_example
)

# set target for checking
add_custom_target(check
    DEPENDS serial
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples check_example
)

