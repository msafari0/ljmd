# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)
project(ljmd VERSION 1.0 LANGUAGES C) 
set(CMAKE_BUILD_TYPE RelWithDebInfo Release) 


# COMPILER
set(CMAKE_C_COMPILER gcc)


# Add executable 
add_executable(ljmd-serial.x src/ljmd.c)
# compile options
target_compile_options(ljmd-serial.x PRIVATE -Wall -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0)
# include directory for header files
include_directories(include)

# add library for force computation
add_library(mdlib src/utilities.c src/verlet_ti.c src/forces.c)
target_link_libraries(mdlib m)

# link libraries
target_link_libraries(ljmd-serial.x m)
target_link_libraries(ljmd-serial.x mdlib)


# add subdirectories
add_subdirectory(examples)


# set target for serial version
add_custom_target(serial
    COMMAND make
)

# set target for cleaning
add_custom_target(clean_data_files
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples clean_example
)

# set target for checking
add_custom_target(check
    DEPENDS serial
    COMMAND ${CMAKE_MAKE_PROGRAM} -C examples check_example
)

