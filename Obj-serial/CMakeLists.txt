# CMake build for P2.1 code
cmake_minimum_required(VERSION 3.10)


# COMPILER
set(CMAKE_C_COMPILER gcc)
# GENERAL FLAGS
set(CMAKE_C_FLAGS "-Wall -g -O3 -ffast-math -fomit-frame-pointer -DLJMD_VERSION=1.0")

# create object file
add_library(ljmd.o OBJECT ../src/ljmd.c)
# Link math library
target_link_libraries(ljmd.o m)

# Add executable 
add_executable(ljmd-serial.x $<TARGET_OBJECTS:ljmd.o>)
# Set output directory for executable
set_target_properties(ljmd-serial.x PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/..)

# Add the clean target (as in make clean)
add_custom_target(clean_obj COMMAND rm -f *.mod *.o ../ljmd-serial.x)
