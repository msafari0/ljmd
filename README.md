This package contains simplified MD code with multi-threading
parallelization for simulating atoms with a Lennard-Jones potential.

The bundled makefiles are set up to compile the executable once
with OpenMP disabled and once with OpenMP enabled with each build
placing the various object files in separate directories.

The examples directory contains 3 sets of example input decks
and the reference directory the corresponding outputs.

Type: make
to compile everything and: make clean
to remove all compiled objects

## CMake implementation
We created three `CMakeLists.txt` file, mimicking the three correspondent `Makefile` respectively in the `02-ljmd/`, `02-ljmd/Obj-serial/`, and `02-ljmd/examples/` folder. In particular

+ The `CMakeLists.txt` file in `02-ljmd/` set the compiler, the subdirectories for the further cmake files and the make targets (`clean` and `check`)
+ The `CMakeLists.txt` file in `02-ljmd/Obj-serial/` build the object file and executable with the appropriate flags, set target for cleaning 
+ The `CMakeLists.txt` file in `02-ljmd/examples/` copies examples and reference files in the respective `build/` subdirectory, set target and commands for `check`. NB the check for the `argon_2916.dat` files are commented out as they fail (but they already fail in the `make` only version)

The commands

    mkdir build; cd build; cmake ..; cmake --build .

create the `build/` folder containing the `ljmd-serial.x` executable, which is a replica of the one obtained via the `Makefile`s. Once inside the folder, `cd build`, one can execute `make check` and `make clean`.