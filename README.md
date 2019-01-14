Code for the physics-based fine-tuning and relaxation steps.

## Compile

Compile this project using the standard cmake routine:

    mkdir build
    cd build
    cmake ..
    make

This should find and build the dependencies and create a `rods_bin` binary.

## Run

From within the root directory just issue:

    ./build/rods_bin

You will need to load in .rod files exported by the relax-field utility.
