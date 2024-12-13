# IB3D-GPU
The solver is based on an Immersed Boundary method capable of simulating flow over fixed and arbitrarily moving non-grid conformal solid bodies.

# Compilation
nvfortran -Mr8 -acc -ta=tesla:cc80 *.f90 -o AppName

# Requires
1. nvidia hpc sdk-
https://developer.nvidia.com/hpc-sdk-downloads

# Case setup
The case set up is that of the flow over a sphere at Reynolds number 100.

# Output
The program output can be found in out directory.
