# IB3D-GPU
The solver is based on an Immersed Boundary method capable of simulating flow over fixed and arbitrarily moving non-grid conformal solid bodies.

Please cite the following article if this code helps your project:
Raj, A., Khan, P. M., Alam, M. I., Prakash, A., & Roy, S. (2023). A GPU-accelerated sharp interface immersed boundary method for versatile geometries. Journal of Computational Physics, 478, 111985.
https://doi.org/10.1016/j.jcp.2023.111985

# Compilation
nvfortran -Mr8 -acc -ta=tesla:cc80 *.f90 -o AppName

# Requires
1. nvidia hpc sdk-
https://developer.nvidia.com/hpc-sdk-downloads

# Case setup
The case set up is that of the flow over a sphere at Reynolds number 100.

# Output
The program output can be found in out directory.
