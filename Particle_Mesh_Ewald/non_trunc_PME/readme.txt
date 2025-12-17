All methods presented above are Particle Mesh Ewald (PME) implementations with truncation errors strictly eliminated. The correctness of the non-truncated PME method can be verified by the following test: when all point charges in the system are precisely located on the grid points, the energy computed by the non-trunc PME algorithm must match the result from the direct Ewald summation exactly, regardless of the chosen grid size. For detailed theory, please refer to my doctoral dissertation.

Three implementation versions are provided:

non_truncPME.f90 (Fortran version)

non_truncPME.cpp (C++ version)

non_truncPME.cu (CUDA version for GPU acceleration)

All three versions satisfy this criterion and produce numerically identical results.
