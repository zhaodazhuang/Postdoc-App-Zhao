# Midtown PME

Particle Mesh Ewald (PME) implementation using Midtown splines interpolation.

## Reference

Method based on: "Midtown splines for particle mesh Ewald method" [J. Chem. Phys. 153, 144115 (2020)](https://doi.org/10.1063/5.0021496)

## Implementations

### C++ Version
- File: `Midtown_PME.cpp`
- CPU implementation using Midtown splines
- Supported interpolation orders: 4 or 6 only
- Higher orders are more complex and not implemented

### CUDA Version  
- File: `Midtown_PME.cu`
- GPU acceleration with CUDA
- Uses cuFFT library and atomic operations
- Compilation command for RTX 4060:
  ```bash
  nvcc -arch=sm_89 Midtown_PME.cu -lcufft -o Midtown_PME   


Validation
Both versions have been verified using the same H2O.xyz configuration and produce identical results.
