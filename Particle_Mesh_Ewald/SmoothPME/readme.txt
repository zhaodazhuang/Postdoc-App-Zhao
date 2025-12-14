# Smooth Particle Mesh Ewald (SPME)

## Reference
"A smooth particle mesh Ewald method"  
*J. Chem. Phys.* **103**, 8577 (1995)  
https://doi.org/10.1063/1.470117

## Implementations
**C++ version**: `SPME.cpp`  
**Fortran version**: `SPME.f90`

- Both implementations are algorithmically identical and have been cross-validated
- Uses B-spline interpolation (analytically differentiable)
- Enforces zero net residual force to achieve momentum conservation

*Note: The original method guarantees energy conservation but not momentum conservation. This modification ensures momentum conservation at the expense of strict energy conservation.*
