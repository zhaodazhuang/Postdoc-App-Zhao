All methods presented above are Particle Mesh Ewald (PME) implementations for computable non-electroneutral systems after eliminating truncation errors. For detailed theoretical foundations, please refer to my doctoral dissertation.

The correctness can be verified via the dual-charge system energy relation:  
E(q₁, q₂) = –E(–q₁, q₂) = –E(q₁, –q₂).

Two versions are provided:  
- `non_neutralPME.cpp` (C++ version)  
- `non_neutralPME.f90` (Fortran version)

Both implementations produce identical numerical results.
