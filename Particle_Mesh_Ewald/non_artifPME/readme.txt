# Non-artif PME

## Reference
PME method without discrete charge artifact errors. For details, see the **PhD dissertation of Yihao Zhao, Shandong University**.

## Implementation
**Fortran version**: `non_artifPME.f90`

- Based on the non-neutral PME framework, applicable to both charge-neutral and non-neutral systems
- Requires pre-calculation of the energy and force of a unit point charge at arbitrary positions within a fine grid
- Subsequently subtracts the discrete charge artifact contribution for each point charge in the system

## Note
This method **does not remove the net residual force** in the energy term of the original PME, and **the net residual force remains** even after subtracting the artifact force contributions. The reasons for this are explained in the relevant section of the doctoral dissertation.
