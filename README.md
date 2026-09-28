# Research-Oriented Computational Kernel (ROCK) - Configuration File Structure
Sample files in ./ioSrc

## Configuration Data
1. **configSolver**: Type of solver to run. (FVMScalar, FVMNS, Spectral)


## FVMScalar
Finite Volume Method single scalar numerical solver.

### Numerical Data
1. **tempScheme**: Temporal interpolation scheme. (explicit, implicit, crank-nicolson)
2. **spatScheme**: Convective interpolation scheme. (CDS, UDS, Hybrid, PowerLaw, Exponential)
3. **solver**: Numerical solver algorithm. (Accepted Values: CG, BiCG)
4. **tolTemporal**: Tolerance for the steady-state convergence check.
5. **tolNumeric**: Tolerance for the numerical solver. 
6. **maxIterations**: Limit of iterations per time-step for the numerical solver.
7. **endTime**: Total duration of the simulation.
8. **timeStep**: Time interval between instants.

### Physical Data
1. **PHI0**: Initial value for Single Scalar map.
2. **V0**: Initial value for the Velocity field.
3. **g**: Gravity
4. **materials**: Registry of all materials with their corresponding properties. (rho, gamma, cp, mu, beta, alpha)
5. **sections**: Registry of geometric regions defining material index and source term.
6. **obstacles**: Registry of geometric regions blocked by an obstacle.

### Mesh Data
1. **N**: Total amount of nodes for each axis.
2. **refinement**: Definition of mesh refinement with number of nodes and ranges defined for each refinement region. Needs to include the refinement algorithm and their corresponding parameters. (Bidirectional [strength, centering], PowerLaw [kappa], HyperSingle [delta], HyperDouble [delta])

### Boundary Conditions
1. **boundariesPhi**: Boundary conditions for Scalar variable.

### Probe Data
1. **probes**: Definition of probe types for data logging, requires specifying the time interval, logging skips and geometric region for each probe. (Accepted Values: Map, Field, Debug) --- This may change significantly

### Medic Data
1. **medicOn**: Boolean to activate the diagnostic tool.


## FVMNS
Finite Volume Method Navier-Stokes numerical solver.

### Numerical Data
Same as FVMScalar.

### Physical Data
1. **P0**: Initial value for the Pressure map.
2. **T0**: Initial value for the Temperature map. (null == No temperature object initialized.)
3. **V0**: Initial value for the Velocity field.
4. **g**: Gravity
5. **materials**: Registry of all materials with their corresponding properties. (Accepted Values: rho, gamma, cp, mu, beta)
6. **sections**: Registry of geometric regions defining material index and source term.
7. **obstacles**: Registry of geometric regions blocked by an obstacle.

### Mesh Data
Same as FVMScalar.

### Boundary Conditions
1. **boundariesPressure**: Boundary conditions for Pressure map.
2. **boundariesVelocity**: Boundary conditions for Velocity field.
3. **boundariesTemperature** : Boundary conditions for Temperature map.

### Probe Data
Same as FVMScalar.

### Medic Data
Same as FVMScalar.


## Spectral
Spectral Methods numerical solver.

### Numerical Data
1. **N**: Truncated mode limit.
1. **solver**: Numerical solver selection. (DNS, LES)
2. **endTime**: Total duration of simulation.
3. **timeStep**: Time interval between instants.

### Physical Data
1. **Re**: Reynolds number.

### Probe Data
Same as FVMScalar.

### Medic Data
Same as FVMScalar.
