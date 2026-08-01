# 2D Conduction Heat and Mass Transfer Model

## Config File Structure (.json)
Summary of the current capabilities of the model and the possible configuration parameters.

### Numerical Data
1. **tempScheme**: Temporal interpolation scheme. (Accepted Values: explicit, implicit, crank-nicolson)
2. **spatScheme**: Convective interpolation scheme. (Accepted Values: CDS, UDS, Hybrid, PowerLaw, Exponential)
3. **solver**: Numerical solver algorithm. (Accepted Values: CG, BiCG)
4. **tolTemporal**: Tolerance for the steady-state convergence check.
5. **tolNumeric**: Tolerance for the numerical solver. 
6. **maxIterations**: Limit of iterations per time-step for the numerical solver.
7. **endTime**: Total duration of the simulation.
8. **timeStep**: Time interval between instants.

> [!NOTE]
> "tempScheme": "implicit", "spatScheme": "Hybrid", "solver": "CG",\
> "tolTemporal": 1e-20, "tolNumeric": 1e-10, "maxIterations": 3000,\
> "endTime": 50, "timeStep": 1e-3,\

### Physical Data
1. **P0**: Initial value for the Pressure map.
2. **T0**: Initial value for the Temperature map.
3. **V0**: Initial value for the Velocity field.
4. **materials**: Registry of all materials with their corresponding properties. (Accepted Values: rho, gamma, mu)

> [!NOTE]
> "P0": 0, "T0": 0, "V0": [1, 0, 0],\
> "materials": [\
>   {"rho": 1, "gamma": 1, "mu": 0.0008}\
> ],

### Geometrical Data
1. **sections**: Definition of geometric regions with the corresponding material index and source term for each one.
2. **obstacles**: Definition of geometric regions blocked by an obstacle.

> [!NOTE]
> "sections": [\
>   {"material": 0, "x0": [0, 0, 0], "x1": [2.5, 1.5, 1.5], "source": 0}\
> ],\
> "obstacles": [\
>   {"x0": [0.4, 0.65, 0.65], "x1": [0.6, 0.85, 0.85]}\
> ],

### Mesh Data
1. **N**: Total amount of nodes for each axis.
2. **refinement**: Definition of mesh refinement with number of nodes and ranges defined for each refinement region. Needs to include the refinement algorithm and their corresponding parameters. (Accepted Values: Bidirectional [strength, centering], Unidirectional [kappa], HyperSingle [delta], HyperDouble [delta])

> [!NOTE]
> "N": [200, 40, 40],\
> "refinement": [\
>   {"axis": 0, "N": 200, "range": [0, 2.5], "algorithm": "Bidirectional",  "strength": 0, "centering": 0.5},\
>   {"axis": 1, "N": 40,  "range": [0, 1.5], "algorithm": "Unidirectional", "kappa": 1},\
>   {"axis": 2, "N": 40,  "range": [0, 1.5], "algorithm": "HyperSingle",     "delta": 0.1}\
> ],

### Boundary Conditions
1. **boundariesPressure**: Boundary conditions for the Pressure map.
2. **boundariesVelocity**: Boundary conditions for the Velocity field.
3. **boundariesTemperature** : Boundary conditions for the Temperature map.

### Probe Data
1. **probes**: Definition of probe types for data logging, requires specifying the time interval and geometric region for each probe. (Accepted Values: Map, Field, Debug)

### Medic Data
1. **medicOn**: Boolean that activates/deactivates the diagnostic tool.
