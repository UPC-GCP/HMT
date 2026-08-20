# 3D Heat and Mass Transfer Model (WIP)
CFD Model capable of solving 1D, 2D and 3D cases. Sample files for the 1D Rod, Smith-Hutton, Lid Driven Cavity, Differentially Heated Cavity, and Square Cylinder cases located in ./ioSrc.

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
> "endTime": 50, "timeStep": 1e-3,


### Physical Data
1. **PHI0**: Initial value for Single Scalar map. (null = Navier-Stokes Solver, otherwise Single Scalar Solver)
2. **P0**: Initial value for the Pressure map. (null = Single Scalar Solver, otherwise Navier-Stokes Solver)
3. **T0**: Initial value for the Temperature map. (null = No temperature component in Navier-Stokes Solver)
4. **V0**: Initial value for the Velocity field. (Needs to be defined in all cases, treated as fixed for SSS, or as transport variables in NSS)
5. **g**: Gravity
5. **materials**: Registry of all materials with their corresponding properties. (Accepted Values: rho, gamma, cp, mu, beta)

> [!NOTE]
> "PHI0": null, "P0": 0, "T0": 0, "V0": [1, 0, 0], "g": null,\
> "materials": [\
>   &ensp;&ensp;&ensp;{"rho": 1, "gamma": 1, "cp": 1, "mu": 0.0008, "beta": 1}\
> ],


### Geometrical Data
1. **sections**: Definition of geometric regions with the corresponding material index and source term for each one.
2. **obstacles**: Definition of geometric regions blocked by an obstacle.

> [!NOTE]
> "sections": [\
>   &ensp;&ensp;&ensp;{"material": 0, "x0": [0, 0, 0], "x1": [2.5, 1.5, 1.5], "source": 0}\
> ],\
> "obstacles": [\
>   &ensp;&ensp;&ensp;{"x0": [0.4, 0.65, 0.65], "x1": [0.6, 0.85, 0.85]}\
> ],


### Mesh Data
1. **N**: Total amount of nodes for each axis.
2. **refinement**: Definition of mesh refinement with number of nodes and ranges defined for each refinement region. Needs to include the refinement algorithm and their corresponding parameters. (Accepted Values: Bidirectional [strength, centering], PowerLaw [kappa], HyperSingle [delta], HyperDouble [delta])

> [!NOTE]
> "N": [200, 40, 40],\
> "refinement": [\
>   &ensp;&ensp;&ensp;{"axis": 0, "N": 200, "range": [0, 2.5], "algorithm": "Bidirectional",  "strength": 0, "centering": 0.5},\
>   &ensp;&ensp;&ensp;{"axis": 1, "N": 40,  "range": [0, 1.5], "algorithm": "Unidirectional", "kappa": 1},\
>   &ensp;&ensp;&ensp;{"axis": 2, "N": 40,  "range": [0, 1.5], "algorithm": "HyperSingle",     "delta": 0.1}\
> ],


### Boundary Conditions
1. **boundariesPressure**: Boundary conditions for the Pressure map.
2. **boundariesVelocity**: Boundary conditions for the Velocity field.
3. **boundariesTemperature** : Boundary conditions for the Temperature map.

> [!NOTE]
> "boundariesPressure": [\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 0],   "x1": [0, 1.5, 1.5],   "side": 0, "value": 0},\
> &ensp;&ensp;&ensp;{"type": "Dirichlet", "x0": [2.5, 0, 0], "x1": [2.5, 1.5, 1.5], "side": 1, "value": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 0],   "x1": [2.5, 0, 1.5],   "side": 0, "value": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 1.5, 0], "x1": [2.5, 1.5, 1.5], "side": 1, "value": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 0],   "x1": [2.5, 1.5, 0],   "side": 0, "value": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 1.5], "x1": [2.5, 1.5, 1.5], "side": 1, "value": 0}\
> ],\
> "boundariesVelocity": [\
> &ensp;&ensp;&ensp;{"type": "Dirichlet", "x0": [0, 0, 0],   "x1": [0, 1.5, 1.5],   "side": 0, "uValue": 1, "vValue": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [2.5, 0, 0], "x1": [2.5, 1.5, 1.5], "side": 1, "uValue": 0, "vValue": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 0],   "x1": [2.5, 0, 1.5],   "side": 0, "uValue": 0, "vValue": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 1.5, 0], "x1": [2.5, 1.5, 1.5], "side": 1, "uValue": 0, "vValue": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 0],   "x1": [2.5, 1.5, 0],   "side": 0, "uValue": 0, "vValue": 0},\
> &ensp;&ensp;&ensp;{"type": "Neumann",   "x0": [0, 0, 1.5], "x1": [2.5, 1.5, 1.5], "side": 1, "uValue": 0, "vValue": 0}\
> ],


### Probe Data
1. **probes**: Definition of probe types for data logging, requires specifying the time interval, logging skips and geometric region for each probe. (Accepted Values: Map, Field, Debug)

> [!NOTE]
> "probes": [\
> &ensp;&ensp;&ensp;{"type": "Map",   "t": [0, 100], "x0": [0, 0, 0], "x1": [2.5, 1.5, 1.5], "nWrite": 100},\
> &ensp;&ensp;&ensp;{"type": "Field", "t": [0, 100], "x0": [0, 0, 0], "x1": [2.5, 1.5, 1.5], "nWrite": 100}\
> ],


### Medic Data
1. **medicOn**: Boolean that activates/deactivates the diagnostic tool.

> [!NOTE]
> "medicOn": false
