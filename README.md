# 2D Conduction Heat and Mass Transfer Model

## Config File Structure (.json)
example.json is also provided for reference.

### Numerical Data
Specify numerical parameters such as scheme, solver, tolerance, time-step and limits for end time and maximum iterations per time-step.

### Physical Data
Specify initial temperature and list all materials with their respective properties. Variable properties can be implemented in a latter patch.

### Geometrical Data
Specify sections of the analyzed geometry, including material and internal heat generation.

### Mesh Data
Specify mesh parameters, mesh algorithm, nodes per axis, and list refinement regions for the mesh.

### Boundary Conditions
Specify the boundary conditions present in the model. Dirichlet boundary conditions also support variable temperature values.

### Probe Data
Specify data acquisition targets.
