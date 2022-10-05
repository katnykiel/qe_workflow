# qe_workflow

## Kat Nykiel, Bobby Appleton, Saswat Mishra

High-throughput density-functional theory calculations of high-pressure equations of state

### TODO:

Kat
- [ ] k-point / ENCUT convergence tests
- [X] script w/ functions
  - electronic + ionic energy convergence graph
  - mp api key
  - view structure objects
  - extract outputs from .stdout files (pymatgen monty regex?)
- [X] add notebook instructions

Bobby
- [ ] investigate change in enforced pressure
- [ ] loop for strain / ionic relaxation (sample enough for equation of state)
- [ ] fitting equations of state

Saswat
- [ ] convex hull calcs?
  - [X] list of what you need from us
    - [ ] Name of compound ($A_xB$)
    - [ ] x-axis as fraction of variable data ($x$)
    - [ ] y-axis as Energy of the compound
- [ ] collect pure states for one example?

General
- [ ] define chemical space
- [X] add slack integration
