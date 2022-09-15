## TODO List

Have in mind that this should all be multithreaded CPU and GPU wise. Maybe cross threading should be possible.

List of things
- build settings structure
  - include: #Lattice points, Spacing, #Particles, Threshold convergence for SCF, Maximum number of iterations for SCF
  - build setters and getter functions for the settings
- create lattice and density space function
  - 1D vector for lattice
  - 2D Matrix for density (1 vector for spin up and one for spin down)
- select and build potentials
  - potential classes: external, hartree, exchange-correlation (LDA to begin with) 
  - create hamiltonian from potentials
    - maybe the hamiltonian(the potentials) should be a class advancing with the density?
- create Kohn-Sham (KS) orbitals (Wavefunctions)
  - make initial guess of KS orbitals
- Self Consistent Field (SCF)
- calculate observables
  - Which?
