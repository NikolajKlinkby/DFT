## TODO List

Have in mind that this should all be multithreaded CPU and GPU wise. Maybe cross threading should be possible.

List of things
- installer
  - build objects for templates: float, double, long double
- build better utils.hpp/ipp
- build settings structure
  - multithreading support
- create Kohn-Sham (KS) orbitals (Wavefunctions)
  - class of orbitals
    - get density functions: spin up, spin down, total, polarization
  - make initial guess of KS orbitals
- select and build potentials
  - DEBUG!
  - potential classes: external 
  - Update function
- Self Consistent Field (SCF)
  - build hamilton class
  - while loop until convergence
- calculate observables
  - Which?
