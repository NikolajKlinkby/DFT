#include "slate/slate.hh"
#include <blas.hh>
#include <mpi.h>

#pragma once

template <typedef scalar_t>
struct settings_struct{
public:
    settings_struct(){};

    // Setters
    set_lattice_points(int lattice_points);
    set_lattice_spacing(scalar_t lattice_spacing);
    set_nr_particles(int nr_particles);
    set_nr_particles(int nr_spin_up, int nr_spin_down);
    set_scf_thres(scalar_t scf_thres);
    set_max_iter(int max_iter);

    // Getters
    int get_lattice_points();
    scalar_t get_lattice_spacing();
    int get_nr_particles();
    int get_nr_spin_down();
    int get_nr_spin_up();
    scalar_t get_scf_thres();
    int get_max_iter();

private:
    int m_lattice_points = 100; // #Lattice points
    scalar_t m_lattice_spacing = 0.05; //Spacing
    int m_nr_particles = 2; // #Particles
    int m_nr_spin_down = 1;
    int m_nr_spin_up = 1;
    scalar_t m_scf_thres = 0.000001; //Threshold convergence for SCF
    int m_max_iter = 100; //Maximum number of iterations for SCF
};

template <typedef scalar_t>
class Space1D {
public:
    // Generate standard settings
    settings_struct<scalar_t> settings;

    // Constructor and destructor
    Space1D(){};
    Space1D(int nr_particles);
    Space1D(int lattice, scalar_t spacing);
    Space1D(int nr_particles, int lattice, scalar_t spacing);
    Space1D(int nr_particles, int lattice, scalar_t spacing, int max_iter, scalar_t scf_thres);
    ~Space1D();


    //TODO create Kohn-Sham (KS) orbitals (Wavefunctions)

        //TODO implicit assume equal amounts of spin up and spin down (for odd number of electrons assume opne extra spin down)

    //TODO make initial guess of KS orbitals

    //TODO create lattice and density space

    //TODO set potentials

        //TODO - create hamiltonian from potentials


    //TODO Self Consistent Field (SCF)

    //TODO calculate observables

};

#include "Space1D.ipp"