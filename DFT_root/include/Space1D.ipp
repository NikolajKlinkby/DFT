#include "utils.hpp"
#include "Space1D.hpp"
#include <slate/slate.hh>
#include <blas.hh>
#include <mpi.h>

/* ---- settings ---- */

template <typedef scalar_t>
settings_struct<scalar_t>::set_lattice_points(int lattice_points){
    n_lattice_points = lattice_points;
};
template <typedef scalar_t>
settings_struct<scalar_t>::set_lattice_spacing(scalar_t lattice_spacing){
    m_lattice_spacing = lattice_spacing;
};
template <typedef scalar_t>
settings_struct<scalar_t>::set_nr_particles(int nr_particles){
    m_nr_particles = nr_particles;
    if (nr_particles % 2 == 0){
        m_nr_spin_up = nr_particles/2;
        m_nr_spin_down = m_nr_spin_up;
    }
    else{
        m_nr_spin_up = (nr_particles-1)/2;
        m_nr_spin_down = nr_spin_up+1;
    }
};
template <typedef scalar_t>
settings_struct<scalar_t>::set_nr_particles(int nr_spin_up, int nr_spin_down){
    m_nr_particles = nr_spin_up+nr_spin_down;
    m_nr_spin_up = nr_spin_up;
    m_nr_spin_down = nr_spin_down;
};
template <typedef scalar_t>
settings_struct<scalar_t>::set_scf_thres(scalar_t scf_thres){
    m_scf_thres = scf_thres;
};
template <typedef scalar_t>
settings_struct<scalar_t>::set_max_iter(int max_iter){
    m_max_iter = max_iter;
};

// Getters
template <typedef scalar_t>
int settings_struct<scalar_t>::get_lattice_points(){
    return m_lattice_points;
};
template <typedef scalar_t>
scalar_t settings_struct<scalar_t>::get_lattice_spacing(){
    return m_lattice_spacing;
};
template <typedef scalar_t>
int settings_struct<scalar_t>::get_nr_particles(){
    return m_nr_particles;
};
template <typedef scalar_t>
int settings_struct<scalar_t>::get_nr_spin_down(){
    return m_nr_spin_down;
};
template <typedef scalar_t>
int settings_struct<scalar_t>::get_nr_spin_up(){
    return m_nr_spin_up;
};
template <typedef scalar_t>
scalar_t settings_struct<scalar_t>::get_scf_thres(){
    return m_scf_thres;
};
template <typedef scalar_t>
int settings_struct<scalar_t>::get_max_iter(){
    return m_max_iter;
};



/* ---- Constructor/Destructor --- */

template <typedef scalar_t>
Space1D<scalar_t>::Space1D(int nr_particles){
    settings.set_nr_particles(nr_particles);
}

template <typedef scalar_t>
Space1D<scalar_t>::Space1D(int lattice, scalar_t spacing){
    settings.set_lattice_points(lattice);
    settings.set_lattice_spacing(spacing);
}

template <typedef scalar_t>
Space1D<scalar_t>::Space1D(int nr_particles, int lattice, scalar_t spacing){
    Space1D(nr_particles);
    Space1D(lattice, spacing);
}

template <typedef scalar_t>
Space1D<scalar_t>::Space1D(int nr_particles, int lattice, scalar_t spacing, int max_iter, scalar_t scf_thres){
    Space1D(nr_particles, lattice, spacing);
    settings.set_scf_thres(scf_thres);
    settings.set_max_iter(max_iter);
};

template <typedef scalar_t>
Space1D<scalar_t>::~Space1D(){};