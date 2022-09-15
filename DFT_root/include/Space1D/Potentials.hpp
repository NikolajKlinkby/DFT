#include "slate/slate.hh"
#include "utils.hpp"
#include <blas.hh>
#include <mpi.h>

#pragma once

template <typedef scaler_t>
class Hamiltonian {
public:
    // Constructor
    Hamiltonian(int size,int nb, int p, int q);
    Hamiltonian(int size,int nb, int p, int q, slate::Target target);

    //Destructor
    Hamiltonian{};

    // Pointer to the matrix
    slate::Matrix<scaler_t>* GetPtr();

    // Update the hamilton with current potentials due to specific density and lattice
    void Update(const slate::Matrix<scaler_t> & density, const Vector<scaler_t> & lattice);

    // Set external potential
    void SetExternal(void * external(const slate::Matrix<scaler_t> &, const Vector<scaler_t> &));

    // Set exchange correlation potential
    void SetXC(const char * name);

private:
    // The matrix
    slate::Matrix<scaler_t> m_hamilton;

    // Functions
    void                    (*m_external_func)(const slate::Matrix<scaler_t> &, const Vector<scaler_t> &); // External potential to be set by user
    void                    (*m_xc)(const slate::Matrix<scaler_t> &, const Vector<scaler_t> &); // Exchange-correlation potential set (LDA by default)
};


#include "Space1D/Potentials.ipp"