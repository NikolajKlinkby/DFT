#include "slate/slate.hh"
#include "utils.hpp"
#include <blas.hh>
#include <mpi.h>

#pragma once

template<typedef scalar_t>
class Potential{
public:
    // Constructor
    Potential(DualVector<scalar_t>* density, Vector<scalar_t>* lattice);

    // Destructor
    ~Potential();

    // Pointer to potential
    DualVector<scalar_t>* GetPtr();

    // Set external function
    void SetExternal(DualVector<scalar_t>* external(const DualVector &, const Vector<scaler_t> &));

    // Update
    void Update();

    // Set exchange correlation potential
    void SetXC(const char * name);

private:
    // Potential
    DualVector<scalar_t>    m_pot;

    // Density and lattice pointers
    DualVector<scalar_t>*   m_density;
    Vector<scalar_t>*       m_lattice;

    // Functions
    DualVector<scalar_t>     (*m_external_func)(const DualVector &, const Vector<scaler_t> &); // External potential to be set by user
    DualVector<scalar_t>     (*m_xc)(const slate::Matrix<scaler_t> &, const Vector<scaler_t> &); // Exchange-correlation potential set (LDA by default)
    bool                      m_external_set=false;
};


#include "Space1D/Potentials.ipp"