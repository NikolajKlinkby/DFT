#include "slate/slate.hh"
#include <blas.hh>
#include <mpi.h>

#pragma once
template <typename scalar_t>
class Vector{
public:
    Vector(int m,int nb, int p);
    Vector(int m,int nb, int p, slate::Target target);

    void Set(int index, scalar_t value);

    scalar_t Get(int index);

    void Multiply(slate::Matrix<scalar_t> matrix, scalar_t a = 1, scalar_t b = 1, slate::Option options = {
            // Set execution target to cpu Devices
            { slate :: Option :: Target , slate :: Target :: HostTask }});

private:
    slate::Matrix<scalar_t> vector_ptr;
};

#include "utils.ipp"

