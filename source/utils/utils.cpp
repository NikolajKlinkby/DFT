#include "utils.h"
#include <slate/slate.hh>
#include <blas.hh>
#include <mpi.h>

template <typename scalar_t>
Vector<scalar_t>::Vector(int m, int nb, int p, slate::Target target) {
    Vector::vector_ptr = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    Vector::vector_ptr.insertLocalTiles(target);
};

template <typename scalar_t>
void Vector<scalar_t>::Set(int index, scalar_t value){
    int mt = Vector::vector_ptr.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+Vector::vector_ptr.tileMb(i) ){
            if (Vector::vector_ptr.tileIsLocal( i, 0)) {
                auto tile = Vector::vector_ptr( i, 0);
                auto tiledata = Vector::vector_ptr.data();
                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
};

template <typename scalar_t>
scalar_t Vector<scalar_t>::Get(int index){
    int mt = Vector::vector_ptr.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+Vector::vector_ptr.tileMb(i) ){
            if (Vector::vector_ptr.tileIsLocal( i, 0)) {
                auto tile = Vector::vector_ptr( i, 0 );
                auto tiledata = Vector::vector_ptr.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
};

template <typename scalar_t>
void Vector<scalar_t>::Multiply(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  Vector::vector_ptr, Vector::vector_ptr, options);
};