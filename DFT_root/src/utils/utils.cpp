#include "utils.h"
#include <slate/slate.hh>
#include <blas.hh>
#include <mpi.h>

#ifdef USE_EXPORT_KEYWORD
export
#endif

template <typename scalar_t>
Vector<scalar_t>::Vector(int m, int nb, int p) {
    vector_ptr = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    vector_ptr.insertLocalTiles(slate::Target::HostTask);
};

Vector<scalar_t>::Vector(int m, int nb, int p, slate::Target target) {
    vector_ptr = slate::Matrix < scalar_t >(m, 1, nb, p, 1, MPI_COMM_WORLD);
    vector_ptr.insertLocalTiles(target);
};

template <typename scalar_t>
void Vector<scalar_t>::Set(int index, scalar_t value){
    int mt = vector_ptr.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_ptr.tileMb(i) ){
            if (vector_ptr.tileIsLocal( i, 0)) {
                auto tile = vector_ptr( i, 0);
                auto tiledata = vector_ptr.data();
                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
};

template <typename scalar_t>
scalar_t Vector<scalar_t>::Get(int index){
    int mt = vector_ptr.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector_ptr.tileMb(i) ){
            if (vector_ptr.tileIsLocal( i, 0)) {
                auto tile = vector_ptr( i, 0 );
                auto tiledata = vector_ptr.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
};

template <typename scalar_t>
void Vector<scalar_t>::Multiply(slate::Matrix<scalar_t> matrix, scalar_t a, scalar_t b, slate::Option options){
    slate::multiply(matrix, a, b,  vector_ptr, vector_ptr, options);
};