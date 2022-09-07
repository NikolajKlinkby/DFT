#include "utils.h"
#include <slate/slate.hh>
#include <blas.hh>
#include <mpi.h>

template <typename scalar_t>
void vector_set(slate::Matrix<scalar_t> vector, int index, scalar_t value){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                std::cout << "hej 3" << std::endl;
                auto tile = vector( i, 0);
                std::cout << "hej 4" << std::endl;
                auto tiledata = tile.data();

                tiledata[index-i*(mt+1)] = value;
                return;
            }
        }
    }
}

template <typename scalar_t>
scalar_t vector_get(slate::Matrix<scalar_t> vector, int index){
    int mt = vector.mt();
    for (int64_t i = 0; i < mt; i++){
        if ( i*(mt+1) <= index && index < i*(mt+1)+vector.tileMb(i) ){
            if (vector.tileIsLocal( i, 0)) {
                auto tile = vector( i, 0 );
                auto tiledata = tile.data();
                return tiledata[index-i*(mt+1)];
            }
        }
    }
}

template <typename scalar_t>
scalar_t vector_get_max(slate::Matrix<scalar_t> vector){
    scalar_t max;
    for (int64_t i = 0; i < vector.mt(); i++){
        if(vector.tileIsLocal(i,0)){
            auto tile = vector(i,0);
            auto tiledata = tile.data();
            for (int64_t ii=0; ii < tile.mb(); ii++){
                if(tiledata[ii]>max){
                    max = tiledata[ii];
                }
            }
        }
    }
    return max;
}

template <typename scalar_t>
scalar_t vector_get_min(slate::Matrix<scalar_t> vector){
    scalar_t min;
    for (int64_t i = 0; i < vector.mt(); i++){
        if(vector.tileIsLocal(i,0)){
            auto tile = vector(i,0);
            auto tiledata = tile.data();
            for (int64_t ii=0; ii < tile.mb(); ii++){
                if(tiledata[ii]<min){
                    min = tiledata[ii];
                }
            }
        }
    }
    return min;
}

template <typename scalar_t>
void interp_deriv(slate::Matrix<scalar_t> x, slate::Matrix<scalar_t> y, slate::Matrix<scalar_t> y_deriv, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){

    /*x_mat and x_mat_deriv need to be equally created*/

    //Setting matrix
    int64_t mt = x_mat.mt();
    int64_t nt = x_mat.nt();
    int64_t n = x_mat.n();
    scalar_t x_val;

    for (int64_t i=0; i < mt; i++){
        for (int64_t j=0; j < nt; j++){



            if (x_mat.tileIsLocal(i,j) && x_mat_deriv.tileIsLocal(i,j)){



                auto tile_x = x_mat(i,j);
                auto tiledata_x = tile_x.data();
                auto tile_d = x_mat_deriv(i,j);
                auto tiledata_d = tile_d.data();

                int64_t mb = tile_x.mb();
                int64_t nb = tile_x.nb();
                int64_t stride = tile_x.stride();

                for (int64_t ii=0; ii < mb; ii++){
                    for (int64_t jj=0; jj < nb; jj++){
                        x_val = vector_get(x,i*(mt+1)+ii);
                        tiledata_x[ii + jj*stride] = pow(x_val , (n-1-j*(nt+1)-jj));
                        if (n-2-j*(nt+1)-jj < 0){
                            tiledata_d[ii + jj*stride] = 0;
                        }
                        else{
                            tiledata_d[ii + jj*stride] = tiledata_x[ii + jj*stride]/x_val;
                        }
                    }
                }
            }
        }
    }




    //Finding derivative
    slate::Pivots pivot;
    slate::gesv(x_mat,pivot,y); //Overwrites y as output

    slate::multiply(1.,x_mat_deriv,y,0.,y_deriv);

}