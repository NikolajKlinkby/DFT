#include "Space1D.hpp"
#include "utils.hpp"
#include "slate/slate.hh"
#include <blas.hh>
#include <mpi.h>


int main(int argc , char ** argv){

    //Initialise MPI
    int err = 0 , mpi_provided = 0;
    err = MPI_Init_thread ( &argc , &argv , MPI_THREAD_MULTIPLE , & mpi_provided);
    assert ( err == 0 && mpi_provided == MPI_THREAD_MULTIPLE );


    int m = 100;
    Vector<double>* vec = new Vector<double>(m, 1, 1);

    //Collect MPI error
    err = MPI_Finalize ();
    assert ( err == 0 );

    return 0;
}