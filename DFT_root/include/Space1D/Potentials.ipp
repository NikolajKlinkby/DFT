#include "Space1D/Potentials.hpp"
#include "slate/slate.hh"
#include "utils.hpp"
#include <blas.hh>
#include <mpi.h>
#include <map>

/* ---- Constructors and Destructors ---- */
template <typedef scalar_t>
Hamiltonian<scalar_t>::Hamiltonian(int size,int nb, int p, int q){
    //TODO fix this as a hermitian matrix!
    m_hamilton = slate::Matrix<scalar_t>(size, size, nb, p, q, slate::HostTask);

};

template <typedef scalar_t>
Hamiltonian<scalar_t>::Hamiltonian(int size,int nb, int p, int q, slate::Target target){
    //TODO fix this as a hermitian matrix!
    m_hamilton = slate::Matrix<scalar_t>(size, size, nb, p, q, target);

};

// Destructor
template <typedef scalar_t>
Hamiltonian<scalar_t>::Hamiltonian{
    m_hamilton.clear();
};

/* ---- Acces functions ---- */
template <typedef scalar_t>
slate::Matrix<scalar_t>* Hamiltonian<scalar_t>::GetPtr(){
    return *m_hamilton;
};

template <typedef scalar_t>
void Hamiltonian::SetExternal(void * external(const slate::Matrix<scaler_t> &, const Vector<scaler_t> &)){
    m_external_func = external;
};

template <typedef scalar_t>
void Hamiltonian::SetXC(const char *name){
    m_xc = functional_dictionary[name];
};

/* ---- Exchange-Correlations dictionary ---- */
template <typedef scalar_t>
std::map<const char *, void> functional_dictionary = {
        // LDA
        {"LDA",
         void LDA(slate::Matrix<scaler_t> density, Vector<scaler_t> lattice){

        }}
};

void ExchangePot(slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vx){
    for (int64_t i = 0; i < vx.mt(); i++){
        if (vx.tileIsLocal(i,0)){
            auto tile = vx(i, 0);
            auto tiledata = tile.data();
            for (int64_t ii; ii < tile.mb(); ii++){
                //index = i*(lattice.mt()+1)+ii
                tiledata[ii] = -pow(3./M_PI*vector_get(density,i*(vx.mt()+1)+ii),1./3.);
            }
        }
    }
}

double X(double rs, double b, double c){
    return rs + b * pow(rs, 1./2.) + c;
}

void CorrelationEnergy(slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc){
    double x0 = -0.10498;
    double b = 3.72744;
    double c = 12.9352;
    double Q = pow(4.*c-pow(b,2.),1./2.);

    for (int i = 0; i < density.m(); i++){
        scalar_t val_i = vector_get(density, i);
        if (val_i == 0){
            vector_set(vc, i, 0.);
        }
        else{
            vector_set(vc, i, val_i * (1.-log(2.))/pow(M_PI,2.) * (
                    -log(val_i * X(1./val_i,b,c)) + 2.*b/Q * atan(Q/(2./pow(val_i,1./2.) + b)) - b*x0/X(pow(x0,2.),b,c) * (
                            log(pow(pow(1./val_i,1./2.)-x0,2.)/X(1./val_i,b,c) + 2.*(2.*x0+b)/Q * atan(Q/(2./pow(val_i,1./2.) + b))  ))  ));
        }
    }
}

void CorrelationPot(slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc, slate::Matrix<scalar_t> buffer, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){
    CorrelationEnergy(density, buffer);

    interp_deriv(lattice, buffer, vc, x_mat, x_mat_deriv);

}

void LDAPot(slate::Matrix<scalar_t> lattice, slate::Matrix<scalar_t> density, slate::Matrix<scalar_t> vc, slate::Matrix<scalar_t> vx, slate::Matrix<scalar_t> vxc, slate::Matrix<scalar_t> buffer, slate::Matrix<scalar_t> x_mat, slate::Matrix<scalar_t> x_mat_deriv){
    ExchangePot(density, vx);
    CorrelationPot(lattice, density, vc, buffer, x_mat, x_mat_deriv);

    slate::copy(vxc, vx);
    slate::add(1.,vc,1.,vxc);
}

