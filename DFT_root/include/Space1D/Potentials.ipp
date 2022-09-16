#include "Space1D/Potentials.hpp"
#include "slate/slate.hh"
#include "utils.hpp"
#include <blas.hh>
#include <mpi.h>
#include <map>

/* -------- Potential funciton --------*/
template <typename scalar_t>
DualVector<scalar_t> HartreePot(const DualVector<scalar_t> & density){
    DualVector<scalar_t> hartree = density.EmptyLike();

    // TODO Make more efficient (code some multithreading!!!)

    // Hartree potential    v(r)=int dr' n(r') / |r' - r|
    //                      n(r) = n_up(r) + n_down(r)
    // Numerical            v(r) = Delta r sum_{r' neq r} n(r') / |r' - r|
    //                      v(i) = sum_{i' neq i} n(i') / |i' - i|
    //

    scalar_t pot;
    int64_t index;
    int64_t index_n;
    int64_t mt_h = hartree.GetUp().mt();
    int64_t mt_n = density.GetUp().mt();
    // Loop over tiles in hartree vector
    for (int64_t j = 0; j < mt_h; ++j){
        // Make sure that the tiles are local
        if(hartree.GetUp().tileIsLocal(j,0) &&
           hartree.GetDown().tileIsLocal(j,0)){
            // Acces the tile data
            auto tile_up = hartree.GetUp()(j,0);
            auto tile_down = hartree.GetDown()(j,0);
            auto tile_up_data = tile_up.data();
            auto tile_down_data = tile_down.data();

            //Loop over the tiles
            for (int64_t jj = 0; jj < tile_up.mb(); ++jj){
                index = j * (mt_h + 1) + jj; pot = 0;

                // Calculate v(index)
                for (int64_t i = 0; i < mt_n; ++i){
                    // Make sure that the tiles are local
                    if(density.GetUp().tileIsLocal(i,0) &&
                       density.GetDown().tileIsLocal(i,0)) {
                        // Acces the tile data
                        auto tile_up_n = density.GetUp()(i,0);
                        auto tile_down_n = density.GetDown()(i,0);
                        auto tile_up_data_n = tile_up.data();
                        auto tile_down_data_n = tile_down.data();

                        // Loop over tile in density
                        for (int64_t ii = 0; ii < tile_up_n.mb(); ++ii){
                            index_n = i* (mt_n + 1) + ii;

                            // Add to the potential v(index)
                            if (index > index_n){
                                pot += (tile_up_data[ii] + tile_down_data[ii]) / (index - index_n);
                            }
                            else if (index < index_n){
                                pot += (tile_up_data[ii] + tile_down_data[ii]) / (index_n - index);
                            }
                        }
                    }
                }

                // Set the potential
                tile_up_data[jj] = pot; tile_down_data[jj] = pot;
            }

        }
    }

    //Return
    return hartree;
}

template <typedef scalar_t>
DualVector<scalar_t> ExchangePot(const DualVector<scalar_t> & density){
    // Copying structure
    DualVector<scalar_t> exchange = density.EmptyLike();

    // Loop over tiles
    int64_t mt = exchange.GetUp().mt();
    for (int64_t i = 0; i < mt; i++){
        // Make sure tile is local
        if (exchange.GetUp().tileIsLocal(i,0) &&
            density.GetUp().tileIsLocal(i,0)){
            // Access tiles
            auto tile_up = exchange.GetUp()(i, 0);
            auto tile_down = exchange.GetDown()(i, 0);
            auto tiledata_up = tile_up.data();
            auto tiledata_down = tile_down.data();

            // Access density tiles
            auto tile_up_n = exchange.GetUp()(i, 0);
            auto tile_down_n = exchange.GetDown()(i, 0);
            auto tiledata_up_n = tile_up.data();
            auto tiledata_down_n = tile_down.data();

            // Loop over tile
            for (int64_t ii; ii < tile_up.mb(); ii++){
                tiledata_up[ii] = -pow( 6. / M_PI * tile_up_n[ii], 1./3.);
                tiledata_down[ii] = -pow( 6. / M_PI * tile_down_n[ii], 1./3.);
            }
        }
    }
}

/* ---- Local Spin Density Approximation (LSDA) ---- */
template <typedef scalar_t>
scalar_t X(scalar_t rs, scalar_t b, scalar_t c){
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


template <typedef scalar_t>
DualVector<scalar_t> LSDA(const DualVector<scalar_t> & density, const Vector<scaler_t> & lattice){
    // Set the exchange potential
    DualVector<scale_t> lsda = ExchangePot(density);

    // Calculate the correlation term

    // Add the correlation term

}

/* ---- Exchange-Correlations dictionary ---- */
template <typedef scalar_t>
using func = void(*)(const slate::Matrix<scalar_t> &, const Vector<scaler_t> &);

template <typedef scalar_t>
std::map<const char *, func<scalar_t>> functional_dictionary = {
        // LSDA
        {"LSDA"     ,    LSDA},
        {"LDA"      ,    LSDA}
};

/* -------- Potential -------- */
/* ---- Constructors and Destructors ---- */
template <typedef scalar_t>
Potential<scalar_t>::Potential(DualVector<scalar_t>* density, Vector<scalar_t>* lattice){
    // Setting pointers
    m_density = density;
    m_lattice = lattice;

    // Making potential dual vector
    m_pot = m_density->EmptyLike();
};

template <typedef scalar_t>
Potential<scalar_t>::~Potential(){};

template <typedef scalar_t>
DualVector<scalar_t>* Potential<scalar_t>::GetPtr(){
    return m_pot;
};

template <typedef scalar_t>
void Potential<scalar_t>::SetExternal(DualVector<scalar_t>* external(const DualVector<scalar_t> &, const Vector<scaler_t> &)){
    m_external_func = external;
    m_external_set = true;
};

template <typedef scalar_t>
void Potential<scalar_t>::Update(){
    // Set the potential to be equal to the hartree potential
    m_pot = HartreePot(m_density);

    // Add the Exchange Correlation potential
    m_pot.Add(m_xc(m_density, m_lattice));

    // Add the external potential (if set)
    if (m_external_set){
        m_pot.Add(m_external_func(m_density, m_lattice))
    }
}

template <typedef scalar_t>
void Potential<scalar_t>::SetXC(const char *name){
    //TODO check availability and spit out exception
    m_xc = functional_dictionary[name];
};



