#include "Space1D/Potentials.hpp"
#include "slate/slate.hh"
#include "utils.hpp"
#include <blas.hh>
#include <mpi.h>
#include <map>
#include <math.h>

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

/* ---- Local Spin Density Approximation (LSDA) ---- */
template <typedef scalar_t>
DualVector<scalar_t> LSDAExchangePot(const DualVector<scalar_t> & density){
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
            auto tile_up_n = density.GetUp()(i, 0);
            auto tile_down_n = density.GetDown()(i, 0);
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

// Emperical functions from VWN
template <typedef scalar_t>
scalar_t X(scalar_t rs, scalar_t b, scalar_t c){
    return rs + b * pow(rs, 1./2.) + c;
}

template <typedef scalar_t>
scalar_t f_z(scalar_t zeta){
    return (pow(1.+zeta,4./3.) + pow(1.-zeta,4./3.) - 2)/(pow(2.,4./3.)-2.);
}

template <typedef scalar_t>
scalar_t alpha_c(scalar_t rs){
    return -log(rs)/(6*pow(M_PI,2.))+0.03547;
}

// para- and ferro-magnetic correlation energies
template <typedef scalar_t>
scalar_t ec(scalar_t rs, scalar_t x0, scalar_t b, scalar_t c, scalar_t Q){
    scalar_t x = pow(rs,1./2.);
    scalar_t Xx = X(rs, b, c);
    scalar_t Qtan = atan(Q / (2 * pow(rs,1./2.) + b));
    return 3. / (4.*pow(M_PI*rs,3.)) * (1.-M_LN2) * (log(rs/Xx) + 2.*b/Q*Qtan - b*x0/X(pow(x0,2.)) * (log(pow(x-x0,2.)/Xx) + 2*(2*x0+b) / Q*qtan));
}

// Spin correlation energy
template <typedef scalar_t>
scalar_t ecs(scalar_t n_up, scalar_t n_down, const scalar_t & para[4], const scalar_t & ferro[4]){
    scalar_t rs = pow(3./(4*M_PI*(n_up + n_down)),1./3.);               //Wigner-Seitz radius
    scalar_t zeta = (n_up - n_down)/(n_up + n_down);                    //Polarization density
    scalar_t ec_para = ec(rs, para[0], para[1], para[2], para[3]);      //Paramagnetic correlation energy
    scalar_t ec_ferro = ec(rs, ferro[0], ferro[1], ferro[2], ferro[3]); //Ferromagnetic correlation energy
    scalar_t f0 = 1.7099209341613656175639627762446829382052199078228271066178937217;
    return ec_para + alpha_c(rs) * f_z(zeta) / f0  * (1. - pow(zeta,4.)) + (ec_ferro - ec_para) * f_z(zeta) * pow(zeta,4.);
}

// Finite difference of Spin correlation energy
template <typedef scalar_t>
scalar_t finite_diff_ecs_up(scalar_t n_up, scalar_t n_down, const scalar_t & para[4], const scalar_t & ferro[4]){
    scalar_t eps = std::numeric_limits<scalar_t>::epsilon()*(1. + n_up);
    return (ecs(n_up+eps, n_down, para, ferro) - ecs(n_up-eps, n_down, para, ferro)) / (2. * eps);
}
template <typedef scalar_t>
scalar_t finite_diff_ecs_down(scalar_t n_up, scalar_t n_down, const scalar_t & para[4], const scalar_t & ferro[4]){
    scalar_t eps = std::numeric_limits<scalar_t>::epsilon()*(1. + n_down);
    return (ecs(n_up, n_down+eps, para, ferro) - ecs(n_up, n_down-eps, para, ferro)) / (2. * eps);
}

// Local Spin Density Approximation correlation potential
template <typedef scalar_t>
DualVector<scalar_t> LSDACorrelationPot(const DualVector<scalar_t> & density){
    // Copying structure
    DualVector<scalar_t> correlation = density.EmptyLike();

    // Parameters {x0, b, c, Q}
    // Paramagnetic parameters
    scalar_t para[4] = {-0.10498, 3.72744, 12.9352, pow(4.*12.9352-pow(3.72744,2.),1./2.)};

    // Ferromagnetic parameters
    scalar_t ferro[4] = {-0.32500, 7.06042, 18.057, pow(4.*18.057-pow(7.06042,2.),1./2.)};

    // Calculate potential as finite difference of spin correlation energy ecs()
    // Loop over tiles
    int64_t mt = correlation.GetUp().mt();
    for (int64_t i = 0; i < mt; i++){
        // Make sure tile is local
        if (correlation.GetUp().tileIsLocal(i,0) &&
            density.GetUp().tileIsLocal(i,0)){
            // Access tiles
            auto tile_up = correlation.GetUp()(i, 0);
            auto tile_down = correlation.GetDown()(i, 0);
            auto tiledata_up = tile_up.data();
            auto tiledata_down = tile_down.data();

            // Access density tiles
            auto tile_up_n = density.GetUp()(i, 0);
            auto tile_down_n = density.GetDown()(i, 0);
            auto tiledata_up_n = tile_up.data();
            auto tiledata_down_n = tile_down.data();

            // Loop over tile
            for (int64_t ii; ii < tile_up.mb(); ii++){
                tiledata_up[ii] = finite_diff_ecs_up(tiledata_up_n[ii], tiledata_down_n[ii], para, ferro);
                tiledata_down[ii] = finite_diff_ecs_down(tiledata_up_n[ii], tiledata_down_n[ii], para, ferro);
            }
        }
    }

    return correlation;
}

template <typedef scalar_t>
DualVector<scalar_t> LSDA(const DualVector<scalar_t> & density, const Vector<scaler_t> & lattice){
    // Set the exchange potential
    DualVector<scale_t> lsda = LSDAExchangePot(density);

    // Calculate the correlation term
    lsda.Add(LSDACorrelationPot(density));

    return lsda;
}

/* ---- Exchange-Correlations dictionary ---- */
template <typedef scalar_t>
using func = DualVector<scalar_t>(*)(const DualVector<scalar_t> &, const Vector<scaler_t> &);

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



