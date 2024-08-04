#ifndef TEMPLATE_H
#define TEMPLATE_H

#include <cmath>
#include <functional>
#include <vector>

#include "Ionization.h"
#include "Tools.h"

// because defining in header
#include "IonizationTables.h"
#include "Particles.h"
#include "Species.h"

class Particles;

struct electricFields
{
    double *x;
    double *y;
    double *z;
    double inv;  // inverse
    double abs;  // absolute value in atomic units
};

template <int Model>
class IonizationTunnel : public Ionization
{
   public:
    inline IonizationTunnel(Params &params, Species *species);

    inline void operator()(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                    int ipart_ref = 0) override;

    //! method for tunnel ionization with tasks
    inline void ionizationTunnelWithTasks(Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *,
                                   int, int, double *b_Jx, double *b_Jy, double *b_Jz, int ipart_ref = 0) override;

   private:
    template <int RateId>
    inline void monteCarloRoutine(Particles *, unsigned int, Patch *, Projector *, const unsigned int,
                                  const electricFields, vector<double> &, vector<double> &);
    template <int RateId>
    inline double ionizationRate(const int Z, const electricFields E);

    // To be conditionally prepared
    // FullPPT
    std::vector<double> Magnetic_quantum_number;

    // Tong&Ling
    double ionization_tl_parameter_;
    std::vector<double> lambda_tunnel;

    // BSI
    const double au_to_eV = 27.2116;
    const double IH = 13.598434005136;
};

template <int Model>
void IonizationTunnel<Model>::ionizationTunnelWithTasks( Particles *particles, unsigned int ipart_min, unsigned int ipart_max, 
                                                  vector<double> *Epart, Patch *patch, Projector *Proj, int ibin, int bin_shift, 
                                                  double *b_Jx, double *b_Jy, double *b_Jz, int ipart_ref )
{

    unsigned int Z, Zp1, newZ, k_times;
    double TotalIonizPot, E, invE, factorJion, delta, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    vector<double> IonizRate_tunnel( atomic_number_ ), Dnom_tunnel( atomic_number_ );
    LocalFields Jion;
    double factorJion_0 = au_to_mec2 * EC_to_au*EC_to_au * invdt;
    
    int nparts = Epart->size()/3;
    double *Ex = &( ( *Epart )[0*nparts] );
    double *Ey = &( ( *Epart )[1*nparts] );
    double *Ez = &( ( *Epart )[2*nparts] );
    
    for( unsigned int ipart=ipart_min ; ipart<ipart_max; ipart++ ) {
    
        // Current charge state of the ion
        Z = ( unsigned int )( particles->charge( ipart ) );
        
        // If ion already fully ionized then skip
        if( Z==atomic_number_ ) {
            continue;
        }
        
        // Absolute value of the electric field normalized in atomic units
        E = EC_to_au * sqrt( pow( *( Ex+ipart-ipart_ref ), 2 )
                             +pow( *( Ey+ipart-ipart_ref ), 2 )
                             +pow( *( Ez+ipart-ipart_ref ), 2 ) );
        if( E<1e-10 ) {
            continue;
        }
        
        // --------------------------------
        // Start of the Monte-Carlo routine
        // --------------------------------
        
        invE = 1./E;
        factorJion = factorJion_0 * invE*invE;
        delta      = gamma_tunnel[Z]*invE;
        ran_p = patch->rand_->uniform();
        IonizRate_tunnel[Z] = beta_tunnel[Z] * exp( -delta*one_third + alpha_tunnel[Z]*log( delta ) );
        
        // Total ionization potential (used to compute the ionization current)
        TotalIonizPot = 0.0;
        
        // k_times will give the nb of ionization events
        k_times = 0;
        Zp1=Z+1;
        
        if( Zp1 == atomic_number_ ) {
            // if ionization of the last electron: single ionization
            // -----------------------------------------------------
            if( ran_p < 1.0 -exp( -IonizRate_tunnel[Z]*dt ) ) {
                TotalIonizPot += Potential[Z];
                k_times        = 1;
            }
        
        } else {
            // else : multiple ionization can occur in one time-step
            //        partial & final ionization are decoupled (see Nuter Phys. Plasmas)
            // -------------------------------------------------------------------------
        
            // initialization
            Mult = 1.0;
            Dnom_tunnel[0]=1.0;
            Pint_tunnel = exp( -IonizRate_tunnel[Z]*dt ); // cummulative prob.
        
            //multiple ionization loop while Pint_tunnel < ran_p and still partial ionization
            while( ( Pint_tunnel < ran_p ) and ( k_times < atomic_number_-Zp1 ) ) {
                newZ = Zp1+k_times;
                delta = gamma_tunnel[newZ]*invE;
                IonizRate_tunnel[newZ] = beta_tunnel[newZ]
                                         *                        exp( -delta*one_third+alpha_tunnel[newZ]*log( delta ) );
                D_sum = 0.0;
                P_sum = 0.0;
                Mult  *= IonizRate_tunnel[Z+k_times];
                for( unsigned int i=0; i<k_times+1; i++ ) {
                    Dnom_tunnel[i]=Dnom_tunnel[i]/( IonizRate_tunnel[newZ]-IonizRate_tunnel[Z+i] );
                    D_sum += Dnom_tunnel[i];
                    P_sum += exp( -IonizRate_tunnel[Z+i]*dt )*Dnom_tunnel[i];
                }
                Dnom_tunnel[k_times+1] -= D_sum;
                P_sum                   = P_sum + Dnom_tunnel[k_times+1]*exp( -IonizRate_tunnel[newZ]*dt );
                Pint_tunnel             = Pint_tunnel + P_sum*Mult;
        
                TotalIonizPot += Potential[Z+k_times];
                k_times++;
            }//END while
        
            // final ionization (of last electron)
            if( ( ( 1.0-Pint_tunnel )>ran_p ) && ( k_times==atomic_number_-Zp1 ) ) {
                TotalIonizPot += Potential[atomic_number_-1];
                k_times++;
            }
        }//END Multiple ionization routine
        
        // Compute ionization current
        if (b_Jx != NULL){  // For the moment ionization current is not accounted for in AM geometry
            factorJion *= TotalIonizPot;
            Jion.x = factorJion * *( Ex+ipart );
            Jion.y = factorJion * *( Ey+ipart );
            Jion.z = factorJion * *( Ez+ipart );
        
            Proj->ionizationCurrentsForTasks( b_Jx, b_Jy, b_Jz, *particles, ipart, Jion, bin_shift );
        }
        
        // Creation of the new electrons
        // (variable weights are used)
        // -----------------------------
        if( k_times !=0 ) {
            new_electrons_per_bin[ibin].createParticle();
            int idNew = new_electrons_per_bin[ibin].size() - 1;//cout<<"ibin "<<ibin<<"size "<<new_electrons_per_bin[ibin].size()<<"capacity "<<new_electrons_per_bin[ibin].capacity()<<"\n"<<endl;
            for( unsigned int i=0; i<new_electrons_per_bin[ibin].dimension(); i++ ) {
                new_electrons_per_bin[ibin].position( i, idNew )=particles->position( i, ipart );
            }
            for( unsigned int i=0; i<3; i++ ) {
                new_electrons_per_bin[ibin].momentum( i, idNew ) = particles->momentum( i, ipart )*ionized_species_invmass;
            }
            new_electrons_per_bin[ibin].weight( idNew )=double( k_times )*particles->weight( ipart );
            new_electrons_per_bin[ibin].charge( idNew )=-1;
            
            if( save_ion_charge_ ) {
                ion_charge_per_bin_[ibin].push_back( particles->charge( ipart ) );
            }
            
            // // Increase the charge of the particle
            particles->charge( ipart ) += k_times;
        }
        
    } // Loop on particles
}

template <int Model>
template <int RateId>
inline void IonizationTunnel<Model>::monteCarloRoutine(Particles *particles, unsigned int ipart, Patch *patch, Projector *Proj,
                                               const unsigned int Z, const electricFields E,
                                               vector<double> &IonizRate_tunnel, vector<double> &Dnom_tunnel)
{
    double TotalIonizPot, factorJion, ran_p, Mult, D_sum, P_sum, Pint_tunnel;
    LocalFields Jion;
    const double factorJion_0 = au_to_mec2 * EC_to_au * EC_to_au * invdt;
    unsigned int newZ, Zp1, k_times;
    factorJion = factorJion_0 * E.inv * E.inv;
    ran_p = patch->rand_->uniform();
    IonizRate_tunnel[Z] = ionizationRate<RateId>(Z, E);

    // Total ionization potential (used to compute the ionization current)
    TotalIonizPot = 0.0;

    // k_times will give the nb of ionization events
    k_times = 0;
    Zp1 = Z + 1;

    if (Zp1 == atomic_number_)
    {
        // if ionization of the last electron: single ionization
        // -----------------------------------------------------
        if (ran_p < 1.0 - exp(-IonizRate_tunnel[Z] * dt))
        {
            TotalIonizPot += Potential[Z];
            k_times = 1;
        }
    }
    else
    {
        // else : multiple ionization can occur in one time-step
        //        partial & final ionization are decoupled (see Nuter Phys.
        //        Plasmas)
        // -------------------------------------------------------------------------

        // initialization
        Mult = 1.0;
        Dnom_tunnel[0] = 1.0;
        Pint_tunnel = exp(-IonizRate_tunnel[Z] * dt);  // cummulative prob.

        // multiple ionization loop while Pint_tunnel < ran_p and still partial
        // ionization
        while ((Pint_tunnel < ran_p) and (k_times < atomic_number_ - Zp1))
        {
            newZ = Zp1 + k_times;
            IonizRate_tunnel[newZ] = ionizationRate<RateId>(newZ, E);
            D_sum = 0.0;
            P_sum = 0.0;
            Mult *= IonizRate_tunnel[Z + k_times];
            for (unsigned int i = 0; i < k_times + 1; i++)
            {
                Dnom_tunnel[i] = Dnom_tunnel[i] / (IonizRate_tunnel[newZ] - IonizRate_tunnel[Z + i]);
                D_sum += Dnom_tunnel[i];
                P_sum += exp(-IonizRate_tunnel[Z + i] * dt) * Dnom_tunnel[i];
            }
            Dnom_tunnel[k_times + 1] -= D_sum;
            P_sum = P_sum + Dnom_tunnel[k_times + 1] * exp(-IonizRate_tunnel[newZ] * dt);
            Pint_tunnel = Pint_tunnel + P_sum * Mult;

            TotalIonizPot += Potential[Z + k_times];
            k_times++;
        }  // END while

        // final ionization (of last electron)
        if (((1.0 - Pint_tunnel) > ran_p) && (k_times == atomic_number_ - Zp1))
        {
            TotalIonizPot += Potential[atomic_number_ - 1];
            k_times++;
        }
    }  // END Multiple ionization routine

    // Compute ionization current
    if (patch->EMfields->Jx_ != NULL)
    {  // For the moment ionization current is
       // not accounted for in AM geometry
        factorJion *= TotalIonizPot;
        Jion.x = factorJion * *(E.x + ipart);
        Jion.y = factorJion * *(E.y + ipart);
        Jion.z = factorJion * *(E.z + ipart);

        Proj->ionizationCurrents(patch->EMfields->Jx_, patch->EMfields->Jy_, patch->EMfields->Jz_, *particles, ipart, Jion);
    }

    // Creation of the new electrons
    // (variable weights are used)
    // -----------------------------

    if (k_times != 0)
    {
        new_electrons.createParticle();
        int idNew = new_electrons.size() - 1;
        for (unsigned int i = 0; i < new_electrons.dimension(); i++)
        {
            new_electrons.position(i, idNew) = particles->position(i, ipart);
        }
        for (unsigned int i = 0; i < 3; i++)
        {
            new_electrons.momentum(i, idNew) = particles->momentum(i, ipart) * ionized_species_invmass;
        }
        new_electrons.weight(idNew) = double(k_times) * particles->weight(ipart);
        new_electrons.charge(idNew) = -1;

        if (save_ion_charge_)
        {
            ion_charge_.push_back(particles->charge(ipart));
        }

        // Increase the charge of the particle
        particles->charge(ipart) += k_times;
    }
}

template <int Model>
inline void IonizationTunnel<Model>::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                                 vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)

{
    unsigned int Z;
    electricFields E;
    vector<double> IonizRate_tunnel(atomic_number_), Dnom_tunnel(atomic_number_, 0.);

    int nparts = Epart->size() / 3;
    E.x = &((*Epart)[0 * nparts]);
    E.y = &((*Epart)[1 * nparts]);
    E.z = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++)
    {
        // Current charge state of the ion
        Z = (unsigned int)(particles->charge(ipart));

        // If ion already fully ionized then skip
        if (Z == atomic_number_)
        {
            continue;
        }

        // Absolute value of the electric field normalized in atomic units
        E.abs = EC_to_au * sqrt(pow(*(E.x + ipart - ipart_ref), 2) + pow(*(E.y + ipart - ipart_ref), 2) +
                                pow(*(E.z + ipart - ipart_ref), 2));
        if (E.abs < 1e-10)
        {
            continue;
        }
        E.inv = 1. / E.abs;

        monteCarloRoutine<0>(particles, ipart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel);
    }  // Loop on particles
}

template <int Model>
template <int RateId>
inline double IonizationTunnel<Model>::ionizationRate(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

// IonizationTunnel : 0
template <>
inline IonizationTunnel<0>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (unsigned int Z = 0; Z < atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

template <>
template <>
inline double IonizationTunnel<0>::ionizationRate<0>(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

// IonizationTunnelFullPPT: 1

template <>
inline IonizationTunnel<1>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel Ionizaton class");

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    Magnetic_quantum_number.resize(atomic_number_);

    for (unsigned int Z = 0; Z < atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);
        Magnetic_quantum_number[Z] = IonizationTables::magnetic_atomic_number(atomic_number_, Z);

        DEBUG("Potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]
                            << " M.q.num: " << Magnetic_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        double abs_m = abs(Magnetic_quantum_number[Z]);
        alpha_tunnel[Z] = cst - 1.0 - abs_m;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0 * tgamma(Azimuthal_quantum_number[Z] + abs_m + 1) /
                         (tgamma(abs_m + 1) * tgamma(Azimuthal_quantum_number[Z] - abs_m + 1));
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }

    DEBUG("Finished Creating the Tunnel Ionizaton class");
}

template <>
template <>
inline double IonizationTunnel<1>::ionizationRate<0>(const int Z, const electricFields E)
{
    double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta));
}

// Tong&Ling: 2

template <>
inline IonizationTunnel<2>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tong-Lin Tunnel Ionizaton class");

    ionization_tl_parameter_ =
        species->ionization_tl_parameter_;  // species->ionization_tl_parameter_ is double
                                            // Varies from 6 to 9. This is the alpha parameter in Tong-Lin exponential, see Eq. (6) in [M F Ciappina and S V Popruzhenko 2020 Laser Phys. Lett. 17 025301 2020].
    lambda_tunnel.resize(atomic_number_);

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (int Z = 0; Z < (int)atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);
        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
        lambda_tunnel[Z] = ionization_tl_parameter_ * pow(cst, 2) / gamma_tunnel[Z];
    }

    DEBUG("Finished Creating the Tong-Lin Tunnel Ionizaton class");
}

template <>
template <>
inline double IonizationTunnel<2>::ionizationRate<0>(const int Z, const electricFields E)
{
    const double delta = gamma_tunnel[Z] * E.inv;
    return beta_tunnel[Z] * exp(-delta * one_third + alpha_tunnel[Z] * log(delta) - E.abs * lambda_tunnel[Z]);
}

// BSI: 3

template <>
inline IonizationTunnel<3>::IonizationTunnel(Params &params, Species *species) : Ionization(params, species)
{
    DEBUG("Creating the Tunnel BSI Ionizaton class");

    // Ionization potential & quantum numbers (all in atomic units 1 au = 27.2116 eV)
    for (int Z = 0; Z < (int)atomic_number_; Z++)
    {
        DEBUG("Z : " << Z);

        Potential[Z] = IonizationTables::ionization_energy(atomic_number_, Z) * eV_to_au;
        Azimuthal_quantum_number[Z] = IonizationTables::azimuthal_atomic_number(atomic_number_, Z);

        DEBUG("potential: " << Potential[Z] << " Az.q.num: " << Azimuthal_quantum_number[Z]);

        double cst = ((double)Z + 1.0) * sqrt(2.0 / Potential[Z]);
        alpha_tunnel[Z] = cst - 1.0;
        beta_tunnel[Z] = pow(2, alpha_tunnel[Z]) * (8. * Azimuthal_quantum_number[Z] + 4.0) / (cst * tgamma(cst)) *
                         Potential[Z] * au_to_w0;
        gamma_tunnel[Z] = 2.0 * pow(2.0 * Potential[Z], 1.5);
    }
    DEBUG("Finished Creating the BSI Tunnel Ionizaton class");
}

// BSI Linear
template <>
template <>
inline double IonizationTunnel<3>::ionizationRate<1>(const int Z, const electricFields E)
{
    const double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
    return au_to_w0 * (0.8 * E.abs * pow(ratio_of_IPs, 0.5));
}

// BSI Linear
template <>
template <>
inline double IonizationTunnel<3>::ionizationRate<2>(const int Z, const electricFields E)
{
    const double ratio_of_IPs = IH / IonizationTables::ionization_energy(atomic_number_, Z);
    return au_to_w0 * (2.4 * (pow(E.abs, 2)) * pow(ratio_of_IPs, 2));
}

template <>
inline void IonizationTunnel<3>::operator()(Particles *particles, unsigned int ipart_min, unsigned int ipart_max,
                             vector<double> *Epart, Patch *patch, Projector *Proj, int ipart_ref)
{
    unsigned int Z;
    electricFields E;
    vector<double> IonizRate_tunnel(atomic_number_),
        Dnom_tunnel(atomic_number_);  // For MonteCarlo procedure using the tunnel
                                      // ioniz. rate.
    vector<double> IonizRate_BSI_linear(atomic_number_),
        Dnom_BSI_linear(atomic_number_);  // For MonteCarlo procedure using the BSI
                                          // ioniz. rate.
    vector<double> IonizRate_BSI_quadratic(atomic_number_), Dnom_BSI_quadratic(atomic_number_);

    int nparts = Epart->size() / 3;  // Over 3 because there are 3 Dims. Epart is like [part1_Ex, part2_Ex,
                                     // ... partnparts_Ex, part1_Ey, part2_Ey,...partnparts_Ey, part1_Ez,
                                     // part2_Ez, partnparts_Ez]}
    E.x = &((*Epart)[0 * nparts]);
    E.y = &((*Epart)[1 * nparts]);
    E.z = &((*Epart)[2 * nparts]);

    for (unsigned int ipart = ipart_min; ipart < ipart_max; ipart++)
    {
        Z = (unsigned int)(particles->charge(ipart));  // current charge state of the particle number ipart of species
                                                       // called ''species'' with atomic no = atomic_number_

        // If ion already fully ionized then skip and go to next particle
        if (Z == atomic_number_)
        {
            continue;
        }  // atom fully ionized already, no electrons remained to be stripped away.

        // Absolute value of the electric field in AtomicUnits (i.e. normalized to
        // (in units of) E_atomic = 5.1422*10^11 V/m)
        E.abs = EC_to_au * sqrt(pow(*(E.x + ipart - ipart_ref),
                                    2)  // (Ex+ipart-ipart_ref) points to the location in the
                                        // container Epart at which E-field_x for particle ipart
                                        // sits.
                                + pow(*(E.y + ipart - ipart_ref),
                                      2)  // Similarly for y. Dereferencing it via *() means we
                                          // access and use the value sitting at that location.
                                + pow(*(E.z + ipart - ipart_ref),
                                      2));  // EC_to_au transforms from SMILEI-normalized units to
                                            // AtomicUnits.
        if (E.abs < 1e-10)
        {  // E is in Atomic Units
            continue;
        }

        // Used to be call to continuity_tool. After some shortening, the call was manually inlined
        double BSI_rate_quadratic = ionizationRate<2>(Z, E);
        double BSI_rate_linear = ionizationRate<1>(Z, E);
        double Tunnel_rate = ionizationRate<0>(Z, E);

        if (BSI_rate_quadratic >= BSI_rate_linear)
        {
            monteCarloRoutine<2>(particles, ipart, patch, Proj, Z, E, IonizRate_BSI_quadratic, Dnom_BSI_quadratic);
        }
        else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic)
        {
            monteCarloRoutine<1>(particles, ipart, patch, Proj, Z, E, IonizRate_BSI_linear, Dnom_BSI_linear);
        }
        else
        {
            monteCarloRoutine<0>(particles, ipart, patch, Proj, Z, E, IonizRate_tunnel, Dnom_tunnel);
        }
    }  // END loop on particles

}  // void IonizationTunnelBSI::operator()(arg1, arg2, ...) scope end.

#endif
