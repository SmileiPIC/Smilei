#ifndef IonizationFactory_H
#define IonizationFactory_H

#include "Ionization.h"
#include "IonizationFromRate.h"
#include "IonizationTunnelEnvelopeAveraged.h"
#include "Params.h"
#include "Species.h"
#include "IonizationTunnel.h"
#include "Tools.h"

//! this class create and associate the right ionization model to species
class IonizationFactory
{
   public:
    static Ionization *create(Params &params, Species *species)
    {
        Ionization *Ionize = NULL;
        std::string model = species->ionization_model_;

        if (species->max_charge_ > (int)species->atomic_number_) {
            ERROR("Charge > atomic_number for species " << species->name_);
        }
        if (species->particles->is_test) {
            ERROR("Cannot ionize test species " << species->name_);
        }

        if (model == "tunnel_envelope_averaged") {
            if (!params.Laser_Envelope_model) {
                ERROR("The ionization model tunnel_envelope_averaged needs a laser envelope");
            }
            Ionize = new IonizationTunnelEnvelopeAveraged(params, species);
        } else if (params.Laser_Envelope_model) {
            ERROR("The ionization model for species interacting with envelope is tunnel_envelope_averaged");
        }

        if (model == "from_rate") {
            Ionize = new IonizationFromRate(params, species);
        } else if (model == "tunnel") {
            Ionize = new IonizationTunnel<0>(params, species); // Tunnel, the original included in Smilei
        } else if (model == "tunnel_full_PPT") {  
            Ionize = new IonizationTunnel<1>(params, species); // FullPPT
        } else if (model == "tunnel_TL") {  
            Ionize = new IonizationTunnel<2>(params, species); // Tong&Ling
        } else if (model == "tunnel_BSI") {  
            Ionize = new IonizationTunnel<3>(params, species); // BSI
        }

        return Ionize;
    }
};

#endif
