#include "IonizationTunnel.h"

IonizationTunnel::IonizationTunnel(PicParams *params, int ispec) : Ionization(params, ispec) {
	DEBUG("Creating the Tunnel Ionizaton class");
}


void IonizationTunnel::operator() (Particle* part, LocalFields Epart) {

}