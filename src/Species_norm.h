
#ifndef SPECIESRRLL_H
#define SPECIESRRLL_H

#include "Species.h"
#include <string>

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PicParams;

class Species_norm : public Species {
public:
  Species_norm(PicParams*, int);
  ~Species_norm();
  //void dynamic(ElectroMagn* Champs, Pusher* ppush, Interpolator* Interp, Projector* proj);
private:
};

#endif

