
#ifndef Species_rrll_H
#define Species_rrll_H

#include "Species.h"
#include <string>

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PicParams;

class Species_rrll : public Species {
public:
  Species_rrll(PicParams*, int, SmileiMPI*);
  ~Species_rrll();
  void dynamic(double, ElectroMagn* Champs, Interpolator* Interp, Projector* proj, SmileiMPI* smpi);
private:
};

#endif

