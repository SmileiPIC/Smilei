#ifndef ELECTROMAGN_H
#define ELECTROMAGN_H

#include <vector>
#include <string>
#include <map>

#include "Tools.h"
#include "LaserParams.h"

class PicParams;
class Species;
class Projector;
class Field;
class Laser;
class SmileiMPI;
class ElectroMagnBC;
class SimWindow;

//! class ElectroMagn: generic class containing all information on the electromagnetic fields and currents

class ElectroMagn
{

public:

    std::vector<unsigned int> dimPrim;
    std::vector<unsigned int> dimDual;

    std::vector<unsigned int> index_bc_min;
    std::vector<unsigned int> index_bc_max;

    //! time-step (from picparams)
    const double timestep;
    
    //! cell length (from picparams)
    const std::vector<double> cell_length;

    //! \todo Generalise this to none-cartersian geometry (e.g rz, MG & JD)

    //! x-component of the electric field
    Field* Ex_;

    //! y-component of the electric field
    Field* Ey_;

    //! z-component of the electric field
    Field* Ez_;

    //! x-component of the magnetic field
    Field* Bx_;

    //! y-component of the magnetic field
    Field* By_;

    //! z-component of the magnetic field
    Field* Bz_;

    //! x-component of the time-centered magnetic field
    Field* Bx_m;

    //! y-component of the time-centered magnetic field
    Field* By_m;

    //! z-component of the time-centered magnetic field
    Field* Bz_m;

    //! x-component of the total charge current
    Field* Jx_;

    //! y-component of the total charge current
    Field* Jy_;

    //! z-component of the total charge current
    Field* Jz_;

    //! Total charge density
    Field* rho_;

    //! time-average x-component of the electric field
    Field* Ex_avg;
    
    //! time-average y-component of the electric field
    Field* Ey_avg;
    
    //! time-average z-component of the electric field
    Field* Ez_avg;
    
    //! time-average x-component of the magnetic field
    Field* Bx_avg;
    
    //! time-average y-component of the magnetic field
    Field* By_avg;
    
    //! time-average z-component of the magnetic field
    Field* Bz_avg;


    //! Vector of charge density and currents for each species
    const unsigned int n_species;
    std::vector<Field*> Jx_s;
    std::vector<Field*> Jy_s;
    std::vector<Field*> Jz_s;
    std::vector<Field*> rho_s;
    //std::vector<double**> nJx_s, nJy_s, nJz_s, nrho_s;
    double ***nJx_s, ***nJy_s, ***nJz_s, ***nrho_s;
    //! Number of bins
    unsigned int nbin;

    //! nDim_field (from params)
    const unsigned int nDim_field;

    //! Volume of the single cell (from params)
    const double cell_volume;

    //! n_space (from params) always 3D
    const std::vector<unsigned int> n_space;

    //! Index of starting elements in arrays without duplicated borders
    //! By constuction 1 element is shared in primal field, 2 in dual
    //! 3 : Number of direction (=1, if dim not defined)
    //! 2 : isPrim/isDual
    unsigned int istart[3][2];
    //! Number of elements in arrays without duplicated borders
    unsigned int bufsize[3][2];

    //!\todo should this be just an integer???
    //! Oversize domain to exchange less particles (from params)
    const std::vector<unsigned int> oversize;

    //! Constructor for Electromagn
    ElectroMagn( PicParams &params, LaserParams &laser_params, SmileiMPI* smpi );

    //! Destructor for Electromagn
    virtual ~ElectroMagn();

    //! Method used to dump data contained in ElectroMagn
    void dump();

    //! Method used to initialize the total charge currents and densities
    virtual void restartRhoJ() = 0;
    //! Method used to initialize the total charge currents and densities of species
    virtual void restartRhoJs(unsigned int ispec, bool currents) = 0;
    void restartRhoJs(unsigned int ispec, unsigned int clrw);

    //! Method used to initialize the total charge density
    void initRhoJ(std::vector<Species*> vecSpecies, Projector* Proj);

    //! Method used to sum all species densities and currents to compute the total charge density and currents
    virtual void computeTotalRhoJ() = 0;
    virtual void addToGlobalRho(int ispec, unsigned int clrw) = 0;
    virtual void computeTotalRhoJs(unsigned int clrw) = 0;
    //! Method used to gather species densities and currents on a single array
    virtual void sumtwins() = 0;

    //! Method used to initialize the Maxwell solver
    virtual void solvePoisson(SmileiMPI* smpi) = 0;

    //! \todo check time_dual or time_prim (MG)
    //! method used to solve Maxwell's equation (takes current time and time-step as input parameter)
    void solveMaxwell(int itime, double time_dual, SmileiMPI* smpi, PicParams &params, SimWindow* simWindow);
    virtual void solveMaxwellAmpere() = 0;
    virtual void solveMaxwellFaraday() = 0;
    virtual void saveMagneticFields() = 0;
    virtual void centerMagneticFields() = 0;

    void movingWindow_x(unsigned int shift, SmileiMPI *smpi);
    
    virtual void incrementAvgFields(unsigned int time_step, unsigned int ntime_step_avg) = 0;
        
    //! compute Poynting on borders
    virtual void computePoynting() = 0;
    
    //! pointing vector on borders 
    //! 1D: poynting[0][0]=left , poynting[1][0]=right
    //! 2D: poynting[0][0]=west , poynting[1][0]=east
    //!     poynting[1][0]=south, poynting[1][0]=north
    std::vector<double> poynting[2];

    //same as above but instantaneous
    std::vector<double> poynting_inst[2];

    //! Check if norm of charge denisty is not null
    bool isRhoNull(SmileiMPI* smpi);

    double computeNRJ(unsigned int shift, SmileiMPI *smpi);
    double getLostNrjMW() const {return nrj_mw_lost;}
    
    double getNewFieldsNRJ() const {return nrj_new_fields;}
    void reinitDiags() {
	nrj_mw_lost = 0.;
	nrj_new_fields = 0.;
    }

private:
    
    //! Vector of boundary-condition per side for the fields
    std::vector<ElectroMagnBC*> emBoundCond;

    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new fields
    double nrj_new_fields;

};

#endif
