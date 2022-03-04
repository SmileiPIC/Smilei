#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    
    //! Velocity of the Center-Of-Mass, expressed in the lab frame
    double COM_vx, COM_vy, COM_vz;
    
    //! Lorentz factor of the COM, expressed in the lab frame
    double COM_gamma;
    
    //! Momentum of the particles expressed in the COM frame
    double px_COM, py_COM, pz_COM, p_COM;
    
    //! Lorentz factors expressed in the COM frame
    double gamma1_COM, gamma2_COM;
};

#endif
