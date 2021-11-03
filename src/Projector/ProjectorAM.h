#ifndef PROJECTORAM_H
#define PROJECTORAM_H

#include "Projector.h"
#include "Params.h"



//----------------------------------------------------------------------------------------------------------------------
//! class ProjectorAM: defines a virtual method for projection in 1d3v simulations
//----------------------------------------------------------------------------------------------------------------------
class ProjectorAM : public Projector
{

public:
    //! Constructor for ProjectorAM
    ProjectorAM( Params &params, Patch *patch );
    virtual ~ProjectorAM() {};
    
    virtual void mv_win( unsigned int shift )
    {
        i_domain_begin_+=shift;
    }
    virtual void setMvWinLimits( unsigned int shift )
    {
        i_domain_begin_ = shift;
    }
    
protected:
    double dr;
    double dt;
    //! Inverse of the spatial step 1/dx
    double dl_inv_;
    double dr_inv_;
    double dl_ov_dt_;
    double dr_ov_dt_;
    double one_ov_dt;
    int i_domain_begin_;
    int j_domain_begin_;
    int nscellr_;
    unsigned int Nmode_;
    unsigned int nprimr_, npriml_;
    int oversize_[2];
    double dq_inv_[2];
    //! Inverse radius
    double *invR, *invRd;
    
};

#endif

