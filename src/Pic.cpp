////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////                                                                                                                ////
////                                                                                                                ////
////                                   PARTICLE-IN-CELL CODE SMILEI                                                 ////
////                    Simulation of Matter Irradiated by Laser at Extreme Intensity                               ////
////                                                                                                                ////
////                          Cooperative OpenSource Object-Oriented Project                                        ////
////                                      from the Plateau de Saclay                                                ////
////                                          started January 2013                                                  ////
////                                                                                                                ////
////                                                                                                                ////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include "Pic.h"
#include "Tools.h"
#include "PicParams.h"
#include "Species.h"
#include "Species_norm.h"
#include "Species_rrll.h"
#include "ElectroMagn.h"
#include "ElectroMagn1D.h"
#include "Pusher.h"
#include "PusherBoris.h"
#include "PusherSklv.h"
#include "Interpolator.h"
#include "Interpolator1D2Order.h"
#include "Interpolator1D3Order.h"
#include "Projector.h"
#include "Projector1D2Order.h"
#include "Field3D.h"
#include "Field1D.h"

#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;



// ------------------------------------------------------------------------------------------------------------------ //
//                                                   MAIN CODE
// ------------------------------------------------------------------------------------------------------------------ //
int main (int argc, char* argv[])
{
    
    
    // -------------------------
    // Simulation Initialisation
    // -------------------------
    
    // Check for namelist (input file)
    if (argc<2) ERROR("No namelists given!");
    string namelist=argv[1];
    
    // Send information on current simulation
	MESSAGE("------------------------------------------");
	MESSAGE(" Version : " << __VERSION DEBUGEXEC(<< " DEBUG") << " Compiled : " << __DATE__ << " " << __TIME__);
	MESSAGE("------------------------------------------");
    MESSAGE(" Namelist  : " << namelist);
    MESSAGE("------------------------------------------");
    
    // Read simulation parameters
    PicParams params(namelist);
    
    // Randomize the seed for simulations running in release mode
    //! \todo{Save the seed in case one wants to re-run the exact same simulation (MG)}
    RELEASEEXEC(srand (time(NULL)));
    
    
    
    // -------------------------------------------
    // Declaration of the main objects & operators
    // -------------------------------------------
    
    // object containing the electromagnetic fields (virtual)
    ElectroMagn* EMfields = NULL;
    
    // interpolation operator (virtual)
    Interpolator* Interp = NULL;
    
    // projection operator (virtual)
    Projector* Proj = NULL;
    
    // vector of Species (virtual)
    vector<Species*> vecSpecies;
    
    // species "dump" file
    //! \todo{Check if we keep going like that (MG)}
    ofstream ofile("dump", ios::out);
    
    
    
    // ------------------------------------------------------------------------------------
    // Initialize the vecSpecies object containing all information of the different Species
    // ------------------------------------------------------------------------------------
    vecSpecies.resize(params.n_species);
    
    for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
        MESSAGE(0,"Initializing Species "<<ispec);
        
        Species* sp = NULL;
        if (params.species_param[ispec].dynamics_type=="norm") {
            // Species with Boris dynamics
            sp = new Species_norm(&params, ispec);
        } else if (params.species_param[ispec].dynamics_type=="rrll") {
            // Species with Boris dynamics + Radiation Back-Reaction (using the Landau-Lifshitz formula)
            sp = new Species_rrLL(&params, ispec);
        }//endif
        
        //dump species at time 0
        sp->dump(ofile); ofile << endl;
        
        //save temporary species sp in vecSpecies
        vecSpecies[ispec] = sp;
        
        MESSAGE(0,sp->getNbrOfParticles() << " Particles of species " << ispec);
    }// END for ispec
    
    
    
    // ----------------------------------------------------------------------------
    // Initialize the electromagnetic fields and interpolation-projection operators
    // according to the simulation geometry
    // ----------------------------------------------------------------------------
    if ( params.geometry == "1d3v" )
    {
        
        // ---------------
        // 1d3v Simulation
        // ---------------
        EMfields   = new ElectroMagn1D(&params);
        if ( params.interpolation_order == 2 )
        {
            // 2nd order interpolation
            Interp = new Interpolator1D2Order(&params);
            Proj   = new Projector1D2Order(&params);
        }
    }
    else {
        ERROR( "Unknwon geometry : " << params.geometry );
    }//endif params.geometry
    
    
    
    
    // -----------------------------------
    // Inialize the electromagnetic fields
    // -----------------------------------
    
    //!\todo{Check & describe what is done here (MG)}
    EMfields->initRho(vecSpecies, Proj);
    EMfields->solvePoisson();
    
    
    
    // ------------------------------------------------------------------------
    // Initialise the simulation times time_prim at n=0 and time_dual at n=-1/2
    // ------------------------------------------------------------------------
    
    // time at integer time-steps (primal grid)
    double time_prim = 0.;
    // time at half-integer time-steps (dual grid)
    double time_dual = -0.5 * params.timestep;
    
    
    // ------------------------------------------------------------------
    //                     HERE STARTS THE PIC LOOP
    // ------------------------------------------------------------------
    MESSAGE(0,"Time-Loop is started: number of time-steps n_time =" << params.n_time);
    for (unsigned int itime=1 ; itime <= params.n_time ; itime++)
    {
        
        // calculate new times
        // -------------------
        time_prim += params.timestep;
        time_dual += params.timestep;
        
        
        // send message at given time-steps
        // --------------------------------
        //!\todo{Introduce a control parameter in PicParams (MG)}
        if (itime % 100 == 0)   MESSAGE(1,"Time: " << time_dual << " " << itime);
        
        
        // put density and currents to 0
        // -----------------------------
        EMfields->initRhoJ();
        
        
        // apply the PIC method
        // --------------------
        
        // for all particles of all species (see dunamic in Species.cpp)
        // (1) interpolate the fields at the particle position
        // (2) move the particle
        // (3) calculate the currents (charge conserving method)
        for (unsigned int ispec=0 ; ispec<params.n_species; ispec++)
        {
            DEBUG(2, "Dynamic Species "<<ispec );
            vecSpecies[ispec]->dynamics(time_dual, EMfields, Interp, Proj);
        }
        
        // calculate the longitudinal current using the charge conservation equation
        // EMfields->chargeConserving();
        
        // solve Maxwell's equations
        EMfields->solveMaxwell(time_dual, params.timestep);
        
        
        // call the various diagnostics
        // ----------------------------
        if (itime % 5000 == 0)
        {
            MESSAGE(1,"diags at time t=" << time_dual);
            EMfields->dump(&params);
/*            for (unsigned int ispec=0 ; ispec<params.n_species ; ispec++) {
                vecSpecies[ispec]->dump(ofile);
                ofile << endl;
            }
 */
            vecSpecies[0]->dump(ofile);
            ofile << endl;

        }
        
    }//END of the time loop
    
    MESSAGE(0,"End time loop");
    // ------------------------------------------------------------------
    //                      HERE ENDS THE PIC LOOP
    // ------------------------------------------------------------------
    
    
    
    // ------------------------------
    //  Cleanup & End the simulation
    // ------------------------------
    
    delete Proj;
    delete Interp;
    delete EMfields;
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    
    MESSAGE("------------------------------------------");
    MESSAGE("END " << namelist);
    MESSAGE("------------------------------------------");
	
    return 0;
    
}//END MAIN 


