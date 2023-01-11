#include "PyTools.h"
#include "DiagnosticScalar.h"
#include "VectorPatch.h"
#include "Tools.h"

#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;


DiagnosticScalar::DiagnosticScalar( Params &params, SmileiMPI *, Patch * = NULL ):
    latest_timestep( -1 )
{
    
    if( PyTools::nComponents( "DiagScalar" ) > 1 ) {
        ERROR( "Only one DiagScalar can be specified" );
    }
    
    if( PyTools::nComponents( "DiagScalar" ) > 0 ) {
    
        // get parameter "every" which describes a timestep selection
        timeSelection = new TimeSelection(
            PyTools::extract_py( "every", "DiagScalar", 0 ),
            "Scalars"
        );
        
        precision = 10;
        PyTools::extract( "precision", precision, "DiagScalar"  );
        PyTools::extractV( "vars", vars, "DiagScalar" );
        
        // copy from params remaining stuff
        res_time       = params.res_time;
        dt             = params.timestep;
        cell_volume    = params.cell_volume;
        patch_size_    = params.patch_size_;
        global_size_   = params.global_size_;
        
        filename = "scalars.txt";
    } else {
        timeSelection = new TimeSelection();
        filename = "";
    }
    
} // END DiagnosticScalar::DiagnosticScalar


DiagnosticScalar::~DiagnosticScalar()
{
    for( unsigned int i=0; i<allScalars.size(); i++ ) {
        delete allScalars[i];
    }
    delete timeSelection;
} // END DiagnosticScalar::#DiagnosticScalar


void DiagnosticScalar::openFile( Params &, SmileiMPI *smpi )
{
    if( !smpi->isMaster() || fout.is_open() ) {
        return;
    }
    
    fout.open( "scalars.txt" );
    
    if( !fout.is_open() ) {
        ERROR( "Can't open scalar.txt file" );
    }
    
} // END openFile


void DiagnosticScalar::closeFile()
{
    if( fout.is_open() ) {
        fout.close();
    }
    
} // END closeFile


unsigned int DiagnosticScalar::calculateWidth( string key )
{
    return 2 + max( ( ( unsigned int )key.length() ), precision+8 );
    // The +8 accounts for the dot and exponent in decimal representation)
}

Scalar_value *DiagnosticScalar::newScalar_SUM( string name )
{
    bool allow = allowedKey( name );
    unsigned int width = calculateWidth( name );
    Scalar_value *scalar = new Scalar_value( name, width, allow, &values_SUM );
    values_SUM.push_back( 0. );
    allScalars.push_back( scalar );
    return scalar;
}
Scalar_value_location *DiagnosticScalar::newScalar_MINLOC( string name )
{
    bool allow = allowedKey( name );
    unsigned int width = calculateWidth( name );
    val_index default_;
    default_.index = -1;
    default_.val = 0.;
    Scalar_value_location *scalar = new Scalar_value_location( name, name+"Cell", width, allow, &values_MINLOC, numeric_limits<double>::max() );
    values_MINLOC.push_back( default_ );
    allScalars.push_back( scalar );
    return scalar;
}
Scalar_value_location *DiagnosticScalar::newScalar_MAXLOC( string name )
{
    bool allow = allowedKey( name );
    unsigned int width = calculateWidth( name );
    val_index default_;
    default_.index = -1;
    default_.val = 0.;
    Scalar_value_location *scalar = new Scalar_value_location( name, name+"Cell", width, allow, &values_MAXLOC, numeric_limits<double>::lowest() );
    values_MAXLOC.push_back( default_ );
    allScalars.push_back( scalar );
    return scalar;
}


void DiagnosticScalar::init( Params &params, SmileiMPI *, VectorPatch &vecPatches )
{

    // Make the list of fields
    ElectroMagn *EMfields = vecPatches( 0 )->EMfields;
    vector<string> fields;
    unsigned int nmodes = 1;
    if( params.geometry != "AMcylindrical" ) {
        fields.push_back( EMfields->Ex_ ->name );
        fields.push_back( EMfields->Ey_ ->name );
        fields.push_back( EMfields->Ez_ ->name );
        fields.push_back( EMfields->Bx_m->name );
        fields.push_back( EMfields->By_m->name );
        fields.push_back( EMfields->Bz_m->name );
        fields.push_back( EMfields->Jx_ ->name );
        fields.push_back( EMfields->Jy_ ->name );
        fields.push_back( EMfields->Jz_ ->name );
        fields.push_back( EMfields->rho_->name );
        // add envelope-related fields
        if( params.Laser_Envelope_model ) {
            fields.push_back( EMfields->Env_A_abs_->name );
            fields.push_back( EMfields->Env_Chi_->name );
            fields.push_back( EMfields->Env_E_abs_->name );
            fields.push_back( EMfields->Env_Ex_abs_->name );
        }
    } else {
        ElectroMagnAM *emfields = static_cast<ElectroMagnAM *>( EMfields );
        nmodes = emfields->El_.size();
        for( unsigned int imode=0; imode < nmodes; imode++ ) {
            fields.push_back( emfields->El_[imode] ->name );
            fields.push_back( emfields->Er_[imode] ->name );
            fields.push_back( emfields->Et_[imode] ->name );
            fields.push_back( emfields->Bl_m[imode]->name );
            fields.push_back( emfields->Br_m[imode]->name );
            fields.push_back( emfields->Bt_m[imode]->name );
        }
        for( unsigned int imode=0; imode < nmodes; imode++ ) {
            fields.push_back( emfields->Jl_[imode]->name );
            fields.push_back( emfields->Jr_[imode]->name );
            fields.push_back( emfields->Jt_[imode]->name );
            fields.push_back( emfields->rho_AM_[imode]->name );
        }
        
    }
    
    // 1 - Prepare the booleans that tell which scalars are necessary to compute
    // -------------------------------------------------------------------------
    
    // General scalars
    necessary_Ubal_norm  = allowedKey( "Ubal_norm" );
    necessary_Ubal       = necessary_Ubal_norm || allowedKey( "Ubal" );
    necessary_Utot       = necessary_Ubal || allowedKey( "Utot" );
    necessary_Uexp       = necessary_Ubal || allowedKey( "Uexp" );
    necessary_Ukin       = necessary_Utot || allowedKey( "Ukin" );
    necessary_Uelm       = necessary_Utot || allowedKey( "Uelm" );
    necessary_Urad       = necessary_Utot || allowedKey( "Urad" );
    necessary_UmBWpairs  = allowedKey( "UmBWpairs" );
    necessary_Ukin_BC = necessary_Uexp || allowedKey( "Ukin_bnd" ) || allowedKey( "Ukin_out_mvw" ) || allowedKey( "Ukin_inj_mvw" );
    necessary_Uelm_BC = necessary_Uexp || allowedKey( "Uelm_bnd" ) || allowedKey( "Uelm_out_mvw" ) || allowedKey( "Uelm_inj_mvw" );
    // Species
    unsigned int nspec = vecPatches( 0 )->vecSpecies.size();
    necessary_species.resize( nspec, false );
    string species_name;
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if( ! vecPatches( 0 )->vecSpecies[ispec]->particles->is_test ) {
            species_name = vecPatches( 0 )->vecSpecies[ispec]->name_;
            necessary_species[ispec] = necessary_Ukin || allowedKey( Tools::merge( "Dens_", species_name ) )
                                       || allowedKey( Tools::merge( "Ntot_", species_name ) )
                                       || allowedKey( Tools::merge( "Zavg_", species_name ) )
                                       || allowedKey( Tools::merge( "Ukin_", species_name ) )
                                       || allowedKey( Tools::merge( "Urad_", species_name ) )
                                       || allowedKey( "UmBWpairs" );
        }
    }
    // Fields
    unsigned int nfield = nmodes * 6;
    necessary_fieldUelm.resize( nfield, false );
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        necessary_fieldUelm[ifield] = necessary_Uelm || allowedKey( Tools::merge( "Uelm_", fields[ifield] ) );
    }
    // Fields min/max
    nfield = fields.size();
    necessary_fieldMinMax.resize( nfield, false );
    necessary_fieldMinMax_any = false;
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        necessary_fieldMinMax[ifield] =
            allowedKey( Tools::merge( fields[ifield], "Min" ) )
            || allowedKey( Tools::merge( fields[ifield], "MinCell" ) )
            || allowedKey( Tools::merge( fields[ifield], "Max" ) )
            || allowedKey( Tools::merge( fields[ifield], "MaxCell" ) );
        if( necessary_fieldMinMax[ifield] ) {
            necessary_fieldMinMax_any = true;
        }
    }
    // Poynting flux
    unsigned int npoy  = EMfields->poynting[0].size() * EMfields->poynting[1].size();
    necessary_poy.resize( npoy );
    string poy_name;
    unsigned int k = 0;
    for( unsigned int j=0; j<2; j++ ) {
        for( unsigned int i=0; i<EMfields->poynting[j].size(); i++ ) {
            //if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            //else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            //else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
            necessary_poy[k] = necessary_Uelm_BC || allowedKey( poy_name ) || allowedKey( poy_name+"Inst" );
            k++;
        }
    }
    
    // 2 - Prepare the Scalar* objects that will contain the data
    // ----------------------------------------------------------
    
    values_SUM   .reserve( 13 + nspec*5 + 6 + 2*npoy );
    
    if( !params.Laser_Envelope_model ) {
        values_MINLOC.reserve( 10 );
        values_MAXLOC.reserve( 10 );
    } else {
        values_MINLOC.reserve( 14 );
        values_MAXLOC.reserve( 14 );
    }
    
    // General scalars
    Ubal_norm    = newScalar_SUM( "Ubal_norm" );
    Ubal         = newScalar_SUM( "Ubal" );
    Utot         = newScalar_SUM( "Utot" );
    Uexp         = newScalar_SUM( "Uexp" );
    Ukin         = newScalar_SUM( "Ukin" );
    Urad         = newScalar_SUM( "Urad" );          // Radiated energy
    UmBWpairs    = newScalar_SUM( "UmBWpairs" );     // multiphoton Breit-Wheeler
    Uelm         = newScalar_SUM( "Uelm" );
    Ukin_bnd     = newScalar_SUM( "Ukin_bnd" );
    Ukin_out_mvw = newScalar_SUM( "Ukin_out_mvw" );
    Ukin_inj_mvw = newScalar_SUM( "Ukin_inj_mvw" );
    Ukin_new     = newScalar_SUM( "Ukin_new" );
    Uelm_bnd     = newScalar_SUM( "Uelm_bnd" );
    Uelm_out_mvw = newScalar_SUM( "Uelm_out_mvw" );
    Uelm_inj_mvw = newScalar_SUM( "Uelm_inj_mvw" );
    
    // Scalars related to species
    sDens.resize( nspec, NULL );
    sNtot.resize( nspec, NULL );
    sZavg.resize( nspec, NULL );
    sUkin.resize( nspec, NULL );
    sUrad.resize( nspec, NULL );
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if( ! vecPatches( 0 )->vecSpecies[ispec]->particles->is_test ) {
            species_name = vecPatches( 0 )->vecSpecies[ispec]->name_;
            sDens[ispec] = newScalar_SUM( Tools::merge( "Dens_", species_name ) );
            sNtot[ispec] = newScalar_SUM( Tools::merge( "Ntot_", species_name ) );
            sZavg[ispec] = newScalar_SUM( Tools::merge( "Zavg_", species_name ) );
            sUkin[ispec] = newScalar_SUM( Tools::merge( "Ukin_", species_name ) );
            sUrad[ispec] = newScalar_SUM( Tools::merge( "Urad_", species_name ) );
        }
    }
    
    // Scalars related to field's electromagnetic energy
    nfield = nmodes * 6;
    fieldUelm.resize( nfield, NULL );
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        fieldUelm[ifield] = newScalar_SUM( Tools::merge( "Uelm_", fields[ifield] ) );
    }
    
    // Scalars related to fields min and max
    nfield = ( params.geometry == "AMcylindrical" ) ? 0 : fields.size();
    fieldMin.resize( nfield, NULL );
    fieldMax.resize( nfield, NULL );
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        if( necessary_fieldMinMax[ifield] ) {
            fieldMin[ifield] = newScalar_MINLOC( Tools::merge( fields[ifield], "Min" ) );
            fieldMax[ifield] = newScalar_MAXLOC( Tools::merge( fields[ifield], "Max" ) );
        }
    }
    
    // Scalars related to the Poynting flux
    poy    .resize( 2*npoy, NULL );
    poyInst.resize( 2*npoy, NULL );
    k = 0;
    for( unsigned int j=0; j<2; j++ ) {
        for( unsigned int i=0; i<EMfields->poynting[j].size(); i++ ) {
            //if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            //else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            //else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
            poy    [k] = newScalar_SUM( poy_name );
            poyInst[k] = newScalar_SUM( Tools::merge( poy_name, "Inst" ) );
            k++;
        }
    }
}



bool DiagnosticScalar::prepare( int itime )
{
    // At the right timestep, reset the scalars
    if( timeSelection->theTimeIsNow( itime ) ) {
        for( unsigned int iscalar=0 ; iscalar<allScalars.size() ; iscalar++ ) {
            allScalars[iscalar]->reset();
        }
    }
    
    // Scalars always run even if they don't dump
    return true;
} // END prepare


void DiagnosticScalar::run( Patch *patch, int itime, SimWindow * )
{

    // Must keep track of Poynting flux even without diag
    if( ! filename.empty() ) {
        patch->computePoynting();
    }
    
    // Compute all scalars when needed
    if( timeSelection->theTimeIsNow( itime ) && itime>latest_timestep ) {
        compute( patch, itime );
    }
    
} // END run


void DiagnosticScalar::write( int itime, SmileiMPI *smpi )
{
    if( smpi->isMaster() ) {
    
        if( timeSelection->theTimeIsNow( itime ) && itime>latest_timestep ) {
        
            unsigned int j, k, s = allScalars.size();
            
            fout << std::scientific << setprecision( precision );
            // At the beginning of the file, we write some headers
            if( fout.tellp()==ifstream::pos_type( 0 ) ) { // file beginning
                // First header: list of scalars, one by line
                fout << "# " << 1 << " time" << endl;
                j = 2;
                for( k=0; k<s; k++ ) {
                    if( allScalars[k]->allowed_ ) {
                        fout << "# " << j << " " << allScalars[k]->name_ << endl;
                        j++;
                        if( ! allScalars[k]->secondname_.empty() ) {
                            fout << "# " << j << " " << allScalars[k]->secondname_ << endl;
                            j++;
                        }
                    }
                }
                // Second header: list of scalars, but all in one line
                fout << "#\n#" << setw( precision+9 ) << "time";
                for( k=0; k<s; k++ ) {
                    if( allScalars[k]->allowed_ ) {
                        fout << setw( allScalars[k]->width_ ) << allScalars[k]->name_;
                        if( ! allScalars[k]->secondname_.empty() ) {
                            fout << setw( allScalars[k]->width_ ) << allScalars[k]->secondname_;
                        }
                    }
                }
                fout << endl;
            }
            // Each requested timestep, the following writes the values of the scalars
            fout << setw( precision+10 ) << itime/res_time;
            for( k=0; k<s; k++ ) {
                if( allScalars[k]->allowed_ ) {
                    fout << setw( allScalars[k]->width_ ) << ( double )*allScalars[k];
                    if( ! allScalars[k]->secondname_.empty() ) {
                        fout << setw( allScalars[k]->width_ ) << ( int )*static_cast<Scalar_value_location *>( allScalars[k] );
                    }
                }
            }
            fout << endl;
            
        }
        
    } // if smpi->isMaster
    
    latest_timestep = itime;
    
} // END write


//! Compute the various scalars when requested
void DiagnosticScalar::compute( Patch *patch, int )
{
    ElectroMagn *EMfields = patch->EMfields;
    std::vector<Species *> &vecSpecies = patch->vecSpecies;
    
    ElectroMagnAM *AMfields = dynamic_cast<ElectroMagnAM *>( patch->EMfields );
    bool AM = AMfields;
    
    // ------------------------
    // SPECIES-related energies
    // ------------------------
    double Ukin_=0.;             // total (kinetic) energy carried by particles (test-particles do not contribute)
    double Urad_=0.;             // total radiated energy by particles
    double UmBWpairs_=0;         // total enegy converted into pairs via the multiphoton Breit-Wheeler
    double Ukin_bnd_=0.;         // total energy lost by particles due to boundary conditions
    double Ukin_out_mvw_=0.;     // total energy lost due to particles being suppressed by the moving-window
    double Ukin_inj_mvw_=0.;     // total energy added due to particles created by the moving-window
    double Ukin_new_=0.;         // total energy added due to new particles (injector)
    
    // Compute scalars for each species
    for( unsigned int ispec=0; ispec<vecSpecies.size(); ispec++ ) {
        if( vecSpecies[ispec]->particles->is_test ) {
            continue;    // No scalar diagnostic for test particles
        }
        
        if( necessary_species[ispec] ) {
            double density=0.0;  // sum of weights of current species ispec
            double charge=0.0;   // sum of charges of current species ispec
            double ener_tot=0.0; // total kinetic energy of current species ispec
            
            unsigned int nPart=vecSpecies[ispec]->getNbrOfParticles(); // number of particles
            
            if( vecSpecies[ispec]->mass_ > 0 ) {
            
                for( unsigned int iPart=0 ; iPart<nPart; iPart++ ) {
                
                    density  += vecSpecies[ispec]->particles->weight( iPart );
                    charge   += vecSpecies[ispec]->particles->weight( iPart )
                                * ( double )vecSpecies[ispec]->particles->charge( iPart );
                    ener_tot += vecSpecies[ispec]->particles->weight( iPart )
                                * ( vecSpecies[ispec]->particles->LorentzFactor( iPart )-1.0 );
                }
                ener_tot *= vecSpecies[ispec]->mass_;
            } else if( vecSpecies[ispec]->mass_ == 0 ) {
                for( unsigned int iPart=0 ; iPart<nPart; iPart++ ) {
                
                    density  += vecSpecies[ispec]->particles->weight( iPart );
                    ener_tot += vecSpecies[ispec]->particles->weight( iPart )
                                * ( vecSpecies[ispec]->particles->momentumNorm( iPart ) );
                }
            }
            
            *sNtot[ispec] += ( double )nPart;
            *sDens[ispec] += density;
            *sZavg[ispec] += charge;
            *sUkin[ispec] += ener_tot;
            
            // incremement the total kinetic energy
            Ukin_ += ener_tot;
            
            // If radiation activated
            if( vecSpecies[ispec]->Radiate ) {
                *sUrad[ispec] += vecSpecies[ispec]->nrj_radiated_;
                Urad_         += vecSpecies[ispec]->nrj_radiated_;
            }
            
            // If multiphoton Breit-Wheeler activated for photons
            // increment the total pair energy from this process
            if( vecSpecies[ispec]->Multiphoton_Breit_Wheeler_process ) {
                UmBWpairs_ += vecSpecies[ispec]->nrj_radiated_;
            }
        }
        
        if( necessary_Ukin_BC ) {
            // particle energy lost due to boundary conditions
            Ukin_bnd_     += vecSpecies[ispec]->nrj_bc_lost;
            // particle energy lost due to moving window
            Ukin_out_mvw_ += vecSpecies[ispec]->nrj_mw_out;
            // particle energy added by moving window
            Ukin_inj_mvw_ += vecSpecies[ispec]->nrj_mw_inj;
        }
        
        // particle energy from new particles
        Ukin_new_ += vecSpecies[ispec]->nrj_new_part_;
        
    } // for ispec
    
    // Add the calculated energies to the data arrays
    if( necessary_Ukin ) {
        *Ukin += Ukin_;
    }
    if( necessary_Urad ) {
        *Urad += Urad_;
    }
    if( necessary_UmBWpairs ) {
        *UmBWpairs += UmBWpairs_;
    }
    if( necessary_Ukin_BC ) {
        *Ukin_bnd     += Ukin_bnd_     ;
        *Ukin_out_mvw += Ukin_out_mvw_ ;
        *Ukin_inj_mvw += Ukin_inj_mvw_ ;
    }
    *Ukin_new += Ukin_new_;
    
    // --------------------------------
    // ELECTROMAGNETIC-related energies
    // --------------------------------
    
    vector<Field *> fields;
    
    if( !AM ) {
        fields.push_back( EMfields->Ex_ );
        fields.push_back( EMfields->Ey_ );
        fields.push_back( EMfields->Ez_ );
        fields.push_back( EMfields->Bx_m );
        fields.push_back( EMfields->By_m );
        fields.push_back( EMfields->Bz_m );
    } else {
        unsigned int nmodes = AMfields->El_.size();
        for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
            fields.push_back( AMfields->El_[imode] );
            fields.push_back( AMfields->Er_[imode] );
            fields.push_back( AMfields->Et_[imode] );
            fields.push_back( AMfields->Bl_m[imode] );
            fields.push_back( AMfields->Br_m[imode] );
            fields.push_back( AMfields->Bt_m[imode] );
        }
    }
    
    double Uelm_ = 0.0; // total electromagnetic energy in the fields
    
    // loop on all electromagnetic fields
    unsigned int nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        if( necessary_fieldUelm[ifield] ) {
            Field *field = fields[ifield];
            
            // total energy in current field
            double Uem = 0.;
            if( ! AM ) {
                Uem = field->norm2( EMfields->istart, EMfields->bufsize );
            } else {
                Uem = 0.5 * AMfields->dr * dynamic_cast<cField2D*>( field )->norm2_cylindrical( AMfields->istart, AMfields->bufsize, AMfields->j_glob_ );
                if( ifield < 6 ) {
                    Uem *= 2.; // Factor 2 for mode 0
                }
            }
            Uem *= 0.5*cell_volume;
            
            *fieldUelm[ifield] += Uem;
            Uelm_ += Uem;
        }
    }
    
    // Total elm energy
    if( necessary_Uelm ) {
        *Uelm += Uelm_;
    }
    
    // Lost/added elm energies through the moving window
    if( necessary_Uelm_BC ) {
        // nrj lost with moving window (fields)
        *Uelm_out_mvw += EMfields->nrj_mw_out * 0.5*cell_volume;
        
        // nrj added due to moving window (fields)
        *Uelm_inj_mvw += EMfields->nrj_mw_inj * 0.5*cell_volume;
    }
    
    
    // ---------------------------
    // ELECTROMAGNETIC MIN and MAX
    // ---------------------------
    
    // add currents and density to fields
    fields.push_back( EMfields->Jx_ );
    fields.push_back( EMfields->Jy_ );
    fields.push_back( EMfields->Jz_ );
    fields.push_back( EMfields->rho_ );
    
    // add envelope-related fields
    if( EMfields->Env_A_abs_ != NULL ) {
        fields.push_back( EMfields->Env_A_abs_ );
        fields.push_back( EMfields->Env_Chi_ );
        fields.push_back( EMfields->Env_E_abs_ );
        fields.push_back( EMfields->Env_Ex_abs_ );
    }
    
    double fieldval;
    unsigned int i_min, j_min, k_min;
    unsigned int i_max, j_max, k_max;
    val_index minloc, maxloc;
    
    nfield = fields.size();

    // if AM, scalar on fields not managed
    if( AM ){
        nfield = 0;
    }

    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
    
        if( necessary_fieldMinMax[ifield] ) {
        
            Field *field = fields[ifield];
            
            vector<unsigned int> iFieldStart( 3, 0 ), iFieldEnd( 3, 1 ), iFieldGlobalSize( 3, 1 );
            for( unsigned int i=0 ; i<field->isDual_.size() ; i++ ) {
                iFieldStart[i] = EMfields->istart[i][field->isDual( i )];
                iFieldEnd [i] = iFieldStart[i] + EMfields->bufsize[i][field->isDual( i )];
                iFieldGlobalSize [i] = field->dims_[i];
            }
            
            unsigned int iifield= iFieldStart[2] + iFieldStart[1]*iFieldGlobalSize[2] +iFieldStart[0]*iFieldGlobalSize[1]*iFieldGlobalSize[2];
            minloc.val = maxloc.val = ( *field )( iifield );
            i_min = iFieldStart[0];
            j_min = iFieldStart[1];
            k_min = iFieldStart[2];
            i_max = iFieldStart[0];
            j_max = iFieldStart[1];
            k_max = iFieldStart[2];
            
            vector<unsigned int> Pcoordinates = patch->Pcoordinates;
            Pcoordinates.resize( 3, 1 );
            
            for( unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++ ) {
                for( unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++ ) {
                    for( unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++ ) {
                        unsigned int ii = k+ ( j + i*iFieldGlobalSize[1] ) *iFieldGlobalSize[2];
                        fieldval = ( *field )( ii );
                        if( minloc.val > fieldval ) {
                            minloc.val = fieldval;
                            i_min=i;
                            j_min=j;
                            k_min=k;
                        }
                        if( maxloc.val < fieldval ) {
                            maxloc.val = fieldval;
                            i_max=i;
                            j_max=j;
                            k_max=k;
                        }
                    }
                }
            }
            
            i_min += Pcoordinates[0]*patch_size_[0] - iFieldStart[0];
            j_min += Pcoordinates[1]*patch_size_[1] - iFieldStart[1];
            k_min += Pcoordinates[2]*patch_size_[2] - iFieldStart[2];
            minloc.index = ( int )( i_min*global_size_[1]*global_size_[2] + j_min*global_size_[2] + k_min );
            i_max += Pcoordinates[0]*patch_size_[0] - iFieldStart[0];
            j_max += Pcoordinates[1]*patch_size_[1] - iFieldStart[1];
            k_max += Pcoordinates[2]*patch_size_[2] - iFieldStart[2];
            maxloc.index = ( int )( i_max*global_size_[1]*global_size_[2] + j_max*global_size_[2] + k_max );
            
            #pragma omp critical
            {
                if( minloc.val < ( double )*fieldMin[ifield] ) {
                    *fieldMin[ifield] = minloc;
                }
                if( maxloc.val > ( double )*fieldMax[ifield] ) {
                    *fieldMax[ifield] = maxloc;
                }
            }
        }
    }
    
    // ------------------------
    // POYNTING-related scalars
    // ------------------------
    
    // electromagnetic energy injected in the simulation (calculated from Poynting fluxes)
    double Uelm_bnd_=0.0;
    unsigned int k=0;
    for( unsigned int j=0; j<2; j++ ) { //directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for( unsigned int i=0; i<EMfields->poynting[j].size(); i++ ) { //axis 0=x, 1=y, 2=z
            if( necessary_poy[k] ) {
                *poy    [k] += EMfields->poynting     [j][i];
                *poyInst[k] += EMfields->poynting_inst[j][i];
                k++;
            }
            
            Uelm_bnd_ += EMfields->poynting[j][i];
        }// i
    }// j
    
    if( necessary_Uelm_BC ) {
        *Uelm_bnd += Uelm_bnd_;
    }
    
} // END compute


double DiagnosticScalar::getScalar( std::string key )
{
    unsigned int k, s=allScalars.size();
    for( k=0; k<s; k++ ) {
        if( allScalars[k]->name_      ==key ) {
            return ( double )*allScalars[k];
        }
        if( allScalars[k]->secondname_==key ) {
            return ( int )   *allScalars[k];
        }
    }
    DEBUG( "key not found " << key );
    return 0.0;
    
} // END getScalar

bool DiagnosticScalar::allowedKey( string key )
{
    if( key.empty() ) {
        return false;
    }
    int s=vars.size();
    if( s==0 ) {
        return true;
    }
    for( int i=0; i<s; i++ ) {
        if( key==vars[i] ) {
            return true;
        }
    }
    return false;
}


bool DiagnosticScalar::needsRhoJs( int itime )
{
    return timeSelection->theTimeIsNow( itime );
}

// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticScalar::getDiskFootPrint( int istart, int istop, Patch *patch )
{
    uint64_t footprint = 0;
    
    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );
    
    if( ndumps == 0 ) {
        return 0;
    }
    
    // Calculate the number of scalars
    // 1 - general scalars
    vector<string> scalars = {"Ubal_norm", "Ubal", "Utot", "Uexp", "Ukin", "Urad", "UmBWpairs", "Uelm", "Ukin_bnd", "Ukin_out_mvw", "Ukin_inj_mvw", "Ukin_new", "Uelm_bnd", "Uelm_out_mvw", "Uelm_inj_mvw"};
    // 2 - species scalars
    unsigned int nspec = patch->vecSpecies.size();
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if( ! patch->vecSpecies[ispec]->particles->is_test ) {
            string species_name = patch->vecSpecies[ispec]->name_;
            scalars.push_back( Tools::merge( "Dens_", species_name ) );
            scalars.push_back( Tools::merge( "Ntot_", species_name ) );
            scalars.push_back( Tools::merge( "Zavg_", species_name ) );
            scalars.push_back( Tools::merge( "Ukin_", species_name ) );
            scalars.push_back( Tools::merge( "Urad_", species_name ) );
        }
    }
    // 3 - Field scalars
    if( !dynamic_cast<ElectroMagnAM *>( patch->EMfields ) ) {
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->Ex_ ->name ) );
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->Ey_ ->name ) );
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->Ez_ ->name ) );
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->Bx_m->name ) );
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->By_m->name ) );
        scalars.push_back( Tools::merge( "Uelm_", patch->EMfields->Bz_m->name ) );
    } else {
        ElectroMagnAM *emfields = static_cast<ElectroMagnAM *>( patch->EMfields );
        unsigned int nmodes = emfields->El_.size();
        for( unsigned int imode=0 ; imode < nmodes ; imode++ ) {
            scalars.push_back( Tools::merge( "Uelm_", emfields->El_[imode] ->name ) );
            scalars.push_back( Tools::merge( "Uelm_", emfields->Er_[imode] ->name ) );
            scalars.push_back( Tools::merge( "Uelm_", emfields->Et_[imode] ->name ) );
            scalars.push_back( Tools::merge( "Uelm_", emfields->Bl_m[imode]->name ) );
            scalars.push_back( Tools::merge( "Uelm_", emfields->Br_m[imode]->name ) );
            scalars.push_back( Tools::merge( "Uelm_", emfields->Bt_m[imode]->name ) );
        }
    }
    // 4 - Scalars related to fields min and max
#ifdef _AM__TODO
    for( unsigned int i=0; i<2; i++ ) {
        string minmax = ( i==0 ) ? "Min" : "Max";
        for( unsigned int j=0; j<2; j++ ) {
            string cell = ( j==0 ) ? "" : "Cell";
            scalars.push_back( Tools::merge( patch->EMfields->Ex_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Ey_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Ez_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Bx_m->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->By_m->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Bz_m->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Jx_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Jy_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->Jz_ ->name, minmax, cell ) );
            scalars.push_back( Tools::merge( patch->EMfields->rho_->name, minmax, cell ) );
            // add envelope-related fields
            if( params.Laser_Envelope_model ) {
                scalars.push_back( Tools::merge( patch->EMfields->Env_A_abs_->name, minmax, cell ) );
                scalars.push_back( Tools::merge( patch->EMfields->Env_Chi_->name, minmax, cell ) );
                scalars.push_back( Tools::merge( patch->EMfields->Env_E_abs_->name, minmax, cell ) );
                scalars.push_back( Tools::merge( patch->EMfields->Env_Ex_abs_->name, minmax, cell ) );
            }
        }
    }
#endif
    // 5 - Scalars related to the Poynting flux
    unsigned int k = 0;
    string poy_name;
    for( unsigned int j=0; j<2; j++ ) {
        for( unsigned int i=0; i<patch->EMfields->poynting[j].size(); i++ ) {
            //if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            //else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            //else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            poy_name = Tools::merge( "Poy", Tools::xyz[i], j==0?"min":"max" );
            scalars.push_back( poy_name );
            scalars.push_back( Tools::merge( poy_name, "Inst" ) );
            k++;
        }
    }
    // Finally, count allowed scalars
    int nscalars = 1;
    for( k=0; k<scalars.size(); k++ ) {
        if( allowedKey( scalars[k] ) ) {
            nscalars ++;
        }
    }
    
    // Calculate the size of one line
    string s = "";
    uint64_t linesize = nscalars * calculateWidth( s );
    
    // Add necessary global headers approximately
    footprint += ( uint64_t )( nscalars * 10 + linesize );
    
    // Add size of each line
    footprint += ndumps * linesize;
    
    return footprint;
}
