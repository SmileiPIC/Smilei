
#include "DiagnosticScalar.h"
#include "VectorPatch.h"

#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

DiagnosticScalar::DiagnosticScalar( Params &params, SmileiMPI* smpi, Patch* patch = NULL ):
latest_timestep(-1)
{
    // patch  == NULL else error
    
    if (PyTools::nComponents("DiagScalar") > 1) {
        ERROR("Only one DiagScalar can be specified");
    }
    
    if (PyTools::nComponents("DiagScalar") > 0 ) {
        
        // get parameter "every" which describes a timestep selection
        timeSelection = new TimeSelection(
            PyTools::extract_py("every", "DiagScalar", 0),
            "Scalars"
        );
        
        precision=10;
        PyTools::extract("precision",precision,"DiagScalar");
        PyTools::extract("vars",vars,"DiagScalar");
        
        // copy from params remaining stuff
        res_time=params.res_time;
        dt=params.timestep;
        cell_volume=params.cell_volume;
    } else {
        timeSelection = new TimeSelection();
    }
    
    // defining default values & reading diagnostic every-parameter
    // ------------------------------------------------------------
    print_every=params.n_time/10;
    PyTools::extract("print_every", print_every, "Main");
    if (!print_every) print_every = 1;
   
    
} // END DiagnosticScalar::DiagnosticScalar


DiagnosticScalar::~DiagnosticScalar()
{
    for(unsigned int i=0; i<allScalars.size(); i++) delete allScalars[i];
    delete timeSelection;
} // END DiagnosticScalar::#DiagnosticScalar


void DiagnosticScalar::openFile( Params& params, SmileiMPI* smpi, bool newfile )
{
    if (!smpi->isMaster()) return;
    
    if (fout.is_open()) return;
    
    //open file scalars.txt
    if ( newfile )
        fout.open("scalars.txt");
    else
        fout.open("scalars.txt", std::ofstream::app);
    
    if (!fout.is_open())
        ERROR("Can't open scalar file");
        
} // END openFile


void DiagnosticScalar::closeFile()
{
    if (fout.is_open()) fout.close();

} // END closeFile


unsigned int DiagnosticScalar::calculateWidth( string key )
{
    return 2 + max(((unsigned int)key.length()),precision+8);
     // The +8 accounts for the dot and exponent in decimal representation)
}

Scalar_value* DiagnosticScalar::newScalar_SUM( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    Scalar_value * scalar = new Scalar_value(name, width, allow, &values_SUM);
    values_SUM.push_back( 0. );
    allScalars.push_back( scalar );
    return scalar;
}
Scalar_value_location* DiagnosticScalar::newScalar_MINLOC( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    val_index default_; default_.index = -1; default_.val = 0.;
    Scalar_value_location * scalar = new Scalar_value_location(name, name+"Cell", width, allow, &values_MINLOC, numeric_limits<double>::max());
    values_MINLOC.push_back( default_ );
    allScalars.push_back( scalar );
    return scalar;
}
Scalar_value_location* DiagnosticScalar::newScalar_MAXLOC( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    val_index default_; default_.index = -1; default_.val = 0.;
    Scalar_value_location * scalar = new Scalar_value_location(name, name+"Cell", width, allow, &values_MAXLOC, numeric_limits<double>::lowest());
    values_MAXLOC.push_back( default_ );
    allScalars.push_back( scalar );
    return scalar;
}


void DiagnosticScalar::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    
    // Make the list of fields
    ElectroMagn* EMfields = vecPatches(0)->EMfields;
    vector<string> fields;
    fields.push_back(EMfields->Ex_ ->name);
    fields.push_back(EMfields->Ey_ ->name);
    fields.push_back(EMfields->Ez_ ->name);
    fields.push_back(EMfields->Bx_m->name);
    fields.push_back(EMfields->By_m->name);
    fields.push_back(EMfields->Bz_m->name);
    fields.push_back(EMfields->Jx_ ->name);
    fields.push_back(EMfields->Jy_ ->name);
    fields.push_back(EMfields->Jz_ ->name);
    fields.push_back(EMfields->rho_->name);
    
    // 1 - Prepare the booleans that tell which scalars are necessary to compute
    // -------------------------------------------------------------------------
    
    // General scalars
    necessary_Ubal_norm = allowedKey("Ubal_norm");
    necessary_Ubal    = necessary_Ubal_norm || allowedKey("Ubal");
    necessary_Utot    = necessary_Ubal || allowedKey("Utot");
    necessary_Uexp    = necessary_Ubal || allowedKey("Uexp");
    necessary_Ukin    = necessary_Utot || allowedKey("Ukin");
    necessary_Uelm    = necessary_Utot || allowedKey("Uelm");
    necessary_Ukin_BC = necessary_Uexp || allowedKey("Ukin_bnd") || allowedKey("Ukin_out_mvw") || allowedKey("Ukin_inj_mvw");
    necessary_Uelm_BC = necessary_Uexp || allowedKey("Uelm_bnd") || allowedKey("Uelm_out_mvw") || allowedKey("Uelm_inj_mvw");
    // Species
    unsigned int nspec = vecPatches(0)->vecSpecies.size();
    necessary_species.resize(nspec, false);
    string species_type;
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if (! vecPatches(0)->vecSpecies[ispec]->particles->isTest) {
            species_type = vecPatches(0)->vecSpecies[ispec]->species_type;
            necessary_species[ispec] = necessary_Ukin || allowedKey("Dens_"+species_type) ||  allowedKey("Ntot_"+species_type) ||  allowedKey("Zavg_"+species_type) || allowedKey("Ukin_"+species_type);
        }
    }
    // Fields
    unsigned int nfield = 6;
    necessary_fieldUelm.resize(nfield, false);
    for( unsigned int ifield=0; ifield<nfield; ifield++ )
        necessary_fieldUelm[ifield] = necessary_Uelm || allowedKey("Uelm_"+fields[ifield]);
    // Fields min/max
    nfield = fields.size();
    necessary_fieldMinMax.resize(nfield, false);
    necessary_fieldMinMax_any = false;
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        necessary_fieldMinMax[ifield] = allowedKey( fields[ifield]+"Min" ) || allowedKey( fields[ifield]+"MinCell" ) || allowedKey( fields[ifield]+"Max" ) || allowedKey( fields[ifield]+"MaxCell" );
        if( necessary_fieldMinMax[ifield] ) necessary_fieldMinMax_any = true;
    }
    // Poynting flux
    unsigned int npoy  = EMfields->poynting[0].size() * EMfields->poynting[1].size();
    necessary_poy.resize(npoy);
    string poy_name;
    unsigned int k = 0;
    for (unsigned int j=0; j<2;j++) {
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {
            if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            necessary_poy[k] = necessary_Uelm_BC || allowedKey(poy_name) || allowedKey(poy_name+"Inst");
            k++;
        }
    }
    
    // 2 - Prepare the Scalar* objects that will contain the data
    // ----------------------------------------------------------
    
    values_SUM   .reserve( 12 + nspec*4 + 6 + 2*npoy);
    values_MINLOC.reserve( 10 );
    values_MAXLOC.reserve( 10 );
    
    // General scalars
    Ubal_norm    = newScalar_SUM( "Ubal_norm"    );
    Ubal         = newScalar_SUM( "Ubal"         );
    Utot         = newScalar_SUM( "Utot"         );
    Uexp         = newScalar_SUM( "Uexp"         );
    Ukin         = newScalar_SUM( "Ukin"         );
    Uelm         = newScalar_SUM( "Uelm"         );
    Ukin_bnd     = newScalar_SUM( "Ukin_bnd"     );
    Ukin_out_mvw = newScalar_SUM( "Ukin_out_mvw" );
    Ukin_inj_mvw = newScalar_SUM( "Ukin_inj_mvw" );
    Uelm_bnd     = newScalar_SUM( "Uelm_bnd"     );
    Uelm_out_mvw = newScalar_SUM( "Uelm_out_mvw" );
    Uelm_inj_mvw = newScalar_SUM( "Uelm_inj_mvw" );
    
    // Scalars related to species
    sDens.resize(nspec, NULL);
    sNtot.resize(nspec, NULL);
    sZavg.resize(nspec, NULL);
    sUkin.resize(nspec, NULL);
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if (! vecPatches(0)->vecSpecies[ispec]->particles->isTest) {
            species_type = vecPatches(0)->vecSpecies[ispec]->species_type;
            sDens[ispec] = newScalar_SUM( "Dens_"+species_type );
            sNtot[ispec] = newScalar_SUM( "Ntot_"+species_type );
            sZavg[ispec] = newScalar_SUM( "Zavg_"+species_type );
            sUkin[ispec] = newScalar_SUM( "Ukin_"+species_type );
        }
    }
    
    // Scalars related to field's electromagnetic energy
    nfield = 6;
    fieldUelm.resize(nfield, NULL);
    for( unsigned int ifield=0; ifield<nfield; ifield++ )
        fieldUelm[ifield] = newScalar_SUM( "Uelm_"+fields[ifield] );
    
    // Scalars related to fields min and max
    nfield = fields.size();
    fieldMin.resize(nfield, NULL);
    fieldMax.resize(nfield, NULL);
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        if( necessary_fieldMinMax[ifield] ) {
            fieldMin[ifield] = newScalar_MINLOC( fields[ifield]+"Min" );
            fieldMax[ifield] = newScalar_MAXLOC( fields[ifield]+"Max" );
        }
    }
    
    // Scalars related to the Poynting flux
    poy    .resize(npoy, NULL);
    poyInst.resize(npoy, NULL);
    k = 0;
    for (unsigned int j=0; j<2;j++) {
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {
            if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            poy    [k] = newScalar_SUM( poy_name        );
            poyInst[k] = newScalar_SUM( poy_name+"Inst" );
            k++;
        }
    }
}



bool DiagnosticScalar::prepare( int timestep )
{
    print_now = printNow(timestep);
    
    // At the right timestep, reset the scalars
    if ( print_now || timeSelection->theTimeIsNow(timestep) )
        for (unsigned int iscalar=0 ; iscalar<allScalars.size() ; iscalar++)
            allScalars[iscalar]->reset();
    
    // Scalars always run even if they don't dump
    return true;
} // END prepare


void DiagnosticScalar::run( Patch* patch, int timestep )
{
    // Must keep track of Poynting flux even without diag
    patch->EMfields->computePoynting(); 
    
    // Compute all scalars when needed
    if ( (print_now || timeSelection->theTimeIsNow(timestep)) && timestep>latest_timestep )
        compute( patch, timestep );

} // END run


void DiagnosticScalar::write(int itime, SmileiMPI* smpi)
{
    if ( smpi->isMaster() ) {
        
        if ( timeSelection->theTimeIsNow(itime) && itime>latest_timestep  ) {
            
            unsigned int j, k, s = allScalars.size();
            
            fout << std::scientific << setprecision(precision);
            // At the beginning of the file, we write some headers
            if (fout.tellp()==ifstream::pos_type(0)) { // file beginning
                // First header: list of scalars, one by line
                fout << "# " << 1 << " time" << endl;
                j = 2;
                for(k=0; k<s; k++) {
                    if( allScalars[k]->allowed ) {
                        fout << "# " << j << " " << allScalars[k]->name << endl;
                        j++;
                        if( ! allScalars[k]->secondname.empty() ) {
                            fout << "# " << j << " " << allScalars[k]->secondname << endl;
                            j++;
                        }
                    }
                }
                // Second header: list of scalars, but all in one line
                fout << "#\n#" << setw(precision+9) << "time";
                for(k=0; k<s; k++) {
                    if( allScalars[k]->allowed ) {
                        fout << setw(allScalars[k]->width) << allScalars[k]->name;
                        if( ! allScalars[k]->secondname.empty() )
                            fout << setw(allScalars[k]->width) << allScalars[k]->secondname;
                    }
                }
                fout << endl;
            }
            // Each requested timestep, the following writes the values of the scalars
            fout << setw(precision+10) << itime/res_time;
            for(k=0; k<s; k++) {
                if( allScalars[k]->allowed ) {
                    fout << setw(allScalars[k]->width) << (double)*allScalars[k];
                    if( ! allScalars[k]->secondname.empty() )
                        fout << setw(allScalars[k]->width) << (int)*static_cast<Scalar_value_location*>(allScalars[k]);
                }
            }
            fout << endl;
        
        }
        
    } // if smpi->isMaster
    
    latest_timestep = itime;
    
} // END write



void DiagnosticScalar::compute( Patch* patch, int timestep )
{
    ElectroMagn* EMfields = patch->EMfields;
    std::vector<Species*>& vecSpecies = patch->vecSpecies;
    
    // ------------------------
    // SPECIES-related energies
    // ------------------------
    double Ukin_=0.;             // total (kinetic) energy carried by particles (test-particles do not contribute)
    double Ukin_bnd_=0.;         // total energy lost by particles due to boundary conditions
    double Ukin_out_mvw_=0.;     // total energy lost due to particles being suppressed by the moving-window
    double Ukin_inj_mvw_=0.;     // total energy added due to particles created by the moving-window
    
    // Compute scalars for each species
    for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
        if (vecSpecies[ispec]->particles->isTest) continue;    // No scalar diagnostic for test particles
        
        if( necessary_species[ispec] || print_now ) {
            double density=0.0;  // sum of weights of current species ispec
            double charge=0.0;   // sum of charges of current species ispec
            double ener_tot=0.0; // total kinetic energy of current species ispec
            
            unsigned int nPart=vecSpecies[ispec]->getNbrOfParticles(); // number of particles
            for (unsigned int iPart=0 ; iPart<nPart; iPart++ ) {
                
                density  += vecSpecies[ispec]->particles->weight(iPart);
                charge   += vecSpecies[ispec]->particles->weight(iPart)
                *          (double)vecSpecies[ispec]->particles->charge(iPart);
                ener_tot += vecSpecies[ispec]->particles->weight(iPart)
                *          (vecSpecies[ispec]->particles->lor_fac(iPart)-1.0);
            }
            ener_tot *= vecSpecies[ispec]->mass;
            
            *sNtot[ispec] += (double)nPart;
            *sDens[ispec] += cell_volume * density;
            *sZavg[ispec] += cell_volume * charge;
            *sUkin[ispec] += cell_volume * ener_tot;
            
            // incremement the total kinetic energy
            Ukin_ += cell_volume * ener_tot;
        }
        
        if( necessary_Ukin_BC || print_now ) {
            // particle energy lost due to boundary conditions
            double ener_lost_bcs=0.0;
            ener_lost_bcs = vecSpecies[ispec]->getLostNrjBC();
            
            // particle energy lost due to moving window
            double ener_lost_mvw=0.0;
            ener_lost_mvw = vecSpecies[ispec]->getLostNrjMW();
            
            // particle energy added due to moving window
            double ener_added_mvw=0.0;
            ener_added_mvw = vecSpecies[ispec]->getNewParticlesNRJ();
            
            // increment all energy loss & energy input
            Ukin_bnd_        += cell_volume * ener_lost_bcs;
            Ukin_out_mvw_    += cell_volume * ener_lost_mvw;
            Ukin_inj_mvw_    += cell_volume * ener_added_mvw;
        }
        
        vecSpecies[ispec]->reinitDiags();
    } // for ispec
    
    // Add the calculated energies to the data arrays
    if( necessary_Ukin || print_now ) {
        *Ukin += Ukin_;
    }
    if( necessary_Ukin_BC || print_now ) {
        *Ukin_bnd     += Ukin_bnd_     ;
        *Ukin_out_mvw += Ukin_out_mvw_ ;
        *Ukin_inj_mvw += Ukin_inj_mvw_ ;
    }
    
    // --------------------------------
    // ELECTROMAGNETIC-related energies
    // --------------------------------
    
    vector<Field*> fields;
    
    fields.push_back(EMfields->Ex_);
    fields.push_back(EMfields->Ey_);
    fields.push_back(EMfields->Ez_);
    fields.push_back(EMfields->Bx_m);
    fields.push_back(EMfields->By_m);
    fields.push_back(EMfields->Bz_m);
    
    double Uelm_=0.0; // total electromagnetic energy in the fields
    
    // loop on all electromagnetic fields
    unsigned int nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        
        if( necessary_fieldUelm[ifield] || print_now ) {
            double Utot_crtField=0.0; // total energy in current field
            Field * field = fields[ifield];
            // compute the starting/ending points of each fields (w/out ghost cells) as well as the field global size
            vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
            for (unsigned int i=0 ; i<field->isDual_.size() ; i++ ) {
                iFieldStart[i]      = EMfields->istart[i][field->isDual(i)];
                iFieldEnd[i]        = iFieldStart[i] + EMfields->bufsize[i][field->isDual(i)];
                iFieldGlobalSize[i] = field->dims_[i];
            }
            
            // loop on all (none-ghost) cells & add-up the squared-field to the energy density
            for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
                for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
                    for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
                        unsigned int ii = k+ (j + i*iFieldGlobalSize[1]) *iFieldGlobalSize[2];
                        Utot_crtField += (*field)(ii) * (*field)(ii);
                    }
                }
            }
            // Utot = Dx^N/2 * Field^2
            Utot_crtField *= 0.5*cell_volume;
            
            *fieldUelm[ifield] += Utot_crtField;
            Uelm_+=Utot_crtField;
        }
    }
    
    // Total elm energy
    if( necessary_Uelm || print_now ) {
        *Uelm += Uelm_;
    }
    
    // Lost/added elm energies through the moving window
    if( necessary_Uelm_BC || print_now ) {
        // nrj lost with moving window (fields)
        double Uelm_out_mvw_ = EMfields->getLostNrjMW();
        Uelm_out_mvw_ *= 0.5*cell_volume;
        
        // nrj added due to moving window (fields)
        double Uelm_inj_mvw_=EMfields->getNewFieldsNRJ();
        Uelm_inj_mvw_ *= 0.5*cell_volume;
        
        *Uelm_out_mvw += Uelm_out_mvw_;
        *Uelm_inj_mvw += Uelm_inj_mvw_ ;
    }
    
    EMfields->reinitDiags();
    
    
    // ---------------------------
    // ELECTROMAGNETIC MIN and MAX
    // ---------------------------
    
    // add currents and density to fields
    fields.push_back(EMfields->Jx_);
    fields.push_back(EMfields->Jy_);
    fields.push_back(EMfields->Jz_);
    fields.push_back(EMfields->rho_);
    
    val_index minloc, maxloc;
    double fieldval;
    
    nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        
        if( necessary_fieldMinMax[ifield] ) {
            
            Field * field = fields[ifield];
            
            vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
            for (unsigned int i=0 ; i<field->isDual_.size() ; i++ ) {
                iFieldStart[i] = EMfields->istart[i][field->isDual(i)];
                iFieldEnd [i] = iFieldStart[i] + EMfields->bufsize[i][field->isDual(i)];
                iFieldGlobalSize [i] = field->dims_[i];
            }
            
            unsigned int ii= iFieldStart[2] + iFieldStart[1]*iFieldGlobalSize[2] +iFieldStart[0]*iFieldGlobalSize[1]*iFieldGlobalSize[2];
            minloc.val = maxloc.val = (*field)(ii);
            minloc.index = maxloc.index = 0;
            
            for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
                for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
                    for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
                        unsigned int ii = k+ (j + i*iFieldGlobalSize[1]) *iFieldGlobalSize[2];
                        fieldval = (*field)(ii);
                        if (minloc.val > fieldval) {
                            minloc.val   = fieldval;
                            minloc.index = ii;
                        }
                        if (maxloc.val < fieldval) {
                            maxloc.val   = fieldval;
                            maxloc.index = ii;
                        }
                    }
                }
            }
            
            #pragma omp critical
            {
                if( minloc.val < (double)*fieldMin[ifield] ) *fieldMin[ifield] = minloc;
                if( maxloc.val > (double)*fieldMax[ifield] ) *fieldMax[ifield] = maxloc;
            }
        }
    }
    
    // ------------------------
    // POYNTING-related scalars
    // ------------------------
    
    // electromagnetic energy injected in the simulation (calculated from Poynting fluxes)
    double Uelm_bnd_=0.0;
    unsigned int k=0;
    for (unsigned int j=0; j<2;j++) {//directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {//axis 0=x, 1=y, 2=z
            if( necessary_poy[k] || print_now ) {
                *poy    [k] += EMfields->poynting     [j][i];
                *poyInst[k] += EMfields->poynting_inst[j][i];
                k++;
            }
            
            Uelm_bnd_ += EMfields->poynting[j][i];
        }// i
    }// j
    
    if( necessary_Uelm_BC || print_now ) {
        *Uelm_bnd += Uelm_bnd_;
    }
    
} // END compute


double DiagnosticScalar::getScalar(std::string key)
{
    unsigned int k, s=allScalars.size();
    for(k=0; k<s; k++) {
        if (allScalars[k]->name      ==key) return (double)*allScalars[k];
        if (allScalars[k]->secondname==key) return (int)   *allScalars[k];
    }
    DEBUG("key not found " << key);
    return 0.0;

} // END getScalar

bool DiagnosticScalar::allowedKey(string key) {
    if( key.empty() ) return false;
    int s=vars.size();
    if (s==0) return true;
    for( int i=0; i<s; i++) {
        if( key==vars[i] ) return true;
    }
    return false;
}


bool DiagnosticScalar::needsRhoJs(int timestep) {
    return printNow(timestep) || timeSelection->theTimeIsNow(timestep);
}
