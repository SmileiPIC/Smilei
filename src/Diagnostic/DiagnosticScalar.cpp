
#include "DiagnosticScalar.h"
#include "VectorPatch.h"

#include <iomanip>
#include <algorithm>
#include <limits>

using namespace std;

DiagnosticScalar::DiagnosticScalar( Params &params, SmileiMPI* smpi, Patch* patch = NULL )
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

Scalar* DiagnosticScalar::newScalar_SUM( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    values_SUM.push_back( 0. );
    Scalar * scalar = new Scalar(name, "", width, allow, &values_SUM.back(), NULL, 0.);
    allScalars.push_back( scalar );
    return scalar;
}
Scalar* DiagnosticScalar::newScalar_MINLOC( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    val_index default_; default_.index = -1; default_.val = 0.;
    values_MINLOC.push_back( default_ );
    Scalar * scalar = new Scalar(name, name+"Cell", width, allow, &(values_MINLOC.back().val), &(values_MINLOC.back().index), numeric_limits<double>::max());
    allScalars.push_back( scalar );
    return scalar;
}
Scalar* DiagnosticScalar::newScalar_MAXLOC( string name )
{
    bool allow = allowedKey(name);
    unsigned int width = calculateWidth( name );
    val_index default_; default_.index = -1; default_.val = 0.;
    values_MAXLOC.push_back( default_ );
    Scalar * scalar = new Scalar(name, name+"Cell", width, allow, &(values_MAXLOC.back().val), &(values_MAXLOC.back().index), numeric_limits<double>::lowest());
    allScalars.push_back( scalar );
    return scalar;
}


void DiagnosticScalar::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{
    
    // Prepare vectors
    ElectroMagn* EMfields = vecPatches(0)->EMfields;
    unsigned int nspec = vecPatches(0)->vecSpecies.size();
    unsigned int npoy  = EMfields->poynting[0].size() * EMfields->poynting[1].size();
    values_SUM   .reserve( 12 + nspec*4 + 6 + 2*npoy);
    values_MINLOC.reserve( 10 );
    values_MAXLOC.reserve( 10 );
    
    // General scalars
    Utot         = newScalar_SUM( "Utot"         );
    Uexp         = newScalar_SUM( "Uexp"         );
    Ubal         = newScalar_SUM( "Ubal"         );
    Ubal_norm    = newScalar_SUM( "Ubal_norm"    );
    Uelm         = newScalar_SUM( "Uelm"         );
    Ukin         = newScalar_SUM( "Ukin"         );
    Uelm_bnd     = newScalar_SUM( "Uelm_bnd"     );
    Ukin_bnd     = newScalar_SUM( "Ukin_bnd"     );
    Ukin_out_mvw = newScalar_SUM( "Ukin_out_mvw" );
    Ukin_inj_mvw = newScalar_SUM( "Ukin_inj_mvw" );
    Uelm_out_mvw = newScalar_SUM( "Uelm_out_mvw" );
    Uelm_inj_mvw = newScalar_SUM( "Uelm_inj_mvw" );
    
    // Scalars related to species
    string species_type;
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
    vector<string> fields;
    fields.push_back(EMfields->Ex_ ->name);
    fields.push_back(EMfields->Ey_ ->name);
    fields.push_back(EMfields->Ez_ ->name);
    fields.push_back(EMfields->Bx_m->name);
    fields.push_back(EMfields->By_m->name);
    fields.push_back(EMfields->Bz_m->name);
    unsigned int nfield = fields.size();
    fieldUelm.resize(nfield);
    for( unsigned int ifield=0; ifield<nfield; ifield++ )
        fieldUelm[ifield] = newScalar_SUM( "Uelm_"+fields[ifield] );
    
    // Scalars related to fields min and max
    fields.push_back(EMfields->Jx_ ->name);
    fields.push_back(EMfields->Jy_ ->name);
    fields.push_back(EMfields->Jz_ ->name);
    fields.push_back(EMfields->rho_->name);
    nfield = fields.size();
    fieldMin.resize(nfield);
    fieldMax.resize(nfield);
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        fieldMin[ifield] = newScalar_MINLOC( fields[ifield]+"Min" );
        fieldMax[ifield] = newScalar_MAXLOC( fields[ifield]+"Max" );
    }
    
    // Scalars related to the Poynting flux
    poy    .resize(npoy);
    poyInst.resize(npoy);
    string poy_name;
    unsigned int k = 0;
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

    // At the right timestep, reset the scalars
    if ( printNow(timestep) || timeSelection->theTimeIsNow(timestep) )
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
    if ( printNow(timestep) || timeSelection->theTimeIsNow(timestep) )
        compute( patch, timestep );

} // END run


void DiagnosticScalar::write(int itime)
{
    unsigned int k, s = allScalars.size();
    
    if ( ! timeSelection->theTimeIsNow(itime) ) return;
    
    fout << std::scientific << setprecision(precision);
    // At the beginning of the file, we write some headers
    if (fout.tellp()==ifstream::pos_type(0)) { // file beginning
        // First header: list of scalars, one by line
        fout << "# " << 1 << " time" << endl;
        for(k=0; k<s; k++) {
            if( allScalars[k]->allowed ) {
                fout << "# " << (k+2) << " " << allScalars[k]->name << endl;
                if( ! allScalars[k]->secondname.empty() )
                    fout << "# " << (k+2) << " " << allScalars[k]->secondname << endl;
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
            fout << setw(allScalars[k]->width) << allScalars[k]->get();
            if( ! allScalars[k]->secondname.empty() )
                fout << setw(allScalars[k]->width) << allScalars[k]->getloc();
        }
    }
    fout << endl;
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
        
        // particle energy lost due to boundary conditions
        double ener_lost_bcs=0.0;
        ener_lost_bcs = vecSpecies[ispec]->getLostNrjBC();
        
        // particle energy lost due to moving window
        double ener_lost_mvw=0.0;
        ener_lost_mvw = vecSpecies[ispec]->getLostNrjMW();
        
        // particle energy added due to moving window
        double ener_added_mvw=0.0;
        ener_added_mvw = vecSpecies[ispec]->getNewParticlesNRJ();
        
        #pragma omp atomic
        *(sNtot[ispec]->value) += nPart;
        #pragma omp atomic
        *(sDens[ispec]->value) += cell_volume * density;
        #pragma omp atomic
        *(sZavg[ispec]->value) += cell_volume * charge;
        #pragma omp atomic
        *(sUkin[ispec]->value) += cell_volume * ener_tot;
        
        // incremement the total kinetic energy
        Ukin_ += cell_volume * ener_tot;
        
        // increment all energy loss & energy input
        Ukin_bnd_        += cell_volume * ener_lost_bcs;
        Ukin_out_mvw_    += cell_volume * ener_lost_mvw;
        Ukin_inj_mvw_    += cell_volume * ener_added_mvw;
        
        vecSpecies[ispec]->reinitDiags();
    } // for ispec
    
    
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
    
    // Compute all electromagnetic energies
    // ------------------------------------
    
    double Uelm_=0.0; // total electromagnetic energy in the fields
    
    // loop on all electromagnetic fields
    unsigned int nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        
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
        
        #pragma omp atomic
        *(fieldUelm[ifield]->value) += Utot_crtField;
        Uelm_+=Utot_crtField;
    }
    
    // nrj lost with moving window (fields)
    double Uelm_out_mvw_ = EMfields->getLostNrjMW();
    Uelm_out_mvw_ *= 0.5*cell_volume;
    
    // nrj added due to moving window (fields)
    double Uelm_inj_mvw_=EMfields->getNewFieldsNRJ();
    Uelm_inj_mvw_ *= 0.5*cell_volume;
    
    EMfields->reinitDiags();
    
    
    // ---------------------------------------------------------------------------------------
    // ALL FIELDS-RELATED SCALARS: Compute all min/max-related scalars (defined on all fields)
    // ---------------------------------------------------------------------------------------
    
    // add currents and density to fields
    fields.push_back(EMfields->Jx_);
    fields.push_back(EMfields->Jy_);
    fields.push_back(EMfields->Jz_);
    fields.push_back(EMfields->rho_);
    
    val_index minloc, maxloc;
    double fieldval;
    
    nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
                
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
            if( minloc.val < *(fieldMin[ifield]->value) ){
                *(fieldMin[ifield]->value) = minloc.val;
                *(fieldMin[ifield]->loc  ) = minloc.index;
            }
            if( maxloc.val > *(fieldMax[ifield]->value) ){
                *(fieldMax[ifield]->value) = maxloc.val;
                *(fieldMax[ifield]->loc  ) = maxloc.index;
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
            #pragma omp atomic
            *(poy    [k]->value) += EMfields->poynting     [j][i];
            #pragma omp atomic
            *(poyInst[k]->value) += EMfields->poynting_inst[j][i];
            k++;
            
            Uelm_bnd_ += EMfields->poynting[j][i];
        }// i
    }// j
    
    
    // -----------
    // FINAL steps
    // -----------
    
    // added & lost energies due to the moving window
    #pragma omp atomic
    *(Ukin_out_mvw->value) += Ukin_out_mvw_;
    #pragma omp atomic
    *(Ukin_inj_mvw->value) += Ukin_inj_mvw_;
    #pragma omp atomic
    *(Uelm_out_mvw->value) += Uelm_out_mvw_;
    #pragma omp atomic
    *(Uelm_inj_mvw->value) += Uelm_inj_mvw_;
    
    // added & lost energies at the boundaries
    #pragma omp atomic
    *(Ukin_bnd->value) += Ukin_bnd_;
    #pragma omp atomic
    *(Uelm_bnd->value) += Uelm_bnd_;
    
    // Total energies
    #pragma omp atomic
    *(Ukin->value) += Ukin_;
    #pragma omp atomic
    *(Uelm->value) += Uelm_;
    
} // END compute


double DiagnosticScalar::getScalar(std::string key)
{
    unsigned int k, s=allScalars.size();
    for(k=0; k<s; k++) {
        if (allScalars[k]->name      ==key) return allScalars[k]->get();
        if (allScalars[k]->secondname==key) return allScalars[k]->getloc();
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
