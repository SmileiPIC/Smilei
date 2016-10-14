
#include "DiagnosticScalar.h"
#include "VectorPatch.h"

#include <iomanip>
#include <algorithm>

using namespace std;

DiagnosticScalar::DiagnosticScalar( Params &params, SmileiMPI* smpi, Patch* patch = NULL )
{
    // patch  == NULL else error
    
    out_width.resize(0);
    
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


int DiagnosticScalar::setKey( string key, int &currentIndex )
{
    int actualIndex = -1;
    if( allowedKey(key) ) {
        out_key.push_back(key);
        out_value.push_back(0.);
        actualIndex = currentIndex;
        currentIndex++;
    }
    return actualIndex;
}

void DiagnosticScalar::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{

    // Define the indices of each scalar in the output array
    int index = 0;
    
    // General scalars
    index_Utot         = setKey( "Utot"         , index );
    index_Uexp         = setKey( "Uexp"         , index );
    index_Ubal         = setKey( "Ubal"         , index );
    index_Ubal_norm    = setKey( "Ubal_norm"    , index );
    index_Uelm         = setKey( "Uelm"         , index );
    index_Ukin         = setKey( "Ukin"         , index );
    index_Uelm_bnd     = setKey( "Uelm_bnd"     , index );
    index_Ukin_bnd     = setKey( "Ukin_bnd"     , index );
    index_Ukin_out_mvw = setKey( "Ukin_out_mvw" , index );
    index_Ukin_inj_mvw = setKey( "Ukin_inj_mvw" , index );
    index_Uelm_out_mvw = setKey( "Uelm_out_mvw" , index );
    index_Uelm_inj_mvw = setKey( "Uelm_inj_mvw" , index );
    
    // Scalars related to species
    unsigned int nspec = vecPatches(0)->vecSpecies.size();
    string species_type;
    index_sNtot.resize(3*nspec);
    index_sZavg.resize(3*nspec);
    index_sUkin.resize(3*nspec);
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if (vecPatches(0)->vecSpecies[ispec]->particles->isTest) continue;
        species_type = vecPatches(0)->vecSpecies[ispec]->species_type;
        index_sNtot[ispec] = setKey( "Ntot_"+species_type , index );
        index_sZavg[ispec] = setKey( "Zavg_"+species_type , index );
        index_sUkin[ispec] = setKey( "Ukin_"+species_type , index );
    }
    
    // Scalars related to field's electromagnetic energy
    vector<string> fields;
    ElectroMagn* EMfields = vecPatches(0)->EMfields;
    fields.push_back(EMfields->Ex_ ->name);
    fields.push_back(EMfields->Ey_ ->name);
    fields.push_back(EMfields->Ez_ ->name);
    fields.push_back(EMfields->Bx_m->name);
    fields.push_back(EMfields->By_m->name);
    fields.push_back(EMfields->Bz_m->name);
    unsigned int nfield = fields.size();
    index_fieldUelm.resize(nfield);
    for( unsigned int ifield=0; ifield<nfield; ifield++ )
        index_fieldUelm[ifield] = setKey( "Uelm_"+fields[ifield] , index );
    
    // Scalars related to fields min and max
    fields.push_back(EMfields->Jx_ ->name);
    fields.push_back(EMfields->Jy_ ->name);
    fields.push_back(EMfields->Jz_ ->name);
    fields.push_back(EMfields->rho_->name);
    nfield = fields.size();
    index_fieldMin    .resize(nfield);
    index_fieldMinCell.resize(nfield);
    index_fieldMax    .resize(nfield);
    index_fieldMaxCell.resize(nfield);
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        index_fieldMin    [ifield] = setKey( fields[ifield]+"Min"     , index );
        index_fieldMinCell[ifield] = setKey( fields[ifield]+"MinCell" , index );
        index_fieldMax    [ifield] = setKey( fields[ifield]+"Max"     , index );
        index_fieldMaxCell[ifield] = setKey( fields[ifield]+"MaxCell" , index );
    }
    
    // Scalars related to the Poynting flux
    unsigned int k=0;
    index_poy    .resize(0);
    index_poyInst.resize(0);
    string poy_name;
    for (unsigned int j=0; j<2;j++) {
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {
            if     (i==0) poy_name = (j==0?"PoyXmin":"PoyXmax");
            else if(i==1) poy_name = (j==0?"PoyYmin":"PoyYmax");
            else if(i==2) poy_name = (j==0?"PoyZmin":"PoyZmax");
            index_poy    .push_back(-1);
            index_poyInst.push_back(-1);
            index_poy    [k] = setKey( poy_name        , index );
            index_poyInst[k] = setKey( poy_name+"Inst" , index );
            k++;
        }
    }
}



bool DiagnosticScalar::prepare( int timestep )
{
    // At the right timestep, zero-out the scalars
    if ( printNow(timestep) || timeSelection->theTimeIsNow(timestep) )
        for (unsigned int iscalar=0 ; iscalar<out_value.size() ; iscalar++)
            out_value[iscalar] = 0.;
    
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
    unsigned int k, s=out_key.size();
    
    if ( ! timeSelection->theTimeIsNow(itime) ) return;
    
    fout << std::scientific << setprecision(precision);
    // At the beginning of the file, we write some headers
    if (fout.tellp()==ifstream::pos_type(0)) { // file beginning
        // First header: list of scalars, one by line
        fout << "# " << 1 << " time" << endl;
        unsigned int i=2;
        for(k=0; k<s; k++) {
            if (allowedKey(out_key[k])) {
                fout << "# " << i << " " << out_key[k] << endl;
                i++;
            }
        }
        // Second header: list of scalars, but all in one line
        fout << "#\n#" << setw(precision+9) << "time";
        for(k=0; k<s; k++) {
            if (allowedKey(out_key[k])) {
                fout << setw(out_width[k]) << out_key[k];
            }
        }
        fout << endl;
    }
    // Each requested timestep, the following writes the values of the scalars
    fout << setw(precision+10) << itime/res_time;
    for(k=0; k<s; k++) {
        if (allowedKey(out_key[k])) {
            fout << setw(out_width[k]) << out_value[k];
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
    double Ukin=0.;             // total (kinetic) energy carried by particles (test-particles do not contribute)
    double Ukin_bnd=0.;         // total energy lost by particles due to boundary conditions
    double Ukin_out_mvw=0.;     // total energy lost due to particles being suppressed by the moving-window
    double Ukin_inj_mvw=0.;     // total energy added due to particles created by the moving-window
    
    // Compute scalars for each species
    for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
        if (vecSpecies[ispec]->particles->isTest) continue;    // No scalar diagnostic for test particles
        
        double charge_avg=0.0;  // average charge of current species ispec
        double ener_tot=0.0;    // total kinetic energy of current species ispec
        
        unsigned int nPart=vecSpecies[ispec]->getNbrOfParticles(); // number of particles
        if (nPart>0) {
            for (unsigned int iPart=0 ; iPart<nPart; iPart++ ) {
                
                charge_avg += (double)vecSpecies[ispec]->particles->charge(iPart);
                ener_tot   += cell_volume * vecSpecies[ispec]->particles->weight(iPart)
                *             (vecSpecies[ispec]->particles->lor_fac(iPart)-1.0);
            }
            ener_tot*=vecSpecies[ispec]->mass;
        } // if
        
        // particle energy lost due to boundary conditions
        double ener_lost_bcs=0.0;
        ener_lost_bcs = vecSpecies[ispec]->getLostNrjBC();
        
        // particle energy lost due to moving window
        double ener_lost_mvw=0.0;
        ener_lost_mvw = vecSpecies[ispec]->getLostNrjMW();
        
        // particle energy added due to moving window
        double ener_added_mvw=0.0;
        ener_added_mvw = vecSpecies[ispec]->getNewParticlesNRJ();
        
        if (nPart!=0) charge_avg /= nPart;
        
        #pragma omp atomic
        out_value[index_sNtot[ispec]] += nPart     ;
        #pragma omp atomic
        out_value[index_sZavg[ispec]] += charge_avg;
        #pragma omp atomic
        out_value[index_sUkin[ispec]] += ener_tot  ;
        
        // incremement the total kinetic energy
        Ukin += ener_tot;
        
        // increment all energy loss & energy input
        Ukin_bnd        += cell_volume*ener_lost_bcs;
        Ukin_out_mvw    += cell_volume*ener_lost_mvw;
        Ukin_inj_mvw    += cell_volume*ener_added_mvw;
        
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
    
    double Uelm=0.0; // total electromagnetic energy in the fields
    
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
        out_value[index_fieldUelm[ifield]] += Utot_crtField;
        Uelm+=Utot_crtField;
    }
    
    // nrj lost with moving window (fields)
    double Uelm_out_mvw = EMfields->getLostNrjMW();
    Uelm_out_mvw *= 0.5*cell_volume;
    
    // nrj added due to moving window (fields)
    double Uelm_inj_mvw=EMfields->getNewFieldsNRJ();
    Uelm_inj_mvw *= 0.5*cell_volume;
    
    EMfields->reinitDiags();
    
    
    // ---------------------------------------------------------------------------------------
    // ALL FIELDS-RELATED SCALARS: Compute all min/max-related scalars (defined on all fields)
    // ---------------------------------------------------------------------------------------
    
    // add currents and density to fields
    fields.push_back(EMfields->Jx_);
    fields.push_back(EMfields->Jy_);
    fields.push_back(EMfields->Jz_);
    fields.push_back(EMfields->rho_);
    
    double minval, maxval, fieldval;
    unsigned int minindex, maxindex;
    
    nfield = fields.size();
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        
        Field * field = fields[ifield];
        minval=maxval=(*field)(0);
        minindex=maxindex=0;
        
        vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
        for (unsigned int i=0 ; i<field->isDual_.size() ; i++ ) {
            iFieldStart[i] = EMfields->istart[i][field->isDual(i)];
            iFieldEnd [i] = iFieldStart[i] + EMfields->bufsize[i][field->isDual(i)];
            iFieldGlobalSize [i] = field->dims_[i];
        }
        for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
            for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
                for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
                    unsigned int ii = k+ (j + i*iFieldGlobalSize[1]) *iFieldGlobalSize[2];
                    fieldval = (*field)(ii);
                    if (minval>fieldval) {
                        minval=fieldval;
                        minindex=ii; // rank encoded
                    }
                    if (maxval<fieldval) {
                        maxval=fieldval;
                        maxindex=ii; // rank encoded
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            if( minval < out_value[index_fieldMin    [ifield]] ) {
                out_value[index_fieldMin    [ifield]] = minval  ;
                out_value[index_fieldMinCell[ifield]] = minindex;
            }
            if( maxval > out_value[index_fieldMin    [ifield]] ) {
                out_value[index_fieldMax    [ifield]] = maxval  ;
                out_value[index_fieldMaxCell[ifield]] = maxindex;
            }
        }
    }
    
    // ------------------------
    
    // POYNTING-related scalars
    // ------------------------
    
    // electromagnetic energy injected in the simulation (calculated from Poynting fluxes)
    double Uelm_bnd=0.0;
    unsigned int k=0;
    for (unsigned int j=0; j<2;j++) {//directions (west/east, south/north, bottom/top)
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {//axis 0=x, 1=y, 2=z
            #pragma omp atomic
            out_value[index_poy    [k]] += EMfields->poynting     [j][i];
            #pragma omp atomic
            out_value[index_poyInst[k]] += EMfields->poynting_inst[j][i];
            k++;
            
            Uelm_bnd += EMfields->poynting[j][i];
        }// i
    }// j
    
    
    // -----------
    // FINAL steps
    // -----------
    
    //// added & lost energies due to the moving window
    //#pragma omp atomic
    //out_value[index_Ukin_out_mvw] += Ukin_out_mvw;
    //#pragma omp atomic
    //out_value[index_Ukin_inj_mvw] += Ukin_inj_mvw;
    //#pragma omp atomic
    //out_value[index_Uelm_out_mvw] += Uelm_out_mvw;
    //#pragma omp atomic
    //out_value[index_Uelm_inj_mvw] += Uelm_inj_mvw;
    //
    //// added & lost energies at the boundaries
    //#pragma omp atomic
    //out_value[index_Ukin_bnd] += Ukin_bnd;
    //#pragma omp atomic
    //out_value[index_Uelm_bnd] += Uelm_bnd;
    //
    //// Total energies
    //#pragma omp atomic
    //out_value[index_Ukin] += Ukin;
    //#pragma omp atomic
    //out_value[index_Uelm] += Uelm;
    
    // Final thing to do: calculate the maximum size of the scalars names
    if (out_width.empty()) { // Only first time
        unsigned int k, l, s=out_key.size();
        out_width.resize(s);
        for(k=0; k<s; k++) {
            l = out_key[k].length();
            out_width[k] = 2 + max(l,precision+8); // The +8 accounts for the dot and exponent in decimal representation
        }
    }
    //MESSAGE ( "compute : Number of diags = " << out_key.size() ) ;

} // END compute


double DiagnosticScalar::getScalar(std::string key)
{
    unsigned int k, s=out_key.size();
    for(k=0; k<s; k++) {
        if (out_key[k]==key) {
            return out_value[k];
        }
    }
    DEBUG("key not found " << key);
    return 0.0;

} // END getScalar


void DiagnosticScalar::setScalar(string my_var, double value){
    for (unsigned int i=0; i< out_key.size(); i++) {
        if (out_key[i]==my_var) {
          out_value[i] = value;
          return;
        }
    }
    DEBUG("key not found " << my_var);
}

bool DiagnosticScalar::allowedKey(string key) {
    int s=vars.size();
    if (s==0) return true;
    for( int i=0; i<s; i++) {
        if( key==vars[i] ) return true;
    }
    return false;
}
