
#include "DiagnosticScalar.h"
#include "VectorPatch.h"

#include <iomanip>
#include <algorithm>
#include <limits>

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


unsigned int DiagnosticScalar::calculateWidth( string key )
{
    return 2 + max(((unsigned int)key.length()),precision+8);
     // The +8 accounts for the dot and exponent in decimal representation)
}

int DiagnosticScalar::setKey( string key, int &currentIndex )
{
    int actualIndex = -1;
    bool allow = allowedKey(key);
    unsigned int width = calculateWidth( key );
    out_key  .push_back(key  );
    out_value.push_back(0.   );
    out_width.push_back(width);
    allowed  .push_back(allow);
    actualIndex = currentIndex;
    currentIndex++;
    return actualIndex;
}
int DiagnosticScalar::setKey_MINLOC( string key, int &currentIndex )
{
    int actualIndex = -1;
    bool allow = allowedKey(key) || allowedKey(key+"Cell");
    unsigned int width = calculateWidth( key );
    val_index default_; default_.index = -1; default_.val = 0.;
    out_key_MINLOC  .push_back(key     );
    out_value_MINLOC.push_back(default_);
    out_width_MINLOC.push_back(width   );
    allowed_MINLOC  .push_back(allow   );
    actualIndex = currentIndex;
    currentIndex++;
    return actualIndex;
}
int DiagnosticScalar::setKey_MAXLOC( string key, int &currentIndex )
{
    int actualIndex = -1;
    bool allow = allowedKey(key) || allowedKey(key+"Cell");
    unsigned int width = calculateWidth( key );
    val_index default_; default_.index = -1; default_.val = 0.;
    out_key_MAXLOC  .push_back(key     );
    out_value_MAXLOC.push_back(default_);
    out_width_MAXLOC.push_back(width   );
    allowed_MAXLOC  .push_back(allow   );
    actualIndex = currentIndex;
    currentIndex++;
    return actualIndex;
}

void DiagnosticScalar::init(Params& params, SmileiMPI* smpi, VectorPatch& vecPatches)
{

    // Define the indices of each scalar in the output array
    int index = 0, index_MINLOC = 0, index_MAXLOC = 0;
    
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
    index_sDens.resize(nspec);
    index_sNtot.resize(nspec);
    index_sZavg.resize(nspec);
    index_sUkin.resize(nspec);
    for( unsigned int ispec=0; ispec<nspec; ispec++ ) {
        if (vecPatches(0)->vecSpecies[ispec]->particles->isTest) continue;
        species_type = vecPatches(0)->vecSpecies[ispec]->species_type;
        index_sDens[ispec] = setKey( "Dens_"+species_type , index );
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
    index_fieldMin.resize(nfield);
    index_fieldMax.resize(nfield);
    for( unsigned int ifield=0; ifield<nfield; ifield++ ) {
        index_fieldMin[ifield] = setKey_MINLOC( fields[ifield]+"Min", index_MINLOC );
        index_fieldMax[ifield] = setKey_MAXLOC( fields[ifield]+"Max", index_MAXLOC );
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
    // At the right timestep, initialize the scalars
    if ( printNow(timestep) || timeSelection->theTimeIsNow(timestep) ) {
        for (unsigned int iscalar=0 ; iscalar<out_value.size() ; iscalar++) {
            out_value[iscalar] = 0.;
        }
        for (unsigned int iscalar=0 ; iscalar<out_value_MINLOC.size() ; iscalar++) {
            out_value_MINLOC[iscalar].index = -1;
            out_value_MINLOC[iscalar].val   = numeric_limits<double>::max();
        }
        for (unsigned int iscalar=0 ; iscalar<out_value_MAXLOC.size() ; iscalar++) {
            out_value_MAXLOC[iscalar].index = -1;
            out_value_MAXLOC[iscalar].val   = numeric_limits<double>::lowest();
        }
    }
    
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
    unsigned int k;
    unsigned int s        = out_key       .size();
    unsigned int s_MINLOC = out_key_MINLOC.size();
    unsigned int s_MAXLOC = out_key_MAXLOC.size();
    
    if ( ! timeSelection->theTimeIsNow(itime) ) return;
    
    fout << std::scientific << setprecision(precision);
    // At the beginning of the file, we write some headers
    if (fout.tellp()==ifstream::pos_type(0)) { // file beginning
        // First header: list of scalars, one by line
        fout << "# " << 1 << " time" << endl;
        unsigned int i=2;
        for(k=0; k<s; k++) {
            if (allowed[k]) {
                fout << "# " << i << " " << out_key[k] << endl;
                i++;
            }
        }
        for(k=0; k<s_MINLOC; k++) {
            if (allowed_MINLOC[k]) {
                fout << "# " << i << " " << out_key_MINLOC[k] << endl;
                fout << "# " << i << " " << out_key_MINLOC[k]+"Cell" << endl;
                i++;
            }
        }
        for(k=0; k<s_MAXLOC; k++) {
            if (allowed_MAXLOC[k]) {
                fout << "# " << i << " " << out_key_MAXLOC[k] << endl;
                fout << "# " << i << " " << out_key_MAXLOC[k]+"Cell" << endl;
                i++;
            }
        }
        // Second header: list of scalars, but all in one line
        fout << "#\n#" << setw(precision+9) << "time";
        for(k=0; k<s; k++)
            if (allowed[k])
                fout << setw(out_width[k]) << out_key[k];
        for(k=0; k<s_MINLOC; k++) {
            if (allowed_MINLOC[k]) {
                fout << setw(out_width_MINLOC[k]  ) << out_key_MINLOC[k];
                fout << setw(out_width_MINLOC[k]+4) << out_key_MINLOC[k]+"Cell";
            }
        }
        for(k=0; k<s_MAXLOC; k++) {
            if (allowed_MAXLOC[k]) {
                fout << setw(out_width_MAXLOC[k]  ) << out_key_MAXLOC[k];
                fout << setw(out_width_MAXLOC[k]+4) << out_key_MAXLOC[k]+"Cell";
            }
        }
        fout << endl;
    }
    // Each requested timestep, the following writes the values of the scalars
    fout << setw(precision+10) << itime/res_time;
    for(k=0; k<s; k++)
        if (allowed[k])
            fout << setw(out_width[k]) << out_value[k];
    for(k=0; k<s_MINLOC; k++) {
        if (allowed_MINLOC[k]) {
            fout << setw(out_width_MINLOC[k]  ) << out_value_MINLOC[k].val;
            fout << setw(out_width_MINLOC[k]+4) << out_value_MINLOC[k].index;
        }
    }
    for(k=0; k<s_MAXLOC; k++) {
        if (allowed_MAXLOC[k]) {
            fout << setw(out_width_MAXLOC[k]  ) << out_value_MAXLOC[k].val;
            fout << setw(out_width_MAXLOC[k]+4) << out_value_MAXLOC[k].index;
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
        density  *= cell_volume;
        charge   *= cell_volume;
        ener_tot *= cell_volume * vecSpecies[ispec]->mass;
        
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
        out_value[index_sNtot[ispec]] += nPart;
        #pragma omp atomic
        out_value[index_sDens[ispec]] += density;
        #pragma omp atomic
        out_value[index_sZavg[ispec]] += charge;
        #pragma omp atomic
        out_value[index_sUkin[ispec]] += ener_tot;
        
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
            if( minloc.val < out_value_MINLOC[index_fieldMin[ifield]].val )
                out_value_MINLOC[index_fieldMin[ifield]] = minloc;
            if( maxloc.val > out_value_MAXLOC[index_fieldMax[ifield]].val )
                out_value_MAXLOC[index_fieldMax[ifield]] = maxloc;
        }
    }
    
    // ------------------------
    // POYNTING-related scalars
    // ------------------------
    
    // electromagnetic energy injected in the simulation (calculated from Poynting fluxes)
    double Uelm_bnd=0.0;
    unsigned int k=0;
    for (unsigned int j=0; j<2;j++) {//directions (xmin/xmax, ymin/ymax, zmin/zmax)
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
    
    // added & lost energies due to the moving window
    #pragma omp atomic
    out_value[index_Ukin_out_mvw] += Ukin_out_mvw;
    #pragma omp atomic
    out_value[index_Ukin_inj_mvw] += Ukin_inj_mvw;
    #pragma omp atomic
    out_value[index_Uelm_out_mvw] += Uelm_out_mvw;
    #pragma omp atomic
    out_value[index_Uelm_inj_mvw] += Uelm_inj_mvw;
    
    // added & lost energies at the boundaries
    #pragma omp atomic
    out_value[index_Ukin_bnd] += Ukin_bnd;
    #pragma omp atomic
    out_value[index_Uelm_bnd] += Uelm_bnd;
    
    // Total energies
    #pragma omp atomic
    out_value[index_Ukin] += Ukin;
    #pragma omp atomic
    out_value[index_Uelm] += Uelm;
    
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


bool DiagnosticScalar::needsRhoJs(int timestep) {
    return printNow(timestep) || timeSelection->theTimeIsNow(timestep);
}
