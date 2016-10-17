
#include "DiagnosticScalar.h"

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
        string nameSpec=vecSpecies[ispec]->species_type;
        
        append("Ntot_"+nameSpec,nPart);
        append("Zavg_"+nameSpec,charge_avg);
        append("Ukin_"+nameSpec,ener_tot);
        
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
    for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {
        
        map<string,val_index> scalars_map;
        
        double Utot_crtField=0.0; // total energy in current field
        
        // compute the starting/ending points of each fields (w/out ghost cells) as well as the field global size
        vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
        for (unsigned int i=0 ; i<(*field)->isDual_.size() ; i++ ) {
            iFieldStart[i]      = EMfields->istart[i][(*field)->isDual(i)];
            iFieldEnd[i]        = iFieldStart[i] + EMfields->bufsize[i][(*field)->isDual(i)];
            iFieldGlobalSize[i] = (*field)->dims_[i];
        }
        
        // loop on all (none-ghost) cells & add-up the squared-field to the energy density
        for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
            for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
                for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
                    unsigned int ii=k+ j*iFieldGlobalSize[2] +i*iFieldGlobalSize[1]*iFieldGlobalSize[2];
                    Utot_crtField+=pow((**field)(ii),2);
                }
            }
        }
        // Utot = Dx^N/2 * Field^2
        Utot_crtField *= 0.5*cell_volume;
        
        append("Uelm_"+(*field)->name,Utot_crtField);
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
    
    vector<val_index> minis, maxis;
    
    for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {
        
        val_index minVal, maxVal;
        
        minVal.val=maxVal.val=(**field)(0);
        minVal.index=maxVal.index=0;
        
        vector<unsigned int> iFieldStart(3,0), iFieldEnd(3,1), iFieldGlobalSize(3,1);
        for (unsigned int i=0 ; i<(*field)->isDual_.size() ; i++ ) {
            iFieldStart[i] = EMfields->istart[i][(*field)->isDual(i)];
            iFieldEnd [i] = iFieldStart[i] + EMfields->bufsize[i][(*field)->isDual(i)];
            iFieldGlobalSize [i] = (*field)->dims_[i];
        }
        for (unsigned int k=iFieldStart[2]; k<iFieldEnd[2]; k++) {
            for (unsigned int j=iFieldStart[1]; j<iFieldEnd[1]; j++) {
                for (unsigned int i=iFieldStart[0]; i<iFieldEnd[0]; i++) {
                    unsigned int ii=k+ j*iFieldGlobalSize[2] +i*iFieldGlobalSize[1]*iFieldGlobalSize[2];
                    if (minVal.val>(**field)(ii)) {
                        minVal.val=(**field)(ii);
                        minVal.index=ii; // rank encoded
                    }
                    if (maxVal.val<(**field)(ii)) {
                        maxVal.val=(**field)(ii);
                        maxVal.index=ii; // rank encoded
                    }
                }
            }
        }
        minis.push_back(minVal);
        maxis.push_back(maxVal);
    }

    if (minis.size() == maxis.size() && minis.size() == fields.size()) {
        unsigned int i=0;
        for (vector<Field*>::iterator field=fields.begin(); field!=fields.end() && i<minis.size(); field++, i++) {
            
            append((*field)->name+"Min", minis[i].val, minis[i].index);
            append((*field)->name+"Max", maxis[i].val, maxis[i].index );
            
        }
    }
    
    // ------------------------
    
    // POYNTING-related scalars
    // ------------------------
    
    // electromagnetic energy injected in the simulation (calculated from Poynting fluxes)
    double Uelm_bnd=0.0;
    
    for (unsigned int j=0; j<2;j++) {//directions (xmin/xmax, ymin/ymax, zmin/zmax)
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {//axis 0=x, 1=y, 2=z
            
            double poy[2]={EMfields->poynting[j][i],EMfields->poynting_inst[j][i]};
            
            string name("Poy");
            switch (i) { // dimension
                case 0:
                    name+=(j==0?"East":"West");
                    break;
                case 1:
                    name+=(j==0?"South":"North");
                    break;
                case 2:
                    name+=(j==0?"Bottom":"Top");
                    break;
                default:
                    break;
            }
            append(name,poy[0]);
            append(name+"Inst",poy[1]);
            
            Uelm_bnd += poy[0];
            
        }// i
    }// j
    
    
    // -----------
    // FINAL steps
    // -----------
    
    // added & lost energies due to the moving window
    prepend("Ukin_out_mvw",Ukin_out_mvw);
    prepend("Ukin_inj_mvw",Ukin_inj_mvw);
    prepend("Uelm_out_mvw",Uelm_out_mvw);
    prepend("Uelm_inj_mvw",Uelm_inj_mvw);
    
    // added & lost energies at the boundaries
    prepend("Ukin_bnd",Ukin_bnd);
    prepend("Uelm_bnd",Uelm_bnd);
    
    // Total energies
    prepend("Ukin",Ukin);
    prepend("Uelm",Uelm);
    
    // total energies & energy balance (set later in SmileiMPI::computeGlobalDiags)
    prepend("Ubal_norm",0.);
    prepend("Ubal"     ,0.);
    prepend("Uexp"     ,0.);
    prepend("Utot"     ,0.);
    
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


void DiagnosticScalar::incrementScalar(string my_var, double value){
    for (unsigned int i=0; i< out_key.size(); i++) {
        if ( (out_key[i]==my_var) && ( (my_var.find("Min")== std::string::npos)&&(my_var.find("Max")== std::string::npos)) ) {
            out_value[i] += value;
            return;
        }
    }
    DEBUG("key not found " << my_var);
}


void DiagnosticScalar::incrementScalar(string my_var, double value, int valIndex){
    for (unsigned int i=0; i< out_key.size(); i++) {
        if  ( (out_key[i]==my_var) && ( my_var.find("Min")!= std::string::npos ) && (my_var.find("Cell")== std::string::npos) ) {
            if ( value <= out_value[i] ) {
                // Retain the smallest valIndex when several min are reached.
                if ( value < out_value[i] ) {
                    out_value[i+1] = valIndex;
                } else {
                    out_value[i+1] = std::min(out_value[i+1], (double)valIndex);
                }
                out_value[i] = value;
            }
            return;
        }
        else if  ( (out_key[i]==my_var) && ( my_var.find("Max")!= std::string::npos ) && (my_var.find("Cell")== std::string::npos) ) {
            if ( value >= out_value[i] ) {
                // Retain the smallest valIndex when several max are reached.
                if ( value > out_value[i] ) {
                    out_value[i+1] = valIndex;
                } else {
                    out_value[i+1] = std::min(out_value[i+1], (double)valIndex);
                }
                out_value[i] = value;
            }
            return;
        }
    }
    DEBUG("key not found " << my_var);
}


void DiagnosticScalar::append(std::string key, double value) {
    #pragma omp critical
    {
        if ( !defined(key) ) {
            out_key.push_back(key  );
            out_value.push_back(value);
        }
        else
            incrementScalar(key, value);
    }

}  // END append


void DiagnosticScalar::append(std::string key, double value, int valIndex) {
    #pragma omp critical
    {
        if ( !defined(key) ) {
            out_key.push_back(key  );
            out_value.push_back(value);
            out_key.push_back(key+"Cell");
            out_value.push_back(valIndex);
        }
        else
            incrementScalar(key, value, valIndex);
    }

}  // END append


void DiagnosticScalar::prepend(std::string key, double value) {
    #pragma omp critical
    {
        if ( !defined(key) ) {
            out_key  .insert(out_key  .begin(), key  );
            out_value.insert(out_value.begin(), value);
        }
        else
            incrementScalar(key, value);
    }

} // END prepend


bool DiagnosticScalar::allowedKey(string key) {
    int s=vars.size();
    if (s==0) return true;
    for( int i=0; i<s; i++) {
        if( key==vars[i] ) return true;
    }
    return false;
}

bool DiagnosticScalar::defined(string key) {
    int s=out_key.size();
    for( int i=0; i<s; i++) {
        if( key==out_key[i] ) return true;
    }
    return false;
}

