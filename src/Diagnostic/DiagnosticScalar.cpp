#include "DiagnosticScalar.h"

#include <string>
#include <iomanip>

#include "PicParams.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "DiagParams.h"

using namespace std;

// constructor
DiagnosticScalar::DiagnosticScalar(PicParams &params, DiagParams &diagParams, SmileiMPI* smpi) :
isMaster(smpi->isMaster()),
cpuSize(smpi->getSize()),
res_time(params.res_time),
every(diagParams.scalar_every),
cell_volume(params.cell_volume),
precision(diagParams.scalar_precision),
vars(diagParams.scalar_vars)
{
    if (isMaster) {
        fout.open("scalars.txt");
        if (!fout.is_open()) ERROR("can't open scalar file");
    }
}

void DiagnosticScalar::close() {
    if (isMaster) {
        fout.close();
    }
}

// wrapper of the methods
void DiagnosticScalar::run(int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, SmileiMPI *smpi) {
    if (timestep==0) {
        compute(EMfields,vecSpecies,smpi);
        Energy_time_zero  = getScalar("Etot");
        EnergyUsedForNorm = Energy_time_zero;
    }
    if (every) {
        EMfields->computePoynting(); // This must be called everytime        
        if (timestep % every == 0) {
            compute(EMfields,vecSpecies,smpi);
            write(timestep);
        }
    }
}


// it contains all to manage the communication of data. It is "transparent" to the user.
void DiagnosticScalar::compute (ElectroMagn* EMfields, vector<Species*>& vecSpecies, SmileiMPI *smpi) {
    out_list.clear();

    ///////////////////////////////////////////////////////////////////////////////////////////
    // SPECIES STUFF
    ///////////////////////////////////////////////////////////////////////////////////////////
    double Etot_part=0;
    for (unsigned int ispec=0; ispec<vecSpecies.size(); ispec++) {
        double charge_tot=0.0;
        double ener_tot=0.0;
        unsigned int nPart=vecSpecies[ispec]->getNbrOfParticles();

	if (!vecSpecies[ispec]->particles.isTestParticles) {

	  if (nPart>0) {
            for (unsigned int iPart=0 ; iPart<nPart; iPart++ ) {
	      charge_tot+=(double)vecSpecies[ispec]->particles.charge(iPart);
	      ener_tot+=cell_volume*vecSpecies[ispec]->particles.weight(iPart)*(vecSpecies[ispec]->particles.lor_fac(iPart)-1.0);
            }
            ener_tot*=vecSpecies[ispec]->species_param.mass;
	  }

	  MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&charge_tot, &charge_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&ener_tot, &ener_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&nPart, &nPart, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	  if (isMaster) {
            if (nPart!=0) charge_tot /= nPart;
            string nameSpec=vecSpecies[ispec]->species_param.species_type;
            append("Z_"+nameSpec,charge_tot);
            append("E_"+nameSpec,ener_tot);
            append("N_"+nameSpec,nPart);
            Etot_part+=ener_tot;
	  }
	}
    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // ELECTROMAGN STUFF
    ///////////////////////////////////////////////////////////////////////////////////////////


    vector<Field*> fields;

    fields.push_back(EMfields->Ex_);
    fields.push_back(EMfields->Ey_);
    fields.push_back(EMfields->Ez_);
    fields.push_back(EMfields->Bx_m);
    fields.push_back(EMfields->By_m);
    fields.push_back(EMfields->Bz_m);

    double Etot_fields=0.0;

    for (vector<Field*>::iterator field=fields.begin(); field!=fields.end(); field++) {

        map<string,val_index> scalars_map;

        double Etot=0.0;

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
                    Etot+=pow((**field)(ii),2);
                }
            }
        }
        Etot*=0.5*cell_volume;
        MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&Etot, &Etot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if (isMaster) {
            append((*field)->name+"_U",Etot);
            Etot_fields+=Etot;
        }
    }

    // now we add currents and density

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
                        minVal.index=ii; // rank encoded
                    }
                }
            }
        }
        minis.push_back(minVal);
        maxis.push_back(maxVal);
    }

    MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&minis[0], &minis[0], minis.size(), MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
    MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:&maxis[0], &maxis[0], maxis.size(), MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);

    if (isMaster) {
        if (minis.size() == maxis.size() && minis.size() == fields.size()) {
            unsigned int i=0;
            for (vector<Field*>::iterator field=fields.begin(); field!=fields.end() && i<minis.size(); field++, i++) {

                append((*field)->name+"Min",minis[i].val);
                append((*field)->name+"MinCell",minis[i].index);

                append((*field)->name+"Max",maxis[i].val);
                append((*field)->name+"MaxCell",maxis[i].index);

            }
        }

    }

    ///////////////////////////////////////////////////////////////////////////////////////////
    // POYNTING STUFF
    ///////////////////////////////////////////////////////////////////////////////////////////
    double poyTot=0.0;
    for (unsigned int j=0; j<2;j++) {
        for (unsigned int i=0; i<EMfields->poynting[j].size();i++) {

            double poy[2]={EMfields->poynting[j][i],EMfields->poynting_inst[j][i]};

            MPI_Reduce(smpi->isMaster()?MPI_IN_PLACE:poy, poy, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (isMaster) {
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

                poyTot+=poy[0];

            }

            
        }
    }


    ///////////////////////////////////////////////////////////////////////////////////////////
    // FINAL STUFF
    ///////////////////////////////////////////////////////////////////////////////////////////

    if (isMaster) {

        double Total_Energy=Etot_part+Etot_fields;

        double Energy_Balance=Total_Energy-(Energy_time_zero+poyTot);
//        double Energy_Bal_norm=Energy_Balance/Total_Energy;
        double Energy_Bal_norm=Energy_Balance/EnergyUsedForNorm;
        EnergyUsedForNorm = Total_Energy;

        prepend("Poynting",poyTot);
        prepend("EFields",Etot_fields);
        prepend("Eparticles",Etot_part);
        prepend("Etot",Total_Energy);
        prepend("Ebalance",Energy_Balance);
        prepend("Ebal_norm",Energy_Bal_norm);
    }


}

bool DiagnosticScalar::allowedKey(string my_var) {
    bool retval=true;
    if (vars.size()) {
        transform(my_var.begin(), my_var.end(), my_var.begin(), ::tolower);
        vector<string>::const_iterator it = find(vars.begin(), vars.end(),my_var);
        retval=(it != vars.end());
    }    
    return retval;
}

void DiagnosticScalar::prepend(std::string my_var, double val) {
    out_list.insert(out_list.begin(),make_pair(my_var,val));
    
}

void DiagnosticScalar::append(std::string my_var, double val) {
    out_list.push_back(make_pair(my_var,val));
}

void DiagnosticScalar::write(int itime) {
    if(isMaster) {
        fout << std::scientific;
        fout.precision(precision);
        if (fout.tellp()==ifstream::pos_type(0)) {
            fout << "# " << 1 << " time" << endl;
            unsigned int i=2;
            for(vector<pair<string,double> >::iterator iter = out_list.begin(); iter !=out_list.end(); iter++) {
                if (allowedKey((*iter).first) == true) {
                    fout << "# " << i << " " << (*iter).first << endl;
                    i++;
                }
            }

            fout << "#\n#" << setw(precision+9) << "time";
            for(vector<pair<string,double> >::iterator iter = out_list.begin(); iter !=out_list.end(); iter++) {
                if (allowedKey((*iter).first) == true) {
                    fout << setw(precision+9) << (*iter).first;
                }
            }
            fout << endl;
        }
        fout << setw(precision+9) << itime/res_time;
        for(vector<pair<string,double> >::iterator iter = out_list.begin(); iter !=out_list.end(); iter++) {
            if (allowedKey((*iter).first) == true) {
                fout << setw(precision+9) << (*iter).second;
            }
        }
        fout << endl;
    }
}

double DiagnosticScalar::getScalar(string my_var){
    for (unsigned int i=0; i< out_list.size(); i++) {
        if (out_list[i].first==my_var) {
            return out_list[i].second;
        }
    }
    DEBUG("key not found " << my_var);
    return 0.0;
}

