#include "Diagnostic.h"

#include <string>
#include <iomanip>

#include <hdf5.h>

#include "Params.h"
#include "SmileiMPI.h"
#include "ElectroMagn.h"
#include "Species.h"

using namespace std;

Diagnostic::Diagnostic(Params& params, vector<Species*>& vecSpecies, SmileiMPI *smpi) :
dtimer(5),
scalars(params,smpi),
probes(params,smpi),
phases(params,smpi)
{
    
    // defining default values & reading diagnostic every-parameter
    // ------------------------------------------------------------
    print_every=params.n_time/10;
    PyTools::extract("print_every", print_every);
    
    if (!PyTools::extract("fieldDump_every", fieldDump_every)) {
        fieldDump_every=params.global_every;
        MESSAGE(1,"Activating all fields to dump");
    }
    
    avgfieldDump_every=params.res_time*10;
    if (!PyTools::extract("avgfieldDump_every", avgfieldDump_every)) avgfieldDump_every=params.global_every;
    
    //!\todo Define default behaviour : 0 or params.res_time
    //ntime_step_avg=params.res_time;
    ntime_step_avg=0;
    PyTools::extract("ntime_step_avg", ntime_step_avg);
    
    particleDump_every=0;
    if (PyTools::extract("particleDump_every", particleDump_every))
        WARNING("Option particleDump_every disabled");
    
    // scalars initialization
    dtimer[0].init(smpi, "scalars");
    
    // probes initialization
    dtimer[1].init(smpi, "probes");
    
    // phasespaces timer initialization
    dtimer[2].init(smpi, "phases");
    
    // particles initialization
    dtimer[3].init(smpi, "particles");
    initParticles(params,vecSpecies);
    
    // test particles initialization
    dtimer[4].init(smpi, "testparticles");
    initTestParticles(params,vecSpecies);
    
}

void Diagnostic::closeAll (SmileiMPI* smpi) {
    
    scalars.closeFile(smpi);
    probes.close();
    phases.close();
    
    for (unsigned int i=0; i<vecDiagnosticParticles.size(); i++) // loop all particle diagnostics
        vecDiagnosticParticles[i]->close();
    
}

void Diagnostic::printTimers (SmileiMPI *smpi, double tottime) {
    
    double coverage(0.);
    if ( smpi->isMaster() ) {
        for (unsigned int i=0 ; i<dtimer.size() ; i++) {
            coverage += dtimer[i].getTime();
        }
    }
    MESSAGE(0, "\n Time in diagnostics : \t"<< tottime <<"\t" << coverage/tottime*100. << "% coverage" );
    if ( smpi->isMaster() ) {
        for (unsigned int i=0 ; i<dtimer.size() ; i++) {
            dtimer[i].print(tottime) ;
        }
    }
}

double Diagnostic::getScalar(string name){
    return scalars.getScalar(name);
}

void Diagnostic::runAllDiags (int timestep, ElectroMagn* EMfields, vector<Species*>& vecSpecies, Interpolator *interp, SmileiMPI *smpi) {
    dtimer[0].restart();
    scalars.run(timestep, EMfields, vecSpecies, smpi);
    dtimer[0].update();
    
    dtimer[1].restart();
    probes.run(timestep, EMfields, interp);
    dtimer[1].update();
    
    dtimer[2].restart();
    phases.run(timestep, vecSpecies);
    dtimer[2].update();
    
    // run all the particle diagnostics
    dtimer[3].restart();
    for (unsigned int i=0; i<vecDiagnosticParticles.size(); i++)
        vecDiagnosticParticles[i]->run(timestep, vecSpecies, smpi);
    dtimer[3].update();
    
    // run all the test particle diagnostics
    dtimer[4].restart();
    for (unsigned int i=0; i<vecDiagnosticTestParticles.size(); i++)
        vecDiagnosticTestParticles[i]->run(timestep, smpi);
    dtimer[4].update();
    
}

void Diagnostic::initParticles(Params& params, vector<Species*> &vecSpecies) {
    unsigned int every, time_average;
    string output;
    vector<string> species;
    vector<unsigned int> species_numbers;
    DiagnosticParticlesAxis  *tmpAxis;
    vector<DiagnosticParticlesAxis*> tmpAxes;
    DiagnosticParticles * tmpDiagParticles;
    vector<PyObject*> allAxes;
    
    bool ok;
    
    unsigned int numDiagParticles=PyTools::nComponents("DiagParticles");
    for (unsigned int n_diag_particles = 0; n_diag_particles < numDiagParticles; n_diag_particles++) {
        
        // get parameter "output" that determines the quantity to sum in the output array
        output = "";
        ok = PyTools::extract("output",output,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `output` required");
        
        // get parameter "every" which is the period (in timesteps) for getting the outputs
        every = 0;
        ok = PyTools::extract("every",every,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `every` required");
        
        // get parameter "time_average" that determines the number of timestep to average the outputs
        time_average = 1;
        PyTools::extract("time_average",time_average,"DiagParticles",n_diag_particles);
        if (time_average > every)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": `time_average` cannot be larger than `every`");
        if (time_average < 1) time_average=1;
        
        // get parameter "species" that determines the species to use (can be a list of species)
        species.clear();
        ok = PyTools::extract("species",species,"DiagParticles",n_diag_particles);
        if (!ok)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": parameter `species` required");
        // verify that the species exist, remove duplicates and sort by number
        species_numbers = params.FindSpecies(vecSpecies, species);
        
        
        // get parameter "axes" that adds axes to the diagnostic
        // Each axis should contain several items:
        //      requested quantity, min value, max value ,number of bins, log (optional), edge_inclusive (optional)
        allAxes=PyTools::extract_pyVec("axes","DiagParticles",n_diag_particles);
        
        if (allAxes.size() == 0)
            ERROR("Diagnotic Particles #" << n_diag_particles << ": axes must contain something");
        
        tmpAxes.clear();
        for (unsigned int iaxis=0; iaxis<allAxes.size(); iaxis++ ) {
            tmpAxis = new DiagnosticParticlesAxis();
            PyObject *oneAxis=allAxes[iaxis];
            if (PyTuple_Check(oneAxis) || PyList_Check(oneAxis)) {
                PyObject* seq = PySequence_Fast(oneAxis, "expected a sequence");
                unsigned int lenAxisArgs=PySequence_Size(seq);
                if (lenAxisArgs<4)
                    ERROR("Diagnotic Particles #" << n_diag_particles << ": axis #" << iaxis << " contain at least 4 arguments");
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 0),tmpAxis->type)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": First item must be a string (axis type)");
                } else {
                    if (   (tmpAxis->type == "z" && params.nDim_particle <3)
                        || (tmpAxis->type == "y" && params.nDim_particle <2) )
                        ERROR("Diagnotic Particles #" << n_diag_particles << ": axis " << tmpAxis->type << " cannot exist in " << params.nDim_particle << "D");
                }
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 1),tmpAxis->min)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Second item must be a double (axis min)");
                }
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 2),tmpAxis->max)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Third item must be a double (axis max)");
                }
                
                
                if (!PyTools::convert(PySequence_Fast_GET_ITEM(seq, 3),tmpAxis->nbins)) {
                    ERROR("Diag Particles #" << n_diag_particles << ", axis #" << iaxis << ": Fourth item must be an int (number of bins)");
                }
                
                // 5 - Check for  other keywords such as "logscale" and "edge_inclusive"
                tmpAxis->logscale = false;
                tmpAxis->edge_inclusive = false;
                for(unsigned int i=4; i<lenAxisArgs; i++) {
                    string my_str("");
                    PyTools::convert(PySequence_Fast_GET_ITEM(seq, i),my_str);
                    if(my_str=="logscale" ||  my_str=="log_scale" || my_str=="log")
                        tmpAxis->logscale = true;
                    else if(my_str=="edges" ||  my_str=="edge" ||  my_str=="edge_inclusive" ||  my_str=="edges_inclusive")
                        tmpAxis->edge_inclusive = true;
                    else
                        ERROR("Diagnotic Particles #" << n_diag_particles << ": keyword `" << my_str << "` not understood");
                }
                
                tmpAxes.push_back(tmpAxis);
                
                Py_DECREF(seq);
            }
            
        }
        // create new diagnostic object
        tmpDiagParticles = new DiagnosticParticles(n_diag_particles, output, every, time_average, species_numbers, tmpAxes);
        // add this object to the list
        vecDiagnosticParticles.push_back(tmpDiagParticles);
    }
}

void Diagnostic::initTestParticles(Params& params, std::vector<Species*>& vecSpecies) {
    DiagnosticTestParticles * tmpDiagTestParticles;
    int n_diag_testparticles=0;
    
    // loop species and make a new diag if test particles
    for(unsigned int i=0; i<vecSpecies.size(); i++) {
        if (vecSpecies[i]->isTest) {
            tmpDiagTestParticles = new DiagnosticTestParticles(n_diag_testparticles, i, params, vecSpecies);
            vecDiagnosticTestParticles.push_back(tmpDiagTestParticles);
            n_diag_testparticles++;
        }
    }
    
}


