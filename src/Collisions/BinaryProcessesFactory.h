#ifndef BINARYPROCESSESFACTORY_H
#define BINARYPROCESSESFACTORY_H

#include "BinaryProcesses.h"
#include "BinaryProcess.h"
#include "Collisions.h"
#include "CollisionalIonization.h"
#include "CollisionalFusionDD.h"
#include "Params.h"
#include "PyTools.h"

class BinaryProcessesFactory
{
public:

    //! Create one BinaryProcesses object from the input file
    static BinaryProcesses *create( Params &params, std::vector<Species *> vecSpecies, unsigned int n_binary_processes, bool &debye_length_required )
    {
        
        // Read the input file by searching for the keywords "species1" and "species2"
        // which are the names of the two species that will collide
        std::vector<std::string> sg1( 0 ), sg2( 0 );
        PyTools::extractV( "species1", sg1, "Collisions", n_binary_processes );
        PyTools::extractV( "species2", sg2, "Collisions", n_binary_processes );
        
        // Obtain the lists of species numbers from the lists of species names.
        std::vector<std::vector<unsigned int>> sgroup( 2 );
        sgroup[0] = params.FindSpecies( vecSpecies, sg1 );
        sgroup[1] = params.FindSpecies( vecSpecies, sg2 );
        
        // Each group of species sgroup[0] and sgroup[1] must not be empty
        if( sgroup[0].size()==0 ) {
            ERROR( "In collisions #" << n_binary_processes << ": No valid `species1`" );
        }
        if( sgroup[1].size()==0 ) {
            ERROR( "In collisions #" << n_binary_processes << ": No valid `species2`" );
        }
        
        // sgroup[0] and sgroup[1] can be equal, but cannot have common species if they are not equal
        bool intra = true;
        if( sgroup[0] != sgroup[1] ) {
            for( unsigned int i0=0; i0<sgroup[0].size(); i0++ ) {
                for( unsigned int i1=0; i1<sgroup[1].size(); i1++ ) {
                    if( sgroup[0][i0] == sgroup[1][i1] ) {
                        ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": species #" << sgroup[0][i0]
                               << " cannot collide with itself",
                            LINK_NAMELIST + std::string("#collisions-reactions") );
                    }
                }
            }
            intra = false;
        }
        
        // Number of timesteps between each collisions
        int every = 1; // default
        PyTools::extract( "every", every, "Collisions", n_binary_processes );
        
        // Number of timesteps between each debug output (if 0 or unset, no debug)
        int debug_every = 0; // default
        PyTools::extract( "debug_every", debug_every, "Collisions", n_binary_processes );
        
        // Now make all the binary processes
        std::vector<BinaryProcess*> processes;
        
        // Nuclear reactions
        PyObject * py_nuclear_reaction = PyTools::extract_py( "nuclear_reaction", "Collisions", n_binary_processes );
        std::string nuclear_reaction_name = "";
        // If fusion, verify parameters
        if( py_nuclear_reaction != Py_None ) {
            
            // Extract the content of the list nuclear_reaction
            std::vector<std::string> nuclear_reaction( 0 );
            if( ! PyTools::py2vector( py_nuclear_reaction, nuclear_reaction ) ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": nuclear_reaction should be a list of strings",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            
            // Verify the atomic number has been set
            unsigned int Z0, Z1;
            if( ! PyTools::extractOrNone( "atomic_number", Z0, "Species", sgroup[0][0] )
             || ! PyTools::extractOrNone( "atomic_number", Z1, "Species", sgroup[1][0] ) ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": nuclear_reaction requires all species have an atomic_number",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            
            // Verify each group has consistent atomic number and mass number
            for( int g=0; g<2; g++ ) { // do sgroup[0], then sgroup[1]
                Species * s0 = vecSpecies[sgroup[g][0]]; // first species of this group
                for( unsigned int i=1; i<sgroup[g].size(); i++ ) { // loop other species of same group
                    Species * s = vecSpecies[sgroup[g][i]]; // current species
                    if( s->mass_ != s0->mass_ ) {
                        ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": nuclear_reaction requires all `species"
                               << g+1 << "` to have equal masses",
                               LINK_NAMELIST + std::string("#collisions-reactions") );
                    }
                    if( s->atomic_number_ != s0->atomic_number_ ) {
                        ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": nuclear_reaction requires all `species"
                               << g+1 << "` to have equal atomic_number",
                               LINK_NAMELIST + std::string("#collisions-reactions") );
                    }
                }
            }
            
            // Rate multiplier
            double rate_multiplier = 0.;
            PyTools::extract( "nuclear_reaction_multiplier", rate_multiplier, "Collisions", n_binary_processes );
            
            // Find products
            std::vector<unsigned int> products = params.FindSpecies( vecSpecies, nuclear_reaction );
            
            // Type of reaction
            unsigned int A0 = round( vecSpecies[sgroup[0][0]]->mass_ / 1822.89 );
            unsigned int A1 = round( vecSpecies[sgroup[1][0]]->mass_ / 1822.89 );
            std::vector<Species *> product_species;
            // D-D fusion
            if( Z0 == 1 && Z1 == 1 && A0 == 2 && A1 == 2 ) {
                
                std::vector<unsigned int> Z = {2, 0};
                std::vector<unsigned int> A = {3 ,1};
                std::vector<std::string> name = {"helium3", "neutron"};
                findProducts( vecSpecies, products, Z, A, name, product_species, n_binary_processes );
                
                processes.push_back( new CollisionalFusionDD( params, product_species, rate_multiplier ) );
                
                nuclear_reaction_name = "D-D fusion";
                
            // Unknown types
            } else {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": nuclear_reaction not available between (Z="<<Z0<<",A="<<A0<<") and (Z="<<Z1<<",A="<<A1<<")",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            
        }
        Py_DECREF( py_nuclear_reaction );
        
        
        // Coulomb logarithm (if negative or unset, then no collisions)
        double clog = -1., clog_factor = 1.; // default
        PyTools::extract( "coulomb_log", clog, "Collisions", n_binary_processes );
        if( clog >= 0. ) {
                
            if( clog == 0. ) {
                debye_length_required = true;    // auto coulomb log requires debye length
            }
        
            // possibility to multiply the Coulomb by a factor
            PyTools::extract( "coulomb_log_factor", clog_factor, "Collisions", n_binary_processes );
            if( clog_factor <= 0. ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": coulomb_log_factor must be strictly positive",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            
            processes.push_back( new Collisions( params, clog, clog_factor ) );
        }
        
        // Collisional ionization
        int Z = 0; // default
        PyObject * ionizing = PyTools::extract_py( "ionizing", "Collisions", n_binary_processes );
        bool ionization = false;
        int ionization_electrons = -1;
        Particles * ionization_particles = NULL;
        
        // If `ionizing` is a species name, then use that one
        std::string ionization_electrons_name = "";
        if( PyTools::py2scalar( ionizing, ionization_electrons_name ) ) {
            
            for( int i=0; i<(int)vecSpecies.size(); i++ ) {
                if( vecSpecies[i]->name_ == ionization_electrons_name ) {
                    ionization_electrons = i;
                    break;
                }
            }
            if( ionization_electrons < 0 ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": ionizing in unknown species `" << ionization_electrons_name << "`",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            if( vecSpecies[ionization_electrons]->atomic_number_ != 0 ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": ionization species are not electrons (atomic_number>0)",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            ionization = true;
        
        } else if( ionizing == Py_True ) {
            
            ionization = true;
        
        } else if( ionizing != Py_False ) {
            
            ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": `ionizing` must be True, False, or the name of an electron species",
                LINK_NAMELIST + std::string("#collisions-reactions") );
            
        }
        
        if( ionization ) {
            
            if( intra ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": cannot ionize with intra-collisions",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            
            for( int g=0; g<2; g++ ) { // do sgroup[0], then sgroup[1]
                Species * s0 = vecSpecies[sgroup[g][0]]; // first species of this group
                for( unsigned int i=1; i<sgroup[g].size(); i++ ) { // loop other species of same group
                    Species * s = vecSpecies[sgroup[g][i]]; // current species
                    if( s->mass_ != s0->mass_ )
                        ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": species in group `species"
                               << g+1 << "` must all have same masses for ionization",
                               LINK_NAMELIST + std::string("#collisions-reactions") );
                               
                    if( s->atomic_number_ != s0->atomic_number_ ) {
                        if( s->atomic_number_ * s0->atomic_number_ ==0 ) {
                            ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": species in group `species"
                                   << g+1 << "` cannot be mixed electrons and ions for ionization",
                                   LINK_NAMELIST + std::string("#collisions-reactions") );
                        } else {
                            ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": species in group `species"
                                   << g+1 << "` must all have same atomic_number for ionization",
                                   LINK_NAMELIST + std::string("#collisions-reactions") );
                        }
                    }
                }
            }
            // atomic number
            unsigned int Z0 = vecSpecies[sgroup[0][0]]->atomic_number_;
            unsigned int Z1 = vecSpecies[sgroup[1][0]]->atomic_number_;
            Z = ( int )( Z0>Z1 ? Z0 : Z1 );
            if( Z0*Z1!=0 ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": ionization requires electrons (no or null atomic_number)",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            if( Z==0 ) {
                ERROR_NAMELIST( "In collisions #" << n_binary_processes << ": ionization requires ions (atomic_number>0)",
                    LINK_NAMELIST + std::string("#collisions-reactions") );
            }
            // If ionizing = True, then select ionization electrons as 1st electron species
            if( ionizing == Py_True ) {
                if( Z0==0 ) {
                    ionization_electrons = sgroup[0][0];
                } else if( Z1==0 ) {
                    ionization_electrons = sgroup[1][0];
                }
            }
            ionization_particles = vecSpecies[ionization_electrons]->particles;
            
            // Create the ionization object
            processes.push_back( new CollisionalIonization( Z, &params, ionization_electrons, ionization_particles ) );
            
        }
        
        // Print Binary processes parameters
        std::ostringstream t;
        t << "(" << sgroup[0][0];
        for( unsigned int rs=1 ; rs<sgroup[0].size() ; rs++ ) {
            t << " " << sgroup[0][rs];
        }
        MESSAGE("");
        if( intra ) {
            MESSAGE( 1, "Binary processes #" << n_binary_processes << " within species " << t.str() << ")" );
        } else {
            t << ") and (" << sgroup[1][0];
            for( unsigned int rs=1 ; rs<sgroup[1].size() ; rs++ ) {
                t << " " << sgroup[1][rs];
            }
            MESSAGE( 1, "Binary processes #" << n_binary_processes << " between species " << t.str() << ")" );
        }
        
        for( unsigned int iBP=0; iBP<processes.size(); iBP++ ) {
            MESSAGE( 2, (iBP+1)<<". "<<processes[iBP]->name() );
        }
        
        if( debug_every>0 ) {
            MESSAGE( 2, "Debug every " << debug_every << " timesteps" );
        }
        
        // If debugging log requested
        std::string filename;
        if( debug_every>0 ) {
            // Build the file name
            t.str( "" );
            t << "BinaryProcesses" << n_binary_processes << ".h5";
            filename = t.str();
            std::ifstream file( filename );
            
            // Check if file exists
            if( ! file ) {
                MPI_Comm comm = MPI_COMM_WORLD;
                H5Write f( filename, &comm );
                // write all parameters as HDF5 attributes
                f.attr( "Version", std::string( __VERSION ) );
                t.str( "" );
                t << sgroup[0][0];
                for( unsigned int i=1; i<sgroup[0].size(); i++ ) {
                    t << "," << sgroup[0][i];
                }
                f.attr( "species1", t.str() );
                t.str( "" );
                t << sgroup[1][0];
                for( unsigned int i=1; i<sgroup[1].size(); i++ ) {
                    t << "," << sgroup[1][i];
                }
                f.attr( "species2", t.str() );
                f.attr( "coulomb_log", clog );
                f.attr( "debug_every", debug_every );
            }
        }
        
        return new BinaryProcesses( 
            params,
            sgroup[0],
            sgroup[1],
            intra,
            processes,
            every,
            debug_every,
            filename
        );
    }
    
    
    //! Creates a vector of BinaryProcesses objects
    static std::vector<BinaryProcesses *> createVector( Params &params, std::vector<Species *> vecSpecies )
    {
        std::vector<BinaryProcesses *> vecBPs;
        bool debye_length_required = false;
        
        // Needs reference_angular_frequency_SI to be defined
        unsigned int numcollisions = PyTools::nComponents( "Collisions" );
        if( numcollisions > 0 ) {
            if( params.reference_angular_frequency_SI <= 0. ) {
                ERROR_NAMELIST( "The parameter `reference_angular_frequency_SI` needs to be defined and positive to compute collisions",
                    LINK_NAMELIST + std::string("#reference_angular_frequency_SI") );
            }
        }
        
        // Loop over each binary processes group and parse info
        for( unsigned int n_binary_processes = 0; n_binary_processes < numcollisions; n_binary_processes++ ) {
            vecBPs.push_back( create( params, vecSpecies, n_binary_processes, debye_length_required ) );
        }
        
        // pass the variable "debye_length_required" into the Collision class
        BinaryProcesses::debye_length_required_ = debye_length_required;
        
        return vecBPs;
    }
    
    
    //! Clone a vector of BinaryProcesses objects
    static std::vector<BinaryProcesses *> cloneVector( std::vector<BinaryProcesses *> vecBPs )
    {
        std::vector<BinaryProcesses *> newVecBPs( 0 );
        
        for( unsigned int i=0; i<vecBPs.size(); i++ ) {
            newVecBPs.push_back( new BinaryProcesses( vecBPs[i] ) );
        }
        
        return newVecBPs;
    }
    
    //! Utility for nuclear reactions: find products in the provided list
    static void findProducts( 
        std::vector<Species*> vecSpecies,
        std::vector<unsigned int> products,
        std::vector<unsigned int> Z,
        std::vector<unsigned int> A,
        std::vector<std::string> name,
        std::vector<Species *> &product_species,
        unsigned int n_coll
        )
    {
        unsigned int n = Z.size();
        product_species.resize( n, NULL );
        
        std::ostringstream list("");
        list << name[0];
        for( unsigned int i=1; i<n; i++ ) {
            list << ", " << name[i];
        }
        
        for( unsigned int j=0; j<products.size(); j++ ) {
            Species * s = vecSpecies[products[j]];
            unsigned int Zj = s->atomic_number_;
            unsigned int Aj = round( s->mass_ / 1822.89 );
            
            bool product_found = false;
            for( unsigned int i=0; i<n; i++ ) {
                if( Zj == Z[i] && Aj == A[i] ) {
                    if( product_species[i] ) {
                        ERROR_NAMELIST( "In collisions #" << n_coll << ", nuclear_reaction : should have only 1 "<<name[i]<<" species",
                            LINK_NAMELIST + std::string("#nuclear_reaction") );
                    }
                    product_species[i] = s;
                    product_found = true;
                }
            }
            if( ! product_found ) {
                ERROR_NAMELIST( "In collisions #" << n_coll << ", nuclear_reaction : species `"<<s->name_<<"` is not one of "<<list.str(),
                    LINK_NAMELIST + std::string("nuclear_reaction") );
            }
        }
    }
    
};

#endif
