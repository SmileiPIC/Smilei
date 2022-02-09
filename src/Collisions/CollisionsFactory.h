#ifndef COLLISIONSFACTORY_H
#define COLLISIONSFACTORY_H

#include "Collisions.h"
#include "CollisionsSingle.h"
#include "Params.h"
#include "PyTools.h"

class CollisionsFactory
{
public:

    //! Create one collision object from the input file
    static Collisions *create( Params &params, Patch *patch, std::vector<Species *> vecSpecies, unsigned int n_collisions, bool &debye_length_required )
    {
    
        MESSAGE( 1, "Parameters for collisions #" << n_collisions << " :" );
        
        // Read the input file by searching for the keywords "species1" and "species2"
        // which are the names of the two species that will collide
        std::vector<std::string> sg1( 0 ), sg2( 0 );
        PyTools::extractV( "species1", sg1, "Collisions", n_collisions );
        PyTools::extractV( "species2", sg2, "Collisions", n_collisions );
        
        // Obtain the lists of species numbers from the lists of species names.
        std::vector<std::vector<unsigned int>> sgroup( 2 );
        sgroup[0] = params.FindSpecies( vecSpecies, sg1 );
        sgroup[1] = params.FindSpecies( vecSpecies, sg2 );
        
        // Each group of species sgroup[0] and sgroup[1] must not be empty
        if( sgroup[0].size()==0 ) {
            ERROR( "In collisions #" << n_collisions << ": No valid `species1`" );
        }
        if( sgroup[1].size()==0 ) {
            ERROR( "In collisions #" << n_collisions << ": No valid `species2`" );
        }
        
        // sgroup[0] and sgroup[1] can be equal, but cannot have common species if they are not equal
        bool intra = true;
        if( sgroup[0] != sgroup[1] ) {
            for( unsigned int i0=0; i0<sgroup[0].size(); i0++ ) {
                for( unsigned int i1=0; i1<sgroup[1].size(); i1++ ) {
                    if( sgroup[0][i0] == sgroup[1][i1] ) {
                        ERROR( "In collisions #" << n_collisions << ": species #" << sgroup[0][i0]
                               << " cannot collide with itself" );
                    }
                }
            }
            intra = false;
        }
        
        // Coulomb logarithm (if negative or unset, then automatically computed)
        double clog = 0.; // default
        PyTools::extract( "coulomb_log", clog, "Collisions", n_collisions );
        if( clog <= 0. ) {
            debye_length_required = true;    // auto coulomb log requires debye length
        }
        
        // possibility to multiply the Coulomb by a factor
        double clog_factor = 1.; // default
        PyTools::extract( "coulomb_log_factor", clog_factor, "Collisions", n_collisions );
        if( clog_factor <= 0. ) {
            ERROR_NAMELIST( "In collisions #" << n_collisions << ": coulomb_log_factor must be strictly positive",
            "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
        }
        
        // Number of timesteps between each collisions
        int every = 1; // default
        PyTools::extract( "every", every, "Collisions", n_collisions );
        
        // Number of timesteps between each debug output (if 0 or unset, no debug)
        int debug_every = 0; // default
        PyTools::extract( "debug_every", debug_every, "Collisions", n_collisions );
        
        // Collisional ionization
        int Z = 0; // default
        PyObject * ionizing = PyTools::extract_py( "ionizing", "Collisions", n_collisions );
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
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": ionizing in unknown species `" << ionization_electrons_name << "`",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            if( vecSpecies[ionization_electrons]->atomic_number_ != 0 ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": ionization species are not electrons (atomic_number>0)",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            ionization = true;
        
        } else if( ionizing == Py_True ) {
            
            ionization = true;
        
        } else if( ionizing != Py_False ) {
            
            ERROR_NAMELIST( "In collisions #" << n_collisions << ": `ionizing` must be True, False, or the name of an electron species",
            "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            
        }
        
        CollisionalIonization *Ionization;
        if( ionization ) {
            
            if( intra ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": cannot ionize with intra-collisions",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            
            for( int g=0; g<2; g++ ) { // do sgroup[0], then sgroup[1]
                Species * s0 = vecSpecies[sgroup[g][0]]; // first species of this group
                for( unsigned int i=1; i<sgroup[g].size(); i++ ) { // loop other species of same group
                    Species * s = vecSpecies[sgroup[g][i]]; // current species
                    if( s->mass_ != s0->mass_ )
                        ERROR_NAMELIST( "In collisions #" << n_collisions << ": species in group `species"
                               << g+1 << "` must all have same masses for ionization",
                           "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
                               
                    if( s->atomic_number_ != s0->atomic_number_ ) {
                        if( s->atomic_number_ * s0->atomic_number_ ==0 ) {
                            ERROR_NAMELIST( "In collisions #" << n_collisions << ": species in group `species"
                                   << g+1 << "` cannot be mixed electrons and ions for ionization",
                               "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
                        } else {
                            ERROR_NAMELIST( "In collisions #" << n_collisions << ": species in group `species"
                                   << g+1 << "` must all have same atomic_number for ionization",
                               "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
                        }
                    }
                }
            }
            // atomic number
            unsigned int Z0 = vecSpecies[sgroup[0][0]]->atomic_number_;
            unsigned int Z1 = vecSpecies[sgroup[1][0]]->atomic_number_;
            Z = ( int )( Z0>Z1 ? Z0 : Z1 );
            if( Z0*Z1!=0 ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": ionization requires electrons (no or null atomic_number)",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            if( Z==0 ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": ionization requires ions (atomic_number>0)",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
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
            Ionization = new CollisionalIonization( Z, &params, ionization_electrons, ionization_particles );
            
        } else {
            
            Ionization = new CollisionalNoIonization();
            
        }
        
        // D-D fusion
        PyObject * py_nuclear_reaction = PyTools::extract_py( "nuclear_reaction", "Collisions", n_collisions );
        CollisionalNuclearReaction *NuclearReaction = NULL;
        // If fusion, verify parameters
        if( py_nuclear_reaction == Py_None ) {
            
            NuclearReaction = new CollisionalNoNuclearReaction();
            
        } else {
            
            // Extract the content of the list nuclear_reaction
            std::vector<std::string> nuclear_reaction( 0 );
            if( ! PyTools::py2vector( py_nuclear_reaction, nuclear_reaction ) ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": nuclear_reaction should be a list of strings",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            
            // Verify the atomic number has been set
            unsigned int Z0, Z1;
            if( ! PyTools::extractOrNone( "atomic_number", Z0, "Species", sgroup[0][0] )
             || ! PyTools::extractOrNone( "atomic_number", Z1, "Species", sgroup[1][0] ) ) {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": nuclear_reaction requires all species have an atomic_number",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions");
            }
            
            // Verify each group has consistent atomic number and mass number
            for( int g=0; g<2; g++ ) { // do sgroup[0], then sgroup[1]
                Species * s0 = vecSpecies[sgroup[g][0]]; // first species of this group
                for( unsigned int i=1; i<sgroup[g].size(); i++ ) { // loop other species of same group
                    Species * s = vecSpecies[sgroup[g][i]]; // current species
                    if( s->mass_ != s0->mass_ ) {
                        ERROR_NAMELIST( "In collisions #" << n_collisions << ": nuclear_reaction requires all `species"
                               << g+1 << "` to have equal masses",
                           "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
                    }
                    if( s->atomic_number_ != s0->atomic_number_ ) {
                        ERROR_NAMELIST( "In collisions #" << n_collisions << ": nuclear_reaction requires all `species"
                               << g+1 << "` to have equal atomic_number",
                           "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
                    }
                }
            }
            
            // Rate multiplier
            double rate_multiplier = 0.;
            PyTools::extract( "nuclear_reaction_multiplier", rate_multiplier, "Collisions", n_collisions );
            
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
                findProducts( vecSpecies, products, Z, A, name, product_species, n_collisions );
                
                NuclearReaction = new CollisionalFusionDD( &params, &product_species, rate_multiplier );
                
            // Unknown types
            } else {
                ERROR_NAMELIST( "In collisions #" << n_collisions << ": nuclear_reaction not available between (Z="<<Z0<<",A="<<A0<<") and (Z="<<Z1<<",A="<<A1<<")",
                "https://smileipic.github.io/Smilei/namelist.html#collisions-reactions" );
            }
            
        }
        Py_DECREF( py_nuclear_reaction );
        
        // Print collisions parameters
        std::ostringstream mystream;
        mystream << "(" << sgroup[0][0];
        for( unsigned int rs=1 ; rs<sgroup[0].size() ; rs++ ) {
            mystream << " " << sgroup[0][rs];
        }
        if( intra ) {
            MESSAGE( 2, "Intra collisions within species " << mystream.str() << ")" );
        } else {
            mystream << ") and (" << sgroup[1][0];
            for( unsigned int rs=1 ; rs<sgroup[1].size() ; rs++ ) {
                mystream << " " << sgroup[1][rs];
            }
            MESSAGE( 2, "Collisions between species " << mystream.str() << ")" );
        }
        MESSAGE( 2, "Coulomb logarithm: " << clog );
        if (clog_factor != 1.) {
            MESSAGE( 2, "Coulomb logarithm is multiplied by a factor " << clog_factor );
        }

        if( debug_every>0 ) {
            MESSAGE( 2, "Debug every " << debug_every << " timesteps" );
        }
        mystream.str( "" ); // clear
        if( ionization_electrons>0 ) {
            MESSAGE( 2, "Collisional ionization with atomic number "<<Z<<" towards species `"<<vecSpecies[ionization_electrons]->name_ << "`" );
        }
        if( py_nuclear_reaction != Py_None ) {
            MESSAGE( 2, "Collisional nuclear reaction "<< NuclearReaction->name() );
        }
        
        // If debugging log requested
        std::string filename;
        if( debug_every>0 ) {
            // Build the file name
            mystream.str( "" );
            mystream << "Collisions" << n_collisions << ".h5";
            filename = mystream.str();
            std::ifstream file( filename );
            MESSAGE( "!" << filename << "!");
            // Check if file exists
            if( ! file ) {
                MPI_Comm comm = MPI_COMM_WORLD;
                H5Write f( filename, &comm );
                // write all parameters as HDF5 attributes
                f.attr( "Version", std::string( __VERSION ) );
                mystream.str( "" );
                mystream << sgroup[0][0];
                for( unsigned int i=1; i<sgroup[0].size(); i++ ) {
                    mystream << "," << sgroup[0][i];
                }
                f.attr( "species1", mystream.str() );
                mystream.str( "" );
                mystream << sgroup[1][0];
                for( unsigned int i=1; i<sgroup[1].size(); i++ ) {
                    mystream << "," << sgroup[1][i];
                }
                f.attr( "species2", mystream.str() );
                f.attr( "coulomb_log", clog );
                f.attr( "debug_every", debug_every );
            }
        }
        
        // new Collisions object
        // if( sgroup[0].size()>1 || sgroup[1].size()>1 ) {
            return new Collisions(
                       params,
                       n_collisions,
                       sgroup[0],
                       sgroup[1],
                       clog, clog_factor, intra,
                       every,
                       debug_every,
                       Ionization,
                       NuclearReaction,
                       filename
                   );
        // } else {
        //     return new CollisionsSingle(
        //                params,
        //                n_collisions,
        //                sgroup[0],
        //                sgroup[1],
        //                clog, clog_factor, intra,
        //                every,
        //                debug_every,
        //                Ionization,
        //                NuclearReaction,
        //                filename
        //            );
        // }
    }
    
    
    //! Creates a vector of collisions objects
    static std::vector<Collisions *> create( Params &params, Patch *patch, std::vector<Species *> vecSpecies )
    {
        std::vector<Collisions *> vecCollisions;
        bool debye_length_required = false;
        
        // Needs reference_angular_frequency_SI to be defined
        unsigned int numcollisions=PyTools::nComponents( "Collisions" );
        if( numcollisions > 0 )
            if( params.reference_angular_frequency_SI <= 0. ) {
                ERROR_NAMELIST( "The parameter `reference_angular_frequency_SI` needs to be defined and positive to compute collisions",
                "https://smileipic.github.io/Smilei/namelist.html#reference_angular_frequency_SI");
            }
            
        // Loop over each binary collisions group and parse info
        for( unsigned int n_collisions = 0; n_collisions < numcollisions; n_collisions++ ) {
            vecCollisions.push_back( create( params, patch, vecSpecies, n_collisions, debye_length_required ) );
        }
        for( unsigned int n_collisions = 0; n_collisions < numcollisions; n_collisions++ ) {
            if( vecCollisions[ n_collisions ]->Ionization ) {
                vecCollisions[ n_collisions ]->Ionization->assignDatabase( vecCollisions[ n_collisions ]->Ionization->dataBaseIndex );
            }
        }
        
        // pass the variable "debye_length_required" into the Collision class
        Collisions::debye_length_required = debye_length_required;
        
        return vecCollisions;
    }
    
    
    //! Clone a vector of Collisions objects
    static std::vector<Collisions *> clone( std::vector<Collisions *> vecCollisions )
    {
        std::vector<Collisions *> newVecCollisions( 0 );
        
        for( unsigned int i=0; i<vecCollisions.size(); i++ ) {
            if( dynamic_cast<CollisionsSingle *>( vecCollisions[i] ) ) {
                newVecCollisions.push_back( new CollisionsSingle( vecCollisions[i] ) );
            } else {
                newVecCollisions.push_back( new Collisions( vecCollisions[i] ) );
            }
        }
        
        return newVecCollisions;
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
                        "https://smileipic.github.io/Smilei/namelist.html#nuclear_reaction" );
                    }
                    product_species[i] = s;
                    product_found = true;
                }
            }
            if( ! product_found ) {
                ERROR_NAMELIST( "In collisions #" << n_coll << ", nuclear_reaction : species `"<<s->name_<<"` is not one of "<<list.str(),
                "https://smileipic.github.io/Smilei/namelist.html#nuclear_reaction");
            }
        }
    }
    
};

#endif
