

#include <sstream>
#include <vector>
#include <limits>

#include "DiagnosticProbes.h"

#include "VectorPatch.h"


using namespace std;





// Calculates the inverse of a square matrix, given row by row
vector<double> matrixInverse( vector<double> A )
{
    unsigned int size = A.size();
    unsigned int dim  = sqrt( size );

    vector<double> invA( size );

    if( dim == 1 ) {
        invA[0] = 1./A[0];
    } else if( dim == 2 ) {
        double det = A[0]*A[3]-A[1]*A[2];
        if( det==0. ) {
            ERROR( "Cannot inverse matrix" );
        }
        double idet = 1./det;
        invA[0] =  A[3]*idet;
        invA[1] = -A[1]*idet;
        invA[2] = -A[2]*idet;
        invA[3] =  A[0]*idet;
    } else if( dim == 3 ) {
        double det = A[0]*A[4]*A[8]+A[1]*A[5]*A[6]+A[2]*A[3]*A[7]-A[2]*A[4]*A[6]-A[1]*A[3]*A[8]-A[0]*A[5]*A[7];
        if( det==0. ) {
            ERROR( "Cannot inverse matrix" );
        }
        double idet = 1./det;
        invA[0] = ( A[4]*A[8]-A[5]*A[7] )*idet;
        invA[1] = ( A[2]*A[7]-A[1]*A[8] )*idet;
        invA[2] = ( A[1]*A[5]-A[2]*A[4] )*idet;
        invA[3] = ( A[5]*A[6]-A[3]*A[8] )*idet;
        invA[4] = ( A[0]*A[8]-A[2]*A[6] )*idet;
        invA[5] = ( A[2]*A[3]-A[0]*A[5] )*idet;
        invA[6] = ( A[3]*A[7]-A[4]*A[6] )*idet;
        invA[7] = ( A[1]*A[6]-A[0]*A[7] )*idet;
        invA[8] = ( A[0]*A[4]-A[1]*A[3] )*idet;
    }

    return invA;
}

// product between a square matrix A and a vector v
vector<double> matrixTimesVector( vector<double> A, vector<double> v )
{
    unsigned int size = A.size();
    unsigned int dim  = sqrt( size );
    if( dim*dim != size ) {
        ERROR( "Matrix is not squared" );
    }
    if( v.size() != dim ) {
        ERROR( "Vector has wrong size" );
    }

    vector<double> w( dim, 0. );
    for( unsigned int i=0; i<dim; i++ ) {
        for( unsigned int j=0; j<dim; j++ ) {
            w[i] += A[i+dim*j] * v[j];
        }
    }
    return w;
}




DiagnosticProbes::DiagnosticProbes( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, int n_probe )
    : Diagnostic( nullptr, "DiagProbe", n_probe )
    , vecNumber( 0 )
    , offset_in_MPI( 0 )
{
    probe_n = n_probe;
    nDim_particle = params.nDim_particle;
    nDim_field = params.nDim_field;
    geometry = params.geometry;
    hasRhoJs = false;
    last_iteration_points_calculated = 0;
    positions_written = false;

    // Extract "every" (time selection)
    ostringstream name( "" );
    name << "Probe #"<<n_probe;
    timeSelection = new TimeSelection(
        PyTools::extract_py( "every", "DiagProbe", n_probe ),
        name.str()
    );

    // Extract "flush_every" (time selection for flushing the file)
    flush_timeSelection = new TimeSelection(
        PyTools::extract_py( "flush_every", "DiagProbe", n_probe ),
        name.str()
    );

    // Extract "number" (number of points you have in each dimension of the probe,
    // which must be smaller than the code dimensions)
    PyTools::extractV( "number", vecNumber, "DiagProbe", n_probe );

    // Dimension of the probe grid
    dimProbe=vecNumber.size();
    if( dimProbe > nDim_particle ) {
        ERROR( "Probe #"<<n_probe<<": probe dimension is greater than simulation dimension" )
    }

    // If there is no "number" argument provided, then it corresponds to
    // a zero-dimensional probe (one point). In this case, we say the probe
    // has actually one dimension with only one point.
    if( dimProbe == 0 ) {
        vecNumber.resize( 1 );
        vecNumber[0]=1;
    }

    // Calculate the total number of points in the grid
    nPart_total=1;
    for( unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++ ) {
        nPart_total *= vecNumber[iDimProbe];
    }

    // Extract origin and corners
    // (positions of the vertices of the grid)
    if( ! PyTools::extractV( "origin", origin, "DiagProbe", n_probe ) ) {
        ERROR( "Probe #"<<n_probe<<": origin missing" );
    }

    if( origin.size() !=nDim_particle ) {
        ERROR( "Probe #"<<n_probe<<": origin size(" << origin.size() << ") != ndim(" << nDim_particle<< ")" );
    }

    corners.resize( 0 );
    bool has_corners = PyTools::extractVV( "corners", corners, "DiagProbe", n_probe );
    bool has_vectors = PyTools::extractVV( "vectors", corners, "DiagProbe", n_probe );
    if( has_corners && has_vectors ) {
        ERROR( "Probe #"<<n_probe<<": cannot define both `corners` and `vectors`" );
    }
    string point_type = has_vectors ? "vectors" : "corners";

    if( corners.size() != dimProbe ) {
        ERROR( "Probe #"<<n_probe<<": needs as many "<<point_type<<" as the size of `number`" );
    }

    for( unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++ ) {
        if( corners[iDimProbe].size() != nDim_particle ) {
            ERROR( "Probe #"<<n_probe<<": "<<point_type<<"["<<iDimProbe<<"] should have "<<nDim_particle<<" elements" );
        }
        if( has_vectors )
            for( unsigned int j=0; j<nDim_particle; j++ ) {
                corners[iDimProbe][j] += origin[j];
            }
    }

    // calculate the coordinate system (base vectors)
    axes.resize( nDim_particle*nDim_particle, 0. );
    vector<int> usedDirections( nDim_particle, 0 );
    double min, max, val;
    unsigned int jmin=0, jmax=0;
    // With available axes, fill the vector, and remember which major direction it occupies
    for( unsigned int i=0; i<dimProbe; i++ ) {
        max = 0.;
        min = numeric_limits<double>::max();
        for( unsigned int j=0; j<nDim_particle; j++ ) {
            axes[j+nDim_particle*i] = corners[i][j] - origin[j];
            val = abs( axes[j+nDim_particle*i] );
            if( val<min ) {
                min=val;
                jmin=j;
            }
            if( val>max ) {
                max=val;
                jmax=j;
            }
        }
        usedDirections[jmax] += 3; // avoid max
        usedDirections[jmin] -= 1; // prefer min
    }
    // Then, complete the probe's coordinate system to have as many axes as the simulation dimension
    for( unsigned int i=dimProbe; i<nDim_particle; i++ ) {
        // find index of the most unused direction
        unsigned int unusedDirectionIndex = min_element( usedDirections.begin(), usedDirections.end() ) - usedDirections.begin();
        // and use that direction as next axis
        axes[unusedDirectionIndex+nDim_particle*i] = 1.;
        for( unsigned int j=0; j<nDim_particle; j++ ) {
            usedDirections[j]--;
        }
        usedDirections[unusedDirectionIndex] += 4;
    }
    // Calculate the inverse matrix of the probe's coordinate system
    axesInverse = matrixInverse( axes );

    // Extract the list of requested fields
    if( !PyTools::extractV( "fields", fieldname, "DiagProbe", n_probe ) ) {
        fieldname.resize( 10 );
        fieldname[0] = "Ex";
        fieldname[1] = "Ey";
        fieldname[2] = "Ez";
        fieldname[3] = "Bx";
        fieldname[4] = "By";
        fieldname[5] = "Bz";
        fieldname[6] = "Jx";
        fieldname[7] = "Jy";
        fieldname[8] = "Jz";
        fieldname[9] = "Rho";
        if( params.Laser_Envelope_model ) {
            fieldname.resize( 14 );
            fieldname[10] = "Env_A_abs";
            fieldname[11] = "Env_Chi";
            fieldname[12] = "Env_E_abs";
            fieldname[13] = "Env_Ex_abs";
        }
    }
    nFields = fieldname.size();
    nBuffers = nFields + 1; // +1 for garbage
    fieldlocation.resize( 17, nFields );
    unsigned int nspec = vecPatches(0)->vecSpecies.size();
    species_field_index.resize( nspec );
    species_field_location.resize( nspec );
    bool has_poynting = false;
    for( unsigned int i=0; i<nFields; i++ ) {
        for( unsigned int j=0; j<i; j++ ) {
            if( fieldname[i]==fieldname[j] ) {
                ERROR( "Probe #"<<n_probe<<": field "<<fieldname[i]<<" appears twice" );
            }
        }
        if( fieldname[i]=="Ex" ) {
            fieldlocation[0] = i;
        } else if( fieldname[i]=="Ey" ) {
            fieldlocation[1] = i;
        } else if( fieldname[i]=="Ez" ) {
            fieldlocation[2] = i;
        } else if( fieldname[i]=="Bx" ) {
            fieldlocation[3] = i;
        } else if( fieldname[i]=="By" ) {
            fieldlocation[4] = i;
        } else if( fieldname[i]=="Bz" ) {
            fieldlocation[5] = i;
        } else if( fieldname[i]=="Jx" ) {
            fieldlocation[6] = i;
            hasRhoJs = true;
        } else if( fieldname[i]=="Jy" ) {
            fieldlocation[7] = i;
            hasRhoJs = true;
        } else if( fieldname[i]=="Jz" ) {
            fieldlocation[8] = i;
            hasRhoJs = true;
        } else if( fieldname[i]=="Rho" ) {
            fieldlocation[9] = i;
            hasRhoJs = true;
        } else if( fieldname[i]=="Env_A_abs" ) {
            fieldlocation[10] = i;
        } else if( fieldname[i]=="Env_Chi" ) {
            fieldlocation[11] = i;
        } else if( fieldname[i]=="Env_E_abs" ) {
            fieldlocation[12] = i;
        } else if( fieldname[i]=="Env_Ex_abs" ) {
            fieldlocation[13] = i;
        } else if( fieldname[i]=="PoyX" ) {
            fieldlocation[14] = i;
            has_poynting = true;
        } else if( fieldname[i]=="PoyY" ) {
            fieldlocation[15] = i;
            has_poynting = true;
        } else if( fieldname[i]=="PoyZ" ) {
            fieldlocation[16] = i;
            has_poynting = true;
        } else {
            // Species-related field
            
            // Find field index
            string f3 = fieldname[i].substr( 0, 3 );
            unsigned int field_index;
            if( fieldname[i].size() < 4 ) {
                ERROR( "Probe #"<<n_probe<<": unknown field `"<<fieldname[i]<<"`" );
            } else if( f3 == "Jx_" ) {
                field_index = 0;
            } else if( f3 == "Jy_" ) {
                field_index = 1;
            } else if( f3 == "Jz_" ) {
                field_index = 2;
            } else if( f3 == "Rho" && fieldname[i][3] == '_' ) {
                field_index = 3;
            } else {
                ERROR( "Probe #"<<n_probe<<": unknown field `"<<fieldname[i]<<"`" );
            }
            
            // Find target species
            size_t i0 = fieldname[i].find( "_" );
            string target_species = fieldname[i].substr( i0+1 );
            
            // Find species number
            unsigned int ispec = 0;
            while( ispec < nspec && target_species != vecPatches( 0 )->vecSpecies[ispec]->name_ ) {
                ispec++;
            }
            if( ispec == nspec ) {
                ERROR( "Probe #"<<n_probe<<": unknown field `"<<fieldname[i]<<"` (unknown species `"<<target_species<<"`)" );
            }
            
            // Allocate fields
            unsigned int start = vecPatches( 0 )->EMfields->species_starts[ispec];
            unsigned int stop = vecPatches( 0 )->EMfields->species_starts[ispec+1];
            if( geometry == "AMcylindrical" ) { // allocate all fields/modes in AM
                if( species_field_index[ispec].size() == 0 ) { // only once
                    for( unsigned int ifield=start; ifield<stop; ifield++ ) {
                        vecPatches.allocateField( ifield, params );
                    }
                }
            } else {
                vecPatches.allocateField( start + field_index, params );
            }
            // Store the input and output locations
            species_field_index[ispec].push_back( field_index );
            species_field_location[ispec].push_back( i );
            hasRhoJs = true;
        }
    }
    
    // If Poynting requested, force temporary storage of fields
    if( has_poynting ) {
        for( unsigned int k=0; k<6; k++ ) {
            if( fieldlocation[k] >= nFields ) {
                fieldlocation[k] = nBuffers;
                nBuffers++;
            }
        }
    }
    
    // Extract time_integral
    PyTools::extract( "time_integral", time_integral, "DiagProbe", n_probe );
    if( time_integral && params.hasWindow ) {
        ERROR( "Probe #"<<n_probe<<": `time_integral` incompatible with the moving window" );
    }
    
    // Pre-calculate patch size
    patch_length.resize( nDim_particle );
    for( unsigned int k=0; k<nDim_particle; k++ ) {
        patch_length[k] = params.patch_size_[k]*params.cell_length[k];
    }

    // Create filename
    ostringstream mystream( "" );
    mystream << "Probes" << n_probe << ".h5";
    filename = mystream.str();
    
    // Create points initially
    if( ! params.restart ) {
        createPoints( smpi, vecPatches, 0. );
    }
    
    // Display info
    MESSAGE( 1, "Probe diagnostic #"<<n_probe<< (time_integral?" (integrated over time)":"") );
    
    ostringstream t( "" );
    t << vecNumber[0];
    for( unsigned int i=1; i<dimProbe; i++ ) {
        t << "x"<<vecNumber[i];
    }
    t << " points";
    if( dimProbe>1 ) {
        t << " (total = "<<nPart_total<<")";
    }
    MESSAGE( 2, t.str() );

    t.str( "" );
    t << "origin : "<<origin[0];
    for( unsigned int k=1; k<nDim_particle; k++ ) {
        t<<", "<<origin[k];
    }
    MESSAGE( 2, t.str() );

    for( unsigned int i=0; i<corners.size(); i++ ) {
        t.str( "" );
        t << "corner "<<i<<" : "<<corners[i][0];
        for( unsigned int k=1; k<nDim_particle; k++ ) {
            t<<", "<<corners[i][k];
        }
        MESSAGE( 2, t.str() );
    }

} // END DiagnosticProbes::DiagnosticProbes



DiagnosticProbes::~DiagnosticProbes()
{
    delete timeSelection;
    delete flush_timeSelection;
}


void DiagnosticProbes::openFile( Params &, SmileiMPI *smpi )
{
    file_ = new H5Write( filename, &smpi->world() );
    
    file_->attr( "name", diag_name_ );
    file_->attr( "Version", string( __VERSION ) );
    file_->attr( "dimension", dimProbe );
    file_->attr( "time_integral", time_integral );
    
    // Add arrays "p0", "p1", ...
    file_->vect( "p0", origin );
    ostringstream pk;
    for( unsigned int iDimProbe=0; iDimProbe<dimProbe; iDimProbe++ ) {
        pk.str( "" );
        pk << "p" << ( iDimProbe+1 );
        file_->vect( pk.str(), corners[iDimProbe] );
    }
    
    // Add array "number"
    file_->vect( "number", vecNumber );
    
    // Add "fields"
    ostringstream fields( "" );
    fields << fieldname[0];
    for( unsigned int i=1; i<fieldname.size(); i++ ) {
        fields << "," << fieldname[i];
    }
    file_->attr( "fields", fields.str() );
    
}


void DiagnosticProbes::closeFile()
{
    if( file_ ) {
        delete file_;
        file_ = NULL;
    }
}


bool DiagnosticProbes::prepare( int itime )
{
    return time_integral || timeSelection->theTimeIsNow( itime );
} 


void DiagnosticProbes::init( Params &params, SmileiMPI *smpi, VectorPatch & )
{
    // create the file
    openFile( params, smpi );
}


void DiagnosticProbes::createPoints( SmileiMPI *smpi, VectorPatch &vecPatches, double x_moved )
{
    if( smpi->test_mode ) return;
    
    nPart_MPI = 0;
    offset_in_MPI .resize( vecPatches.size() );
    offset_in_file.resize( vecPatches.size() );
    unsigned int numCorners = 1<<nDim_particle; // number of patch corners
    unsigned int ntot, IP, ipart_local, i, k, iDim;
    bool is_in_domain;
    vector<double> point( nDim_particle ), mins( nDim_particle ), maxs( nDim_particle ); //warning, works only if nDim_particle >= nDim_field
    vector<double> patchMin( nDim_particle ), patchMax( nDim_particle );
    vector<unsigned int> minI( nDim_particle ), maxI( nDim_particle ), nI( nDim_particle );

    // Loop patches to create particles
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {

        // The first step is to reduce the area of the probe to search in this patch
        if( geometry == "AMcylindrical" ) {
            mins[0] = numeric_limits<double>::max();
            maxs[0] = numeric_limits<double>::lowest();
            patchMin[0] = ( vecPatches( ipatch )->Pcoordinates[0] )*patch_length[0];
            patchMax[0] = ( vecPatches( ipatch )->Pcoordinates[0]+1 )*patch_length[0];
            for( k=1; k<3; k++ ) {
                mins[k] = numeric_limits<double>::max();
                maxs[k] = -mins[k];
                patchMax[k] = ( vecPatches( ipatch )->Pcoordinates[1]+1 )*patch_length[1];
                patchMin[k] = - patchMax[k] ; //patchMin = -rmax for the first filter
            }
        } else {
            for( k=0; k<nDim_particle; k++ ) {
                mins[k] = numeric_limits<double>::max();
                maxs[k] =  numeric_limits<double>::lowest();
                patchMin[k] = ( vecPatches( ipatch )->Pcoordinates[k] )*patch_length[k];
                patchMax[k] = ( vecPatches( ipatch )->Pcoordinates[k]+1 )*patch_length[k];
            }
        }
        // loop patch corners
        for( i=0; i<numCorners; i++ ) {
            // Get coordinates of the current corner in terms of x,y,...
            for( k=0; k<nDim_particle; k++ ) {
                point[k] = ( ( ( ( i>>k )&1 )==0 ) ? patchMin[k] : patchMax[k] ) - origin[k];
            }
            // Get position of the current corner in the probe's coordinate system
            point = matrixTimesVector( axesInverse, point );
            // Store mins and maxs
            for( k=0; k<nDim_particle; k++ ) {
                if( point[k]<mins[k] ) {
                    mins[k]=point[k];
                }
                if( point[k]>maxs[k] ) {
                    maxs[k]=point[k];
                }
            }
        }
        // Loop directions to figure out the range of useful indices
        for( i=0; i<dimProbe; i++ ) {
            if( mins[i]<0. ) {
                mins[i]=0.;
            }
            if( mins[i]>1. ) {
                mins[i]=1.;
            }
            minI[i] = ( ( unsigned int ) floor( mins[i]*( ( double )( vecNumber[i]-1 ) ) ) );
            if( maxs[i]<0. ) {
                maxs[i]=0.;
            }
            if( maxs[i]>1. ) {
                maxs[i]=1.;
            }
            maxI[i] = ( ( unsigned int ) ceil( maxs[i]*( ( double )( vecNumber[i]-1 ) ) ) ) + 1;
        }
        for( i=dimProbe; i<nDim_particle; i++ ) {
            minI[i] = 0;
            maxI[i] = ( mins[i]*maxs[i]<=1.e-8 ) ? 1:0;
        }
        // Now, minI and maxI contain the min and max indexes of the probe, useful for this patch
        // Calculate total number of useful points
        ntot = 1;
        for( i=0; i<nDim_particle; i++ ) {
            nI[i] = maxI[i]-minI[i];
            ntot *= nI[i];
        }
        if( ntot > 1000000000 ) {
            ERROR( "Probe too large" );
        }
        // Initialize the list of "fake" particles (points) just as actual macro-particles
        Particles *particles = &( vecPatches( ipatch )->probes[probe_n]->particles );
        particles->initialize( ntot, nDim_particle, false );
        // In AM, redefine patchmin as rmin and not -rmax anymore
        if( geometry == "AMcylindrical" ) {
            patchMin[1] = patchMax[1] - patch_length[1];
        }
        // Loop useful probe points
        ipart_local=0;
        for( unsigned int ip=0; ip<ntot; ip++ ) {
            // Find the coordinates of this point in the global probe array
            IP = ip;
            for( i=0; i<dimProbe; i++ ) {
                point[i] = ( ( double )( IP % nI[i] + minI[i] ) ) / ( ( double )( vecNumber[i]-1 ) );
                IP /= nI[i];
            }
            for( i=dimProbe; i<nDim_particle; i++ ) {
                point[i] = 0.;
            }
            // Compute this point's coordinates in terms of x, y, ...
            point = matrixTimesVector( axes, point );
            // Check if point is in patch
            is_in_domain = true;
            if( geometry == "AMcylindrical" ) {
                point[0] += origin[0];
                point[1] += origin[1];
                point[2] += origin[2];
                double r = sqrt( point[1]*point[1] + point[2]*point[2] );
                if( point[0] < patchMin[0] || point[0] >= patchMax[0]
                 || r < patchMin[1] || r>= patchMax[1] ) {
                    is_in_domain = false;
                }
            } else {
                for( i=0; i<nDim_field; i++ ) {
                    point[i] += origin[i];
                    if( point[i] < patchMin[i] || point[i] >= patchMax[i] ) {
                        is_in_domain = false;
                        break;
                    }
                }
            }
            if( is_in_domain ) {
                point[0] += x_moved;
                for( iDim=0; iDim<nDim_particle; iDim++ ) {
                    particles->position( iDim, ipart_local ) = point[iDim];
                }
                ipart_local++;
            }
        }

        // Resize the array with only particles in this patch
        particles->resize( ipart_local, nDim_particle, false );
        particles->shrinkToFit();

        // Add the local offset
        offset_in_MPI[ipatch] = nPart_MPI;
        nPart_MPI += ipart_local;
    }

    // Calculate the offset of each MPI in the final file
    unsigned int global_offset;
    MPI_Scan( &nPart_MPI, &global_offset, 1, MPI_UNSIGNED, MPI_SUM, MPI_COMM_WORLD );

    // Broadcast the global number of points
    nPart_total_actual = global_offset;
    MPI_Bcast( &nPart_total_actual, 1, MPI_UNSIGNED, smpi->getSize()-1, MPI_COMM_WORLD );

    // For each patch, calculate its global offset
    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
        offset_in_file[ipatch] = global_offset - nPart_MPI + offset_in_MPI[ipatch];
        vecPatches( ipatch )->probes[probe_n]->offset_in_file = offset_in_file[ipatch];
    }

    if( nPart_total_actual==0 ) {
        ERROR( "Probe has no points in the box" );
    }
}



void DiagnosticProbes::run( SmileiMPI *smpi, VectorPatch &vecPatches, int itime, SimWindow *simWindow, Timers & )
{
    ostringstream name_t;
    
    unsigned int nPatches( vecPatches.size() );
    double x_moved = simWindow ? simWindow->getXmoved() : 0.;
    
    // Leave if this timestep has already been written
    #pragma omp master
    {
        name_t.str( "" );
        name_t << "/" << setfill( '0' ) << setw( 10 ) << itime;
        dataset_name = name_t.str();
        has_dataset = file_->has( dataset_name );
    }
    #pragma omp barrier
    if( has_dataset ) {
        return;
    }
    
    #pragma omp master
    {
        // If the patches have been moved (moving window or load balancing) we must re-compute the probes positions
        if( !positions_written || last_iteration_points_calculated <= vecPatches.lastIterationPatchesMoved ) {
            createPoints( smpi, vecPatches, x_moved );
            last_iteration_points_calculated = itime;

            // Store the positions of all particles, unless done already
            if( !positions_written ) {
                vector<unsigned int> posArraySize( 2 );
                posArraySize[0] = nPart_MPI;
                posArraySize[1] = nDim_particle;
                Field2D *posArray;
                if( nPart_MPI > 0 ) {
                    posArray = new Field2D( posArraySize );
                    unsigned int ipart = 0;
                    for( unsigned int ipatch=0 ; ipatch<vecPatches.size() ; ipatch++ ) {
                        if( ipart>=nPart_MPI ) {
                            break;
                        }
                        Particles *particles = &( vecPatches( ipatch )->probes[probe_n]->particles );
                        for( unsigned int ip=0 ; ip<particles->hostVectorSize() ; ip++ ) {
                            for( unsigned int idim=0 ; idim<nDim_particle  ; idim++ ) {
                                ( *posArray )( ipart, idim ) = particles->position( idim, ip );
                            }
                            ( *posArray )( ipart, 0 ) -= x_moved;
                            ipart++;
                        }
                    }
                } else {
                    posArray = new Field2D();
                }
                
                // Define spaces
                H5Space memspace( {nPart_MPI, nDim_particle}, {}, {} );
                H5Space filespace( {nPart_total_actual, nDim_particle}, {offset_in_file[0], 0}, {nPart_MPI, nDim_particle} );
                // Create dataset
                file_->array( "positions", *(posArray->data_), &filespace, &memspace );
                file_->flush();
                
                delete posArray;
                positions_written = true;
            }
        }
        
        // Make the array that will contain the data
        vector<unsigned int> probesArraySize( 2 );
        probesArraySize[1] = nPart_MPI; // number of particles
        probesArraySize[0] = nBuffers;
        probesArray = new Field2D( probesArraySize );
    }
    #pragma omp barrier

    // Loop patches to fill the array
    #pragma omp for schedule(runtime)
    for( unsigned int ipatch=0 ; ipatch<nPatches ; ipatch++ ) {
        Patch * patch = vecPatches( ipatch );
        unsigned int npart = patch->probes[probe_n]->particles.hostVectorSize();
        
        LocalFields Jloc_fields;
        double Rloc_fields;
        
        int ithread = 0;
#ifdef _OPENMP
        ithread = omp_get_thread_num();
#endif
        
        // Interpolate all usual fields on probe ("fake") particles of current patch
        unsigned int iPart_MPI = offset_in_MPI[ipatch];
        unsigned int maxPart_MPI = offset_in_MPI[ipatch] + npart;
        smpi->resizeBuffers( ithread, nDim_particle, npart, false );
        for( unsigned int ipart=0; ipart<npart; ipart++ ) {
            int iparticle( ipart ); // Compatibility
            int false_idx( 0 );   // Use in classical interp for now, not for probes
            patch->probesInterp->fieldsAndCurrents(
                patch->EMfields,
                patch->probes[probe_n]->particles, smpi,
                &iparticle, &false_idx, ithread,
                &Jloc_fields, &Rloc_fields
            );
            //! here we fill the probe data!!!
            ( *probesArray )( fieldlocation[0], iPart_MPI )=smpi->dynamics_Epart[ithread][ipart+0*npart];
            ( *probesArray )( fieldlocation[1], iPart_MPI )=smpi->dynamics_Epart[ithread][ipart+1*npart];
            ( *probesArray )( fieldlocation[2], iPart_MPI )=smpi->dynamics_Epart[ithread][ipart+2*npart];
            ( *probesArray )( fieldlocation[3], iPart_MPI )=smpi->dynamics_Bpart[ithread][ipart+0*npart];
            ( *probesArray )( fieldlocation[4], iPart_MPI )=smpi->dynamics_Bpart[ithread][ipart+1*npart];
            ( *probesArray )( fieldlocation[5], iPart_MPI )=smpi->dynamics_Bpart[ithread][ipart+2*npart];
            ( *probesArray )( fieldlocation[6], iPart_MPI )=Jloc_fields.x;
            ( *probesArray )( fieldlocation[7], iPart_MPI )=Jloc_fields.y;
            ( *probesArray )( fieldlocation[8], iPart_MPI )=Jloc_fields.z;
            ( *probesArray )( fieldlocation[9], iPart_MPI )=Rloc_fields;
            iPart_MPI++;
        }
        
        // Calculate Poynting flux on each point if needed
        // PoyX
        if( fieldlocation[14] < nFields ) {
            for( unsigned int iPart_MPI = offset_in_MPI[ipatch]; iPart_MPI<maxPart_MPI; iPart_MPI++ ) {
                ( *probesArray )( fieldlocation[14], iPart_MPI )
                    = ( *probesArray )( fieldlocation[1], iPart_MPI ) * ( *probesArray )( fieldlocation[5], iPart_MPI )
                    - ( *probesArray )( fieldlocation[2], iPart_MPI ) * ( *probesArray )( fieldlocation[4], iPart_MPI );
            }
        }
        // PoyY
        if( fieldlocation[15] < nFields ) {
            for( unsigned int iPart_MPI = offset_in_MPI[ipatch]; iPart_MPI<maxPart_MPI; iPart_MPI++ ) {
                ( *probesArray )( fieldlocation[15], iPart_MPI )
                    = ( *probesArray )( fieldlocation[2], iPart_MPI ) * ( *probesArray )( fieldlocation[3], iPart_MPI )
                    - ( *probesArray )( fieldlocation[0], iPart_MPI ) * ( *probesArray )( fieldlocation[5], iPart_MPI );
            }
        }
        // PoyZ
        if( fieldlocation[16] < nFields ) {
            for( unsigned int iPart_MPI = offset_in_MPI[ipatch]; iPart_MPI<maxPart_MPI; iPart_MPI++ ) {
                ( *probesArray )( fieldlocation[16], iPart_MPI )
                    = ( *probesArray )( fieldlocation[0], iPart_MPI ) * ( *probesArray )( fieldlocation[4], iPart_MPI )
                    - ( *probesArray )( fieldlocation[1], iPart_MPI ) * ( *probesArray )( fieldlocation[3], iPart_MPI );
            }
        }
        
        // Interpolate the species-related fields
        for( unsigned int ispec=0; ispec<species_field_index.size(); ispec++ ) {
            unsigned int start = patch->EMfields->species_starts[ispec];
            // In cylindrical geometry, all fields + all modes are interpolated
            // The unecessary results are discarded
            if( geometry == "AMcylindrical" ) {
                if( species_field_index[ispec].size() > 0 ) {
                    int istart( 0 ), iend( npart );
                    vector<double> dummy( npart );
                    vector<double *> loc( 4, &dummy[0] );
                    for( unsigned int j=0; j<species_field_index[ispec].size(); j++ ) {
                        unsigned int ifield = species_field_index[ispec][j];
                        unsigned int iloc = species_field_location[ispec][j];
                        loc[ifield] = &( ( *probesArray )( iloc, offset_in_MPI[ipatch] ) );
                    }
                    patch->probesInterp->oneField(
                        &patch->EMfields->allFields[start],
                        patch->probes[probe_n]->particles,
                        &istart, &iend,
                        loc[0], loc[1], loc[2], loc[3]
                    );
                }
            } else {
                for( unsigned int j=0; j<species_field_index[ispec].size(); j++ ) {
                    unsigned int ifield = species_field_index[ispec][j];
                    unsigned int iloc = species_field_location[ispec][j];
                    int istart( 0 ), iend( npart );
                    double *FieldLoc = &( ( *probesArray )( iloc, offset_in_MPI[ipatch] ) );
                    patch->probesInterp->oneField(
                        &patch->EMfields->allFields[start+ifield],
                        patch->probes[probe_n]->particles,
                        &istart, &iend,
                        FieldLoc
                    );
                }
            }
        }
        
        // Probes for envelope
        if( patch->EMfields->envelope != NULL ) {
            iPart_MPI = offset_in_MPI[ipatch];
            double Env_AabsLoc_fields, Env_ChiLoc_fields, Env_EabsLoc_fields, Env_ExabsLoc_fields;
            for( unsigned int ipart=0; ipart<npart; ipart++ ) {
                int iparticle( ipart ); // Compatibility
                patch->probesInterp->envelopeAndSusceptibility(
                    patch->EMfields,
                    patch->probes[probe_n]->particles,
                    iparticle,
                    &Env_AabsLoc_fields, &Env_ChiLoc_fields, &Env_EabsLoc_fields, &Env_ExabsLoc_fields
                );
                //! here we fill the probe data!!!
                ( *probesArray )( fieldlocation[10], iPart_MPI )=Env_AabsLoc_fields;
                ( *probesArray )( fieldlocation[11], iPart_MPI )=Env_ChiLoc_fields;
                ( *probesArray )( fieldlocation[12], iPart_MPI )=Env_EabsLoc_fields;
                ( *probesArray )( fieldlocation[13], iPart_MPI )=Env_ExabsLoc_fields;
                iPart_MPI++;
            } // END for ipart
        } // END if envelope
        
        // Accumulate the integral
        if( time_integral ) {
            patch->probes[probe_n]->integrated_data.resize( nFields );
            double dt = patch->EMfields->timestep;
            for( unsigned int i = 0; i < nFields; i++ ) {
                patch->probes[probe_n]->integrated_data[i].resize( npart, 0. );
                unsigned int iPart_MPI = offset_in_MPI[ipatch];
                for( unsigned int ipart=0; ipart<npart; ipart++ ) {
                    patch->probes[probe_n]->integrated_data[i][ipart] += dt * ( *probesArray )( i, iPart_MPI );
                    ( *probesArray )( i, iPart_MPI ) =  patch->probes[probe_n]->integrated_data[i][ipart];
                    iPart_MPI++;
                }
            }
        }
        
    } // END for ipatch
    
    #pragma omp master
    {
        if( timeSelection->theTimeIsNow( itime ) ) {
            // Define spaces
            H5Space memspace( {(hsize_t)nFields, nPart_MPI}, {}, {} );
            H5Space filespace( {(hsize_t)nFields, nPart_total_actual}, {0, offset_in_file[0]}, {(hsize_t)nFields, nPart_MPI} );
            // Create new dataset for this timestep
            H5Write d = file_->array( dataset_name, *(probesArray->data_), &filespace, &memspace, true );
            // Write x_moved
            d.attr( "x_moved", x_moved );
            
            delete probesArray;
            if( flush_timeSelection->theTimeIsNow( itime ) ) {
                file_->flush();
            }
        }
    }
    #pragma omp barrier
}

bool DiagnosticProbes::needsRhoJs( int itime )
{
    return hasRhoJs && timeSelection->theTimeIsNow( itime );
}

// SUPPOSED TO BE EXECUTED ONLY BY MASTER MPI
uint64_t DiagnosticProbes::getDiskFootPrint( int istart, int istop, Patch * )
{
    uint64_t footprint = 0;

    // Calculate the number of dumps between istart and istop
    uint64_t ndumps = timeSelection->howManyTimesBefore( istop ) - timeSelection->howManyTimesBefore( istart );

    // Add necessary global headers approximately
    footprint += 2200;
    if( dimProbe>0 ) {
        footprint += 2400;
    }
    if( ndumps>0 ) {
        footprint += ( uint64_t )( nDim_particle * nPart_total * 8 );
    }

    // Add local headers
    footprint += ndumps * ( uint64_t )( 480 + nFields * 6 );

    // Add size of each field
    footprint += ndumps * ( uint64_t )( nFields * nPart_total ) * 8;

    return footprint;
}
