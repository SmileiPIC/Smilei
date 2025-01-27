#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "PyTools.h"
#include "Species.h"
#include "ParticleData.h"
#include "Patch.h"
#include "SimWindow.h"
#include <algorithm>

// Class for each axis of the particle diags
class HistogramAxis
{
public:
    HistogramAxis() {};
    virtual ~HistogramAxis() {};
    
    void init( std::string, double, double, int, bool, bool, std::vector<double> );
    
    //! Function that goes through the particles and find where they should go in the axis
    virtual void calculate_locations( Species *, double *, int *, unsigned int, SimWindow * ) {};
    
    //! Print some info about the axis
    std::string info( std::string title = "" ) {
        std::ostringstream mystream( "" );
        if( title.empty() ) {
            title = type;
        }
        mystream << "Axis " << title;
        if( std::isnan(min) ) {
            mystream << " from auto";
        } else {
            mystream << " from " << min;
        }
        if( std::isnan(max) ) {
            mystream << " to auto";
        } else {
            mystream << " to " << max;
        }
        mystream << " in " << nbins << " steps";
        if( logscale ) {
            mystream << " [LOGSCALE] ";
        }
        if( edge_inclusive ) {
            mystream << " [EDGE INCLUSIVE]";
        }
        return mystream.str();
    }
    
    //! quantity of the axis (e.g. 'x', 'px', ...)
    std::string type;
    
    //! starting/ending point for the axis binning
    double min, max;
    //! starting/ending point for the axis binning, for all MPIs
    double global_min, global_max;
    //! number of bins for the axis binning
    int nbins;
    
    //! determines whether linear scale or log scale
    bool logscale;
    
    //! determines whether particles beyond min and max are counted in the first and last bin
    bool edge_inclusive;
    
    //! List of coefficients for some axes types
    std::vector<double> coefficients;
};


// Class for making a histogram of particle data
class Histogram
{
public:
    Histogram() {};
    virtual ~Histogram() {
        for( unsigned int iaxe=0; iaxe<axes.size(); iaxe++ ) {
            delete axes[iaxe];
        }
    };
    
    //! Compute the index of each particle in the final histogram
    void digitize( std::vector<Species *>, std::vector<double> &, std::vector<int> &, SimWindow * );
    //! Calculate the quantity of each particle to be summed in the histogram
    virtual void valuate( Species *, double *, int * ) {
        ERROR( "`deposited_quantity` should not be empty" );
    };
    //! Same as `valuate` for all species
    void valuate( std::vector<Species *> species, std::vector<double> &double_buffer, std::vector<int> &int_buffer ) {
        unsigned int istart = 0;
        for( unsigned int ispec=0; ispec < species.size(); ispec++ ) {
            valuate( species[ispec], &double_buffer[istart], &int_buffer[istart] );
            istart += species[ispec]->getNbrOfParticles();
        }
    };
    //! Add the contribution of each particle in the histogram
    void distribute( std::vector<double> &, std::vector<int> &, std::vector<double> & );

    std::string deposited_quantity;

    std::vector<HistogramAxis *> axes;
};



class HistogramAxis_x : public HistogramAxis
{
    ~HistogramAxis_x() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Position[0][ipart];
        }
    };
};
class HistogramAxis_moving_x : public HistogramAxis
{
    ~HistogramAxis_moving_x() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow *simWindow )
    {
        double x_moved = simWindow->getXmoved();
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Position[0][ipart]-x_moved;
        }
    };
};
class HistogramAxis_y : public HistogramAxis
{
    ~HistogramAxis_y() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Position[1][ipart];
        }
    };
};
class HistogramAxis_z : public HistogramAxis
{
    ~HistogramAxis_z() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Position[2][ipart];
        }
    };
};
class HistogramAxis_vector : public HistogramAxis
{
    ~HistogramAxis_vector() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        unsigned int idim, ndim = coefficients.size()/2;
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = 0.;
            for( idim=0; idim<ndim; idim++ ) {
                array[ipart] += ( s->particles->Position[idim][ipart] - coefficients[idim] ) * coefficients[idim+ndim];
            }
        }
    };
};
class HistogramAxis_theta2D : public HistogramAxis
{
    ~HistogramAxis_theta2D() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        double X, Y;
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            X = s->particles->Position[0][ipart] - coefficients[0];
            Y = s->particles->Position[1][ipart] - coefficients[1];
            array[ipart] = atan2( coefficients[2]*Y - coefficients[3]*X, coefficients[2]*X + coefficients[3]*Y );
        }
    };
};
class HistogramAxis_theta3D : public HistogramAxis
{
    ~HistogramAxis_theta3D() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = ( s->particles->Position[0][ipart] - coefficients[0] ) * coefficients[3]
                           + ( s->particles->Position[1][ipart] - coefficients[1] ) * coefficients[4]
                           + ( s->particles->Position[2][ipart] - coefficients[2] ) * coefficients[5];
            if( array[ipart]> 1. ) {
                array[ipart] = 0.;
            } else if( array[ipart]<-1. ) {
                array[ipart] = M_PI;
            } else {
                array[ipart] = acos( array[ipart] );
            }
        }
    };
};
class HistogramAxis_phi : public HistogramAxis
{
    ~HistogramAxis_phi() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        unsigned int idim;
        double a, b;
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            a = 0.;
            b = 0.;
            for( idim=0; idim<3; idim++ ) {
                a += ( s->particles->Position[idim][ipart] - coefficients[idim] ) * coefficients[idim+3];
                b += ( s->particles->Position[idim][ipart] - coefficients[idim] ) * coefficients[idim+6];
            }
            array[ipart] = atan2( b, a );
        }
    };
};
class HistogramAxis_px : public HistogramAxis
{
    ~HistogramAxis_px() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Momentum[0][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[0][ipart];
            }
        }
    };
};
class HistogramAxis_py : public HistogramAxis
{
    ~HistogramAxis_py() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Momentum[1][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[1][ipart];
            }
        }
    };
};
class HistogramAxis_pz : public HistogramAxis
{
    ~HistogramAxis_pz() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Momentum[2][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[2][ipart];
            }
        }
    };
};
class HistogramAxis_p : public HistogramAxis
{
    ~HistogramAxis_p() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                               + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                               + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_gamma : public HistogramAxis
{
    ~HistogramAxis_gamma() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_ekin : public HistogramAxis
{
    ~HistogramAxis_ekin() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * ( sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] ) - 1. );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_vx : public HistogramAxis
{
    ~HistogramAxis_vx() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[0][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[0][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_vy : public HistogramAxis
{
    ~HistogramAxis_vy() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[1][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[1][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_vz : public HistogramAxis
{
    ~HistogramAxis_vz() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[2][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Momentum[2][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_v : public HistogramAxis
{
    ~HistogramAxis_v() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = 1.0 / sqrt( 1. + 1./( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] ) );
        }
    };
};
class HistogramAxis_vperp2 : public HistogramAxis
{
    ~HistogramAxis_vperp2() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = ( s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                 + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart]
                               ) / ( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = ( s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                 + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart]
                               ) / ( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class HistogramAxis_charge : public HistogramAxis
{
    ~HistogramAxis_charge() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = ( double ) s->particles->Charge[ipart];
        }
    };
};
class HistogramAxis_chi : public HistogramAxis
{
    ~HistogramAxis_chi() {};
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Chi[ipart];
        }
    };
};
#ifdef SMILEI_USE_NUMPY
class HistogramAxis_user_function : public HistogramAxis
{
public:
    HistogramAxis_user_function( PyObject *type_object ) :
        HistogramAxis(),
        function( type_object )
    {};
    ~HistogramAxis_user_function()
    {
        SMILEI_PY_ACQUIRE_GIL
        Py_DECREF( function );
        SMILEI_PY_RELEASE_GIL
    };
private:
    void calculate_locations( Species *s, double *array, int *index, unsigned int npart, SimWindow * )
    {
        PyArrayObject *ret;
        SMILEI_PY_ACQUIRE_GIL
        {
            // Expose particle data as numpy arrays
            ParticleData particleData( npart );
            particleData.set( s->particles );
            // run the function
            ret = ( PyArrayObject * )PyObject_CallFunctionObjArgs( function, particleData.get(), NULL );
        }
        // Copy the result to "array"
        double *arr = ( double * ) PyArray_GETPTR1( ret, 0 );
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ )
        {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = arr[ipart];
        }
        Py_DECREF( ret );
        SMILEI_PY_RELEASE_GIL
    };

    PyObject *function;
};
#endif

//! Children classes, for various manners to fill the histogram
class Histogram_number : public Histogram
{
    ~Histogram_number() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Weight[ipart];
        }
    };
};
class Histogram_charge : public Histogram
{
    ~Histogram_charge() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Weight[ipart] * ( double )( s->particles->Charge[ipart] );
        }
    };
};
class Histogram_jx : public Histogram
{
    ~Histogram_jx() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * ( double )( s->particles->Charge[ipart] )
                               * s->particles->Momentum[0][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_jy : public Histogram
{
    ~Histogram_jy() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * ( double )( s->particles->Charge[ipart] )
                               * s->particles->Momentum[1][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[1][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_jz : public Histogram
{
    ~Histogram_jz() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * ( double )( s->particles->Charge[ipart] )
                               * s->particles->Momentum[2][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[2][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_ekin : public Histogram
{
    ~Histogram_ekin() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * ( sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] ) - 1. );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * ( sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] ) );
            }
        }
    };
};
//! Children class of Histogram: for the quantum parameter
//! of the radiating particles
class Histogram_chi : public Histogram
{
    ~Histogram_chi() {};
public:
    Histogram_chi( Patch *patch, std::vector<unsigned int> &species, std::string errorPrefix )
        : Histogram()
    {
        // The requested species must be radiating
        for( unsigned int ispec=0 ; ispec < species.size() ; ispec++ )
            if( ! patch->vecSpecies[species[ispec]]->particles->has_quantum_parameter ) {
                ERROR( errorPrefix << ": 'weight_chi' requires all species to be radiating" );
            }
    };
private:
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = s->particles->Weight[ipart]
                           * s->particles->Chi[ipart];
        }
    };
};
class Histogram_p : public Histogram
{
    ~Histogram_p() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_px : public Histogram
{
    ~Histogram_px() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart] * s->particles->Momentum[0][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * s->particles->Momentum[0][ipart];
            }
        }
    };
};
class Histogram_py : public Histogram
{
    ~Histogram_py() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart] * s->particles->Momentum[1][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * s->particles->Momentum[1][ipart];
            }
        }
    };
};
class Histogram_pz : public Histogram
{
    ~Histogram_pz() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart] * s->particles->Momentum[2][ipart];
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart] * s->particles->Momentum[2][ipart];
            }
        }
    };
};
class Histogram_pressure_xx : public Histogram
{
    ~Histogram_pressure_xx() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * ( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart] )
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * ( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart] )
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_pressure_yy : public Histogram
{
    ~Histogram_pressure_yy() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * ( s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart] )
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * ( s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart] )
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_pressure_zz : public Histogram
{
    ~Histogram_pressure_zz() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * ( s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] )
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * ( s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] )
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_pressure_xy : public Histogram
{
    ~Histogram_pressure_xy() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               * s->particles->Momentum[1][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               * s->particles->Momentum[1][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_pressure_xz : public Histogram
{
    ~Histogram_pressure_xz() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               * s->particles->Momentum[2][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               * s->particles->Momentum[2][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_pressure_yz : public Histogram
{
    ~Histogram_pressure_yz() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * s->particles->Momentum[1][ipart]
                               * s->particles->Momentum[2][ipart]
                               / sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[1][ipart]
                               * s->particles->Momentum[2][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};
class Histogram_ekin_vx : public Histogram
{
    ~Histogram_ekin_vx() {};
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        // Matter Particles
        if( s->mass_ > 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->mass_ * s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               * ( 1. - 1./sqrt( 1. + s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] ) );
            }
        }
        // Photons
        else if( s->mass_ == 0 ) {
            for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
                if( index[ipart]<0 ) {
                    continue;
                }
                array[ipart] = s->particles->Weight[ipart]
                               * s->particles->Momentum[0][ipart]
                               / sqrt( s->particles->Momentum[0][ipart] * s->particles->Momentum[0][ipart]
                                     + s->particles->Momentum[1][ipart] * s->particles->Momentum[1][ipart]
                                     + s->particles->Momentum[2][ipart] * s->particles->Momentum[2][ipart] );
            }
        }
    };
};

#ifdef SMILEI_USE_NUMPY
class Histogram_user_function : public Histogram
{
public:
    Histogram_user_function( PyObject *deposited_quantity_object ) :
        Histogram(),
        function( deposited_quantity_object )
    {};
    ~Histogram_user_function()
    {
        SMILEI_PY_ACQUIRE_GIL
        Py_DECREF( function );
        SMILEI_PY_RELEASE_GIL
    };
private:
    void valuate( Species *s, double *array, int *index )
    {
        unsigned int npart = s->getNbrOfParticles();
        PyArrayObject *ret;
        SMILEI_PY_ACQUIRE_GIL
        // Expose particle data as numpy arrays
        {
            ParticleData particleData( npart );
            particleData.set( s->particles );
            // run the function
            ret = ( PyArrayObject * )PyObject_CallFunctionObjArgs( function, particleData.get(), NULL );
        }
        // Copy the result to "array"
        double *arr = ( double * ) PyArray_GETPTR1( ret, 0 );
        for( unsigned int ipart = 0 ; ipart < npart ; ipart++ ) {
            if( index[ipart]<0 ) {
                continue;
            }
            array[ipart] = arr[ipart];
        }
        Py_DECREF( ret );
        SMILEI_PY_RELEASE_GIL
    };
    
    PyObject *function;
};
#endif

#endif
