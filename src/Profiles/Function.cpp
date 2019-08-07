#include "Function.h"
#include <complex>
#include <cmath>

using namespace std;

// Functions to evaluate a python function with various numbers of arguments
// 1D
double Function_Python1D::valueAt( double time )
{
    return PyTools::runPyFunction( py_profile, time );
}
double Function_Python1D::valueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction( py_profile, time );
}
double Function_Python1D::valueAt( vector<double> x_cell )
{
    return PyTools::runPyFunction( py_profile, x_cell[0] );
}

// 2D
double Function_Python2D::valueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction( py_profile, x_cell[0], time );
}
double Function_Python2D::valueAt( vector<double> x_cell )
{
    return PyTools::runPyFunction( py_profile, x_cell[0], x_cell[1] );
}
// 2D complex
std::complex<double> Function_Python2D::complexValueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction_complex( py_profile, x_cell[0], time );
}
std::complex<double> Function_Python2D::complexValueAt( vector<double> x_cell )
{
    return PyTools::runPyFunction_complex( py_profile, x_cell[0], x_cell[1] );
}

// 3D
double Function_Python3D::valueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction( py_profile, x_cell[0], x_cell[1], time );
}
double Function_Python3D::valueAt( vector<double> x_cell )
{
    return PyTools::runPyFunction( py_profile, x_cell[0], x_cell[1], x_cell[2] );
}
// 3D complex
std::complex<double> Function_Python3D::complexValueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction_complex( py_profile, x_cell[0], x_cell[1], time );
}

// 4D
double Function_Python4D::valueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction( py_profile, x_cell[0], x_cell[1], x_cell[2], time );
}
// 4D complex
std::complex<double> Function_Python4D::complexValueAt( vector<double> x_cell, double time )
{
    return PyTools::runPyFunction_complex( py_profile, x_cell[0], x_cell[1], x_cell[2], time );
}

// Special cases for locations specified in numpy arrays
#ifdef SMILEI_USE_NUMPY
PyArrayObject *Function_Python1D::valueAt( std::vector<PyArrayObject *> x )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], NULL );
}
PyArrayObject *Function_Python2D::valueAt( std::vector<PyArrayObject *> x )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], NULL );
}
PyArrayObject *Function_Python3D::valueAt( std::vector<PyArrayObject *> x )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], x[2], NULL );
}
PyArrayObject *Function_Python4D::complexValueAt( std::vector<PyArrayObject *> x, PyArrayObject *t )
{
    PyObject *values = PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], x[2], t, NULL );
    PyArrayObject *cvalues = ( PyArrayObject * )PyObject_CallMethod( values, "astype", "s", "complex", NULL );
    Py_DECREF( values );
    return cvalues;
}
PyArrayObject *Function_Python2D::complexValueAt( std::vector<PyArrayObject *> x )
{
    PyObject *values = PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], NULL );
    PyArrayObject *cvalues = ( PyArrayObject * )PyObject_CallMethod( values, "astype", "s", "complex", NULL );
    Py_DECREF( values );
    return cvalues;
}
// Time dependent

PyArrayObject *Function_Python1D::valueAt( std::vector<PyArrayObject *> x, double time )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], time, NULL );
}
PyArrayObject *Function_Python2D::valueAt( std::vector<PyArrayObject *> x, double time  )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], time, NULL );
}
PyArrayObject *Function_Python3D::valueAt( std::vector<PyArrayObject *> x, double time  )
{
    return ( PyArrayObject * )PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], x[2], time, NULL );
}
PyArrayObject *Function_Python4D::complexValueAt( std::vector<PyArrayObject *> x, double time )
{
    PyObject *values = PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], x[2], time, NULL );
    PyArrayObject *cvalues = ( PyArrayObject * )PyObject_CallMethod( values, "astype", "s", "complex", NULL );
    Py_DECREF( values );
    return cvalues;
}
// PyArrayObject *Function_Python2D::complexValueAt( std::vector<PyArrayObject *> x )
// {
//     PyObject *values = PyObject_CallFunctionObjArgs( py_profile, x[0], x[1], NULL );
//     PyArrayObject *cvalues = ( PyArrayObject * )PyObject_CallMethod( values, "astype", "s", "complex", NULL );
//     Py_DECREF( values );
//     return cvalues;
// }

#endif


// Constant profiles
double Function_Constant1D::valueAt( vector<double> x_cell )
{
    return ( x_cell[0]>xvacuum ) ? value : 0.;
}
double Function_Constant2D::valueAt( vector<double> x_cell )
{
    return ( ( x_cell[0]>xvacuum ) && ( x_cell[1]>yvacuum ) ) ? value : 0.;
}
double Function_Constant3D::valueAt( vector<double> x_cell )
{
    return ( ( x_cell[0]>xvacuum ) && ( x_cell[1]>yvacuum ) && ( x_cell[2]>zvacuum ) ) ? value : 0.;
}

// Trapezoidal profiles
inline double trapeze( double x, double plateau, double slope1, double slope2, double invslope1, double invslope2 )
{
    double result = 0.;
    if( x > 0. ) {
        if( x < slope1 ) {
            result = invslope1 * x;
        } else {
            x -= slope1 + plateau;
            if( x < 0. ) {
                result = 1.;
            } else {
                x -= slope2;
                if( x < 0. ) {
                    result = - x * invslope2;
                }
            }
        }
    }
    return result;
}
double Function_Trapezoidal1D::valueAt( vector<double> x_cell )
{
    return value * trapeze( x_cell[0]-xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2 );
}
double Function_Trapezoidal2D::valueAt( vector<double> x_cell )
{
    return value
           * trapeze( x_cell[0]-xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2 )
           * trapeze( x_cell[1]-yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2 );
}
double Function_Trapezoidal3D::valueAt( vector<double> x_cell )
{
    return value
           * trapeze( x_cell[0]-xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2 )
           * trapeze( x_cell[1]-yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2 )
           * trapeze( x_cell[2]-zvacuum, zplateau, zslope1, zslope2, invzslope1, invzslope2 );
}

// Gaussian profiles
double Function_Gaussian1D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0], xfactor=0.;
    if( x > xvacuum  && x < xvacuum+xlength ) {
        xfactor = exp( -pow( x-xcenter, xorder ) * invxsigma );
    }
    return value * xfactor;
}
double Function_Gaussian2D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0], xfactor=0.;
    double y = x_cell[1], yfactor=0.;
    if( x > xvacuum  && x < xvacuum+xlength ) {
        xfactor = exp( -pow( x-xcenter, xorder ) * invxsigma );
    }
    if( y > yvacuum  && y < yvacuum+ylength ) {
        yfactor = exp( -pow( y-ycenter, yorder ) * invysigma );
    }
    return value * xfactor * yfactor;
}
double Function_Gaussian3D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0], xfactor=0.;
    double y = x_cell[1], yfactor=0.;
    double z = x_cell[2], zfactor=0.;
    if( x > xvacuum  && x < xvacuum+xlength ) {
        xfactor = exp( -pow( x-xcenter, xorder ) * invxsigma );
    }
    if( y > yvacuum  && y < yvacuum+ylength ) {
        yfactor = exp( -pow( y-ycenter, yorder ) * invysigma );
    }
    if( z > zvacuum  && z < zvacuum+zlength ) {
        zfactor = exp( -pow( z-zcenter, zorder ) * invzsigma );
    }
    return value * xfactor * yfactor * zfactor;
}

// Polygonal profiles
double Function_Polygonal1D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0];
    if( x < xpoints[0] ) {
        return 0.;
    }
    for( int i=1; i<npoints; i++ )
        if( x < xpoints[i] ) {
            return xvalues[i-1] + xslopes[i-1] * ( x - xpoints[i-1] );
        }
    return 0.;
}
double Function_Polygonal2D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0];
    if( x < xpoints[0] ) {
        return 0.;
    }
    for( int i=1; i<npoints; i++ )
        if( x < xpoints[i] ) {
            return xvalues[i-1] + xslopes[i-1] * ( x - xpoints[i-1] );
        }
    return 0.;
}
double Function_Polygonal3D::valueAt( vector<double> x_cell )
{
    double x = x_cell[0];
    if( x < xpoints[0] ) {
        return 0.;
    }
    for( int i=1; i<npoints; i++ )
        if( x < xpoints[i] ) {
            return xvalues[i-1] + xslopes[i-1] * ( x - xpoints[i-1] );
        }
    return 0.;
}

// Cosine profiles
double Function_Cosine1D::valueAt( vector<double> x_cell )
{
    double x = ( x_cell[0] - xvacuum ) * invxlength, xfactor = 0.;
    if( x > 0. && x < 1. ) {
        xfactor = base + xamplitude * cos( xphi + xnumber2pi * x );
    }
    return xfactor;
}
double Function_Cosine2D::valueAt( vector<double> x_cell )
{
    double x = ( x_cell[0] - xvacuum ) * invxlength, xfactor = 0.;
    double y = ( x_cell[1] - yvacuum ) * invylength, yfactor = 0.;
    if( x > 0. && x < 1. ) {
        xfactor = base + xamplitude * cos( xphi + xnumber2pi * x );
    }
    if( y > 0. && y < 1. ) {
        yfactor = base + yamplitude * cos( yphi + ynumber2pi * y );
    }
    return xfactor * yfactor;
}
double Function_Cosine3D::valueAt( vector<double> x_cell )
{
    double x = ( x_cell[0] - xvacuum ) * invxlength, xfactor = 0.;
    double y = ( x_cell[1] - yvacuum ) * invylength, yfactor = 0.;
    double z = ( x_cell[2] - zvacuum ) * invzlength, zfactor = 0.;
    if( x > 0. && x < 1. ) {
        xfactor = base + xamplitude * cos( xphi + xnumber2pi * x );
    }
    if( y > 0. && y < 1. ) {
        yfactor = base + yamplitude * cos( yphi + ynumber2pi * y );
    }
    if( z > 0. && z < 1. ) {
        zfactor = base + zamplitude * cos( zphi + znumber2pi * z );
    }
    return xfactor * yfactor * zfactor;
}

// Polynomial profiles
double Function_Polynomial1D::valueAt( vector<double> x_cell )
{
    double r = 0., xx0 = x_cell[0]-x0, xx = 1.;
    unsigned int currentOrder = 0;
    for( unsigned int i=0; i<n_orders; i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            xx *= xx0;
        }
        r += coeffs[i][0] * xx;
    }
    return r;
}
double Function_Polynomial2D::valueAt( vector<double> x_cell )
{
    double r = 0., xx0 = x_cell[0]-x0, yy0 = x_cell[1]-y0;
    vector<double> xx;
    unsigned int currentOrder = 0, j;
    xx.resize( n_coeffs );
    xx[0] = 1.;
    for( unsigned int i=0; i<n_orders; i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            j = currentOrder;
            xx[j] = xx[j-1] * yy0;
            do {
                j--;
                xx[j] *= xx0;
            } while( j>0 );
        }
        for( j=0; j<=orders[i]; j++ ) {
            r += coeffs[i][j] * xx[j];
        }
    }
    return r;
}
double Function_Polynomial3D::valueAt( vector<double> x_cell )
{
    double r = 0., xx0 = x_cell[0]-x0, yy0 = x_cell[1]-y0, zz0 = x_cell[2]-z0;
    vector<double> xx;
    unsigned int currentOrder = 0, current_n_coeffs = 1, j, k;
    xx.resize( n_coeffs );
    xx[0] = 1.;
    for( unsigned int i=0; i<n_orders; i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            k = current_n_coeffs-1;
            j = current_n_coeffs+currentOrder;
            xx[j] = xx[k] * zz0;
            do {
                j--;
                xx[j] = xx[k] * yy0;
                k--;
            } while( j>current_n_coeffs );
            do {
                j--;
                xx[j] = xx[j] * xx0;
            } while( j>0 );
            current_n_coeffs += currentOrder+1;
        }
        for( j=0; j<current_n_coeffs; j++ ) {
            r += coeffs[i][j] * xx[j];
        }
    }
    return r;
}

// Time constant profile
double Function_TimeConstant::valueAt( double time )
{
    if( time > start ) {
        return 1.;
    } else {
        return 0.;
    }
}

// Time trapezoidal profile
double Function_TimeTrapezoidal::valueAt( double time )
{
    return trapeze( time-start, plateau, slope1, slope2, invslope1, invslope2 );
}

// Time gaussian profile
double Function_TimeGaussian::valueAt( double time )
{
    if( time < start || time > end ) {
        return 0.;
    } else {
        return exp( -pow( time-center, order ) * invsigma );
    }
}

// Time polygonal profile
double Function_TimePolygonal::valueAt( double time )
{
    if( time < points[0] ) {
        return 0.;
    }
    for( int i=1; i<npoints; i++ )
        if( time<points[i] ) {
            return values[i-1] + slopes[i-1] * ( time-points[i-1] );
        }
    return 0.;
}

// Time cosine profile
double Function_TimeCosine::valueAt( double time )
{
    if( time < start || time > end ) {
        return 0.;
    } else {
        return base + amplitude * cos( phi + freq*( time-start ) );
    }
}

// Time polynomial profile
double Function_TimePolynomial::valueAt( double time )
{
    double r = 0., tt0 = time-t0, tt = 1.;
    int currentOrder = 0;
    for( unsigned int i=0; i<orders.size(); i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            tt *= tt0;
        }
        r += coeffs[i] * tt;
    }
    return r;
}
// Time sin2 profile with a plateau the funtion and the derivative are continuos
double Function_TimeSin2Plateau::valueAt( double time )
{
    if( time<start ) {
        return 0.;
    } else if( ( time < start+slope1 )&&( slope1!=0. ) ) {
        return pow( sin( 0.5*M_PI*( time-start )/slope1 ), 2 );
    } else if( time < start+slope1+plateau ) {
        return 1.;
    } else if( ( time < start+slope1+plateau+slope2 )&&( slope2!=0. ) ) {
        return pow( cos( 0.5*M_PI*( time-start-slope1-plateau )/slope2 ), 2 );
    } else {
        return 0.;
    }
}
