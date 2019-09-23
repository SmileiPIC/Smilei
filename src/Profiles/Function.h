#ifndef Function_H
#define Function_H

#include "PyTools.h"
#include <vector>
#include <string>
#include <complex>

class Function
{
public:
    //! Default constructor
    Function() {};
    //! Default destructor
    virtual ~Function() {};
    
    //! Gets the value of a N-D function at a point located by its coordinates in a vector
    virtual double valueAt( std::vector<double> )
    {
        return 0.; // virtual => will be redefined
    };
    
    //! Gets the value of a 1-D function at a point located by a double
    virtual double valueAt( double x )
    {
        // just in case someone uses 1D space profile instead of time
        std::vector<double> v( 1 );
        v[0] = x;
        return valueAt( v );
    };
    
    //! Gets the value of a N-D function from both a vector and a double. The double is the last argument.
    virtual double valueAt( std::vector<double>, double )
    {
        return 0.; // virtual => will be redefined
    };
    
    //! Gets the complex value of a N-D function from both a vector and a double. The double is the last argument.
    virtual std::complex<double> complexValueAt( std::vector<double>, double )
    {
        return 0.; // virtual => will be redefined
    };
    
    //! Gets the complex value of a N-D function from a vector.
    virtual std::complex<double> complexValueAt( std::vector<double> )
    {
        return 0.; // virtual => will be redefined
    };
    
    //! Provide information about the function
    virtual std::string getInfo()
    {
        std::string info = "";
        return info; // virtual => will be redefined
    };
    
#ifdef SMILEI_USE_NUMPY
    //! Gets the value of an N-D function at points specified as numpy arrays
    virtual PyArrayObject *valueAt( std::vector<PyArrayObject *> )
    {
        return NULL;
    };
    //! Gets the value of an N-D function at points specified as numpy arrays
    virtual PyArrayObject *complexValueAt( std::vector<PyArrayObject *> )
    {
        return NULL;
    };
    //! Gets the value of an N-D function at points specified as numpy arrays
    virtual PyArrayObject *complexValueAt( std::vector<PyArrayObject *>, PyArrayObject * )
    {
        return NULL;
    };
#endif
};


// Children classes for python functions

class Function_Python1D : public Function
{
public:
    Function_Python1D( PyObject *pp ) : py_profile( pp ) {};
    Function_Python1D( Function_Python1D *f ) : py_profile( f->py_profile ) {};
    double valueAt( double ); // time
    double valueAt( std::vector<double>, double ); // time (space discarded)
    double valueAt( std::vector<double> ); // space
#ifdef SMILEI_USE_NUMPY
    PyArrayObject *valueAt( std::vector<PyArrayObject *> ); // numpy
#endif
private:
    PyObject *py_profile;
};


class Function_Python2D : public Function
{
public:
    Function_Python2D( PyObject *pp ) : py_profile( pp ) {};
    Function_Python2D( Function_Python2D *f ) : py_profile( f->py_profile ) {};
    double valueAt( std::vector<double>, double ); // space + time
    double valueAt( std::vector<double> ); // space
    std::complex<double> complexValueAt( std::vector<double>, double ); // space + time
    std::complex<double> complexValueAt( std::vector<double> ); // space
#ifdef SMILEI_USE_NUMPY
    PyArrayObject *valueAt( std::vector<PyArrayObject *> ); // numpy
    PyArrayObject *complexValueAt( std::vector<PyArrayObject *> ); // numpy
#endif
private:
    PyObject *py_profile;
};

class Function_Python2D_Complex : public Function
{
public:
    Function_Python2D_Complex( PyObject *pp ) : py_profile( pp ) {};
    Function_Python2D_Complex( Function_Python2D_Complex *f ) : py_profile( f->py_profile ) {};
    std::complex<double> valueAtComplex( std::vector<double>, double ); // space + time
private:
    PyObject *py_profile;
};

class Function_Python3D : public Function
{
public:
    Function_Python3D( PyObject *pp ) : py_profile( pp ) {};
    Function_Python3D( Function_Python3D *f ) : py_profile( f->py_profile ) {};
    double valueAt( std::vector<double>, double ); // space + time
    double valueAt( std::vector<double> ); // space
    std::complex<double> complexValueAt( std::vector<double>, double ); // space + time
#ifdef SMILEI_USE_NUMPY
    PyArrayObject *valueAt( std::vector<PyArrayObject *> ); // numpy
#endif
private:
    PyObject *py_profile;
};

class Function_Python3D_Complex : public Function
{
public:
    Function_Python3D_Complex( PyObject *pp ) : py_profile( pp ) {};
    Function_Python3D_Complex( Function_Python3D_Complex *f ) : py_profile( f->py_profile ) {};
    std::complex<double> valueAtComplex( std::vector<double>, double ); // space + time
private:
    PyObject *py_profile;
};

class Function_Python4D : public Function
{
public:
    Function_Python4D( PyObject *pp ) : py_profile( pp ) {};
    Function_Python4D( Function_Python4D *f ) : py_profile( f->py_profile ) {};
    double valueAt( std::vector<double>, double ); // space + time
    std::complex<double> complexValueAt( std::vector<double>, double ); // space + time
#ifdef SMILEI_USE_NUMPY
    PyArrayObject *complexValueAt( std::vector<PyArrayObject *>, PyArrayObject * ); // numpy
#endif
private:
    PyObject *py_profile;
};
/*
class Function_Python4D_Complex : public Function
{
public:
    Function_Python4D_Complex(PyObject *pp) : py_profile(pp) {};
    Function_Python4D_Complex(Function_Python4D_Complex *f) : py_profile(f->py_profile) {};
    std::complex<double> valueAtComplex(std::vector<double>, double); // space + time
//#ifdef SMILEI_USE_NUMPY
    //PyArrayObject* complexValueAt(std::vector<PyArrayObject*>,std::vector<PyArrayObject*>); // numpy
//#endif
private:
    PyObject *py_profile;
};
*/


// Children classes for hard-coded functions

class Function_Constant1D : public Function
{
public:
    Function_Constant1D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
    };
    Function_Constant1D( Function_Constant1D *f )
    {
        value   = f->value  ;
        xvacuum = f->xvacuum;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value) + ")";
        return info;
    };
private:
    double value, xvacuum;
};


class Function_Constant2D : public Function
{
public:
    Function_Constant2D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
    };
    Function_Constant2D( Function_Constant2D *f )
    {
        value   = f->value  ;
        xvacuum = f->xvacuum;
        yvacuum = f->yvacuum;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value) + ")";
        return info;
    };
private:
    double value, xvacuum, yvacuum;
};


class Function_Constant3D : public Function
{
public:
    Function_Constant3D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "zvacuum", zvacuum );
    };
    Function_Constant3D( Function_Constant3D *f )
    {
        value   = f->value  ;
        xvacuum = f->xvacuum;
        yvacuum = f->yvacuum;
        zvacuum = f->zvacuum;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value) + ")";
        return info;
    };
private:
    double value, xvacuum, yvacuum, zvacuum;
};


class Function_Trapezoidal1D : public Function
{
public:
    Function_Trapezoidal1D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xplateau", xplateau );
        PyTools::getAttr( py_profile, "xslope1", xslope1 );
        PyTools::getAttr( py_profile, "xslope2", xslope2 );
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
    };
    Function_Trapezoidal1D( Function_Trapezoidal1D *f )
    {
        value    = f->value   ;
        xvacuum  = f->xvacuum ;
        xplateau = f->xplateau;
        xslope1  = f->xslope1 ;
        xslope2  = f->xslope2 ;
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", xplateau: " + std::to_string(xplateau)
                         + ", xslope1: " + std::to_string(xslope1)
                         + ", xslope2: " + std::to_string(xslope2)
                         + ")";
        return info;
    };
private:
    double value, xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2;
};


class Function_Trapezoidal2D : public Function
{
public:
    Function_Trapezoidal2D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xplateau", xplateau );
        PyTools::getAttr( py_profile, "xslope1", xslope1 );
        PyTools::getAttr( py_profile, "xslope2", xslope2 );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "yplateau", yplateau );
        PyTools::getAttr( py_profile, "yslope1", yslope1 );
        PyTools::getAttr( py_profile, "yslope2", yslope2 );
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
        invyslope1 = 1./yslope1;
        invyslope2 = 1./yslope2;
    };
    Function_Trapezoidal2D( Function_Trapezoidal2D *f )
    {
        value    = f->value   ;
        xvacuum  = f->xvacuum ;
        xplateau = f->xplateau;
        xslope1  = f->xslope1 ;
        xslope2  = f->xslope2 ;
        yvacuum  = f->yvacuum ;
        yplateau = f->yplateau;
        yslope1  = f->yslope1 ;
        yslope2  = f->yslope2 ;
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
        invyslope1 = 1./yslope1;
        invyslope2 = 1./yslope2;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", yvacuum: " + std::to_string(yvacuum)
                         + ", xplateau: " + std::to_string(xplateau)
                         + ", yplateau: " + std::to_string(yplateau)
                         + ", xslope1: " + std::to_string(xslope1)
                         + ", xslope2: " + std::to_string(xslope2)
                         + ", yslope1: " + std::to_string(yslope1)
                         + ", yslope2: " + std::to_string(yslope2)
                         + ")";
        return info;
    };
private:
    double value,
           xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2,
           yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2;
};


class Function_Trapezoidal3D : public Function
{
public:
    Function_Trapezoidal3D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xplateau", xplateau );
        PyTools::getAttr( py_profile, "xslope1", xslope1 );
        PyTools::getAttr( py_profile, "xslope2", xslope2 );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "yplateau", yplateau );
        PyTools::getAttr( py_profile, "yslope1", yslope1 );
        PyTools::getAttr( py_profile, "yslope2", yslope2 );
        PyTools::getAttr( py_profile, "zvacuum", zvacuum );
        PyTools::getAttr( py_profile, "zplateau", zplateau );
        PyTools::getAttr( py_profile, "zslope1", zslope1 );
        PyTools::getAttr( py_profile, "zslope2", zslope2 );
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
        invyslope1 = 1./yslope1;
        invyslope2 = 1./yslope2;
        invzslope1 = 1./zslope1;
        invzslope2 = 1./zslope2;
    };
    Function_Trapezoidal3D( Function_Trapezoidal3D *f )
    {
        value    = f->value   ;
        xvacuum  = f->xvacuum ;
        xplateau = f->xplateau;
        xslope1  = f->xslope1 ;
        xslope2  = f->xslope2 ;
        yvacuum  = f->yvacuum ;
        yplateau = f->yplateau;
        yslope1  = f->yslope1 ;
        yslope2  = f->yslope2 ;
        zvacuum  = f->zvacuum ;
        zplateau = f->zplateau;
        zslope1  = f->zslope1 ;
        zslope2  = f->zslope2 ;
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
        invyslope1 = 1./yslope1;
        invyslope2 = 1./yslope2;
        invzslope1 = 1./zslope1;
        invzslope2 = 1./zslope2;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", yvacuum: " + std::to_string(yvacuum)
                         + ", zvacuum: " + std::to_string(zvacuum)
                         + ", xplateau: " + std::to_string(xplateau)
                         + ", yplateau: " + std::to_string(yplateau)
                         + ", zplateau: " + std::to_string(zplateau)
                         + ", xslope1: " + std::to_string(xslope1)
                         + ", xslope2: " + std::to_string(xslope2)
                         + ", yslope1: " + std::to_string(yslope1)
                         + ", yslope2: " + std::to_string(yslope2)
                         + ", zslope1: " + std::to_string(zslope1)
                         + ", zslope2: " + std::to_string(zslope2)
                         + ")";
        return info;
    };
private:
    double value,
           xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2,
           yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2,
           zvacuum, zplateau, zslope1, zslope2, invzslope1, invzslope2;
};


class Function_Gaussian1D : public Function
{
public:
    Function_Gaussian1D( PyObject *py_profile )
    {
        double xsigma( 0 );
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xsigma", xsigma );
        PyTools::getAttr( py_profile, "xcenter", xcenter );
        PyTools::getAttr( py_profile, "xorder", xorder );
        invxsigma = 1./xsigma;
    };
    Function_Gaussian1D( Function_Gaussian1D *f )
    {
        value     = f->value  ;
        xvacuum   = f->xvacuum;
        xlength   = f->xlength;
        invxsigma = f->invxsigma ;
        xcenter   = f->xcenter;
        xorder    = f->xorder ;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", xlength: " + std::to_string(xlength)
                         + ", xsigma: " + std::to_string(1./invxsigma)
                         + ", xcenter: " + std::to_string(xcenter)
                         + ", xorder: " + std::to_string(xorder)
                         + ")";
        return info;
    };
private:
    double value, xvacuum, xlength, invxsigma, xcenter;
    int xorder;
};


class Function_Gaussian2D : public Function
{
public:
    Function_Gaussian2D( PyObject *py_profile )
    {
        double xsigma( 0 ), ysigma( 0 );
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xsigma", xsigma );
        PyTools::getAttr( py_profile, "xcenter", xcenter );
        PyTools::getAttr( py_profile, "xorder", xorder );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "ylength", ylength );
        PyTools::getAttr( py_profile, "ysigma", ysigma );
        PyTools::getAttr( py_profile, "ycenter", ycenter );
        PyTools::getAttr( py_profile, "yorder", yorder );
        invxsigma = 1./xsigma;
        invysigma = 1./ysigma;
    };
    Function_Gaussian2D( Function_Gaussian2D *f )
    {
        value     = f->value  ;
        xvacuum   = f->xvacuum;
        xlength   = f->xlength;
        invxsigma = f->invxsigma ;
        xcenter   = f->xcenter;
        xorder    = f->xorder ;
        yvacuum   = f->yvacuum;
        ylength   = f->ylength;
        invysigma = f->invysigma ;
        ycenter   = f->ycenter;
        yorder    = f->yorder ;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", xlength: " + std::to_string(xlength)
                         + ", xsigma: " + std::to_string(1./invxsigma)
                         + ", xcenter: " + std::to_string(xcenter)
                         + ", xorder: " + std::to_string(xorder)
                         + ")";
        return info;
    };
private:
    double value,
           xvacuum, xlength, invxsigma, xcenter,
           yvacuum, ylength, invysigma, ycenter;
    int xorder, yorder;
};


class Function_Gaussian3D : public Function
{
public:
    Function_Gaussian3D( PyObject *py_profile )
    {
        double xsigma, ysigma, zsigma;
        PyTools::getAttr( py_profile, "value", value );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xsigma", xsigma );
        PyTools::getAttr( py_profile, "xcenter", xcenter );
        PyTools::getAttr( py_profile, "xorder", xorder );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "ylength", ylength );
        PyTools::getAttr( py_profile, "ysigma", ysigma );
        PyTools::getAttr( py_profile, "ycenter", ycenter );
        PyTools::getAttr( py_profile, "yorder", yorder );
        PyTools::getAttr( py_profile, "zvacuum", zvacuum );
        PyTools::getAttr( py_profile, "zlength", zlength );
        PyTools::getAttr( py_profile, "zsigma", zsigma );
        PyTools::getAttr( py_profile, "zcenter", zcenter );
        PyTools::getAttr( py_profile, "zorder", zorder );
        invxsigma = 1./xsigma;
        invysigma = 1./ysigma;
        invzsigma = 1./zsigma;
    };
    Function_Gaussian3D( Function_Gaussian3D *f )
    {
        value     = f->value  ;
        xvacuum   = f->xvacuum;
        xlength   = f->xlength;
        invxsigma = f->invxsigma ;
        xcenter   = f->xcenter;
        xorder    = f->xorder ;
        yvacuum   = f->yvacuum;
        ylength   = f->ylength;
        invysigma = f->invysigma ;
        ycenter   = f->ycenter;
        yorder    = f->yorder ;
        zvacuum   = f->zvacuum;
        zlength   = f->zlength;
        invzsigma = f->invzsigma ;
        zcenter   = f->zcenter;
        zorder    = f->zorder ;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (value: " + std::to_string(value)
                         + ", xvacuum: " + std::to_string(xvacuum)
                         + ", xlength: " + std::to_string(xlength)
                         + ", xsigma: " + std::to_string(1./invxsigma)
                         + ", xcenter: " + std::to_string(xcenter)
                         + ", xorder: " + std::to_string(xorder)
                         + ")";
        return info;
    };
private:
    double value,
           xvacuum, xlength, invxsigma, xcenter,
           yvacuum, ylength, invysigma, ycenter,
           zvacuum, zlength, invzsigma, zcenter;
    int xorder, yorder, zorder;
};


class Function_Polygonal1D : public Function
{
public:
    Function_Polygonal1D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "xpoints", xpoints );
        PyTools::getAttr( py_profile, "xvalues", xvalues );
        PyTools::getAttr( py_profile, "xslopes", xslopes );
        npoints = xpoints.size();
    };
    Function_Polygonal1D( Function_Polygonal1D *f )
    {
        xpoints = f->xpoints;
        xvalues = f->xvalues;
        xslopes = f->xslopes;
        npoints = xpoints.size();
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (xpoints: [";
        for(std::vector<double>::iterator point = xpoints.begin(); point != xpoints.end(); ++point) {
            info += " " + std::to_string(*point) ;
        }
        info += " ]";
        info += ", values: [";
        for(std::vector<double>::iterator value = xvalues.begin(); value != xvalues.end(); ++value) {
            info += " " + std::to_string(*value) ;
        }
        info += " ]";
        info += ", xslopes: [";
        for(std::vector<double>::iterator slope = xslopes.begin(); slope != xslopes.end(); ++slope) {
            info += " " + std::to_string(*slope) ;
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    std::vector<double> xpoints, xvalues, xslopes;
    int npoints;
};


class Function_Polygonal2D : public Function
{
public:
    Function_Polygonal2D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "xpoints", xpoints );
        PyTools::getAttr( py_profile, "xvalues", xvalues );
        PyTools::getAttr( py_profile, "xslopes", xslopes );
        npoints = xpoints.size();
    };
    Function_Polygonal2D( Function_Polygonal2D *f )
    {
        xpoints = f->xpoints;
        xvalues = f->xvalues;
        xslopes = f->xslopes;
        npoints = xpoints.size();
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (xpoints: [";
        for(std::vector<double>::iterator point = xpoints.begin(); point != xpoints.end(); ++point) {
            info += " " + std::to_string(*point) ;
        }
        info += " ]";
        info += ", values: [";
        for(std::vector<double>::iterator value = xvalues.begin(); value != xvalues.end(); ++value) {
            info += " " + std::to_string(*value) ;
        }
        info += " ]";
        info += ", xslopes: [";
        for(std::vector<double>::iterator slope = xslopes.begin(); slope != xslopes.end(); ++slope) {
            info += " " + std::to_string(*slope) ;
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    std::vector<double> xpoints, xvalues, xslopes;
    int npoints;
};


class Function_Polygonal3D : public Function
{
public:
    Function_Polygonal3D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "xpoints", xpoints );
        PyTools::getAttr( py_profile, "xvalues", xvalues );
        PyTools::getAttr( py_profile, "xslopes", xslopes );
        npoints = xpoints.size();
    };
    Function_Polygonal3D( Function_Polygonal3D *f )
    {
        xpoints = f->xpoints;
        xvalues = f->xvalues;
        xslopes = f->xslopes;
        npoints = xpoints.size();
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (xpoints: [";
        for(std::vector<double>::iterator point = xpoints.begin(); point != xpoints.end(); ++point) {
            info += " " + std::to_string(*point) ;
        }
        info += " ]";
        info += ", values: [";
        for(std::vector<double>::iterator value = xvalues.begin(); value != xvalues.end(); ++value) {
            info += " " + std::to_string(*value) ;
        }
        info += " ]";
        info += ", xslopes: [";
        for(std::vector<double>::iterator slope = xslopes.begin(); slope != xslopes.end(); ++slope) {
            info += " " + std::to_string(*slope) ;
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    std::vector<double> xpoints, xvalues, xslopes;
    int npoints;
};


class Function_Cosine1D : public Function
{
public:
    Function_Cosine1D( PyObject *py_profile )
    {
        double xlength, xnumber;
        PyTools::getAttr( py_profile, "base", base );
        PyTools::getAttr( py_profile, "xamplitude", xamplitude );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xphi", xphi );
        PyTools::getAttr( py_profile, "xnumber", xnumber );
        invxlength = 1./xlength;
        xnumber2pi = 2.*M_PI*xnumber;
    };
    Function_Cosine1D( Function_Cosine1D *f )
    {
        base       = f->base      ;
        xamplitude = f->xamplitude;
        xvacuum    = f->xvacuum   ;
        invxlength = f->invxlength;
        xphi       = f->xphi      ;
        xnumber2pi = f->xnumber2pi     ;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = "";
        info += " (base: " + std::to_string(base);
        info += ", xamplitude: " + std::to_string(xamplitude);
        info += ", xvacuum: " + std::to_string(xvacuum);
        info += ", xphi: " + std::to_string(xphi);
        info += ", xlength: " + std::to_string(1./invxlength);
        info += ", xnumber: " + std::to_string(xnumber2pi/(2.*M_PI));
        info += ")";
        return info;
    };
private:
    double base, xamplitude, xvacuum, invxlength, xphi, xnumber2pi;
};


class Function_Cosine2D : public Function
{
public:
    Function_Cosine2D( PyObject *py_profile )
    {
        double xlength, xnumber, ylength, ynumber;
        PyTools::getAttr( py_profile, "base", base );
        PyTools::getAttr( py_profile, "xamplitude", xamplitude );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xphi", xphi );
        PyTools::getAttr( py_profile, "xnumber", xnumber );
        PyTools::getAttr( py_profile, "yamplitude", yamplitude );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "ylength", ylength );
        PyTools::getAttr( py_profile, "yphi", yphi );
        PyTools::getAttr( py_profile, "ynumber", ynumber );
        invxlength = 1./xlength;
        xnumber2pi = 2.*M_PI*xnumber;
        invylength = 1./ylength;
        ynumber2pi = 2.*M_PI*ynumber;
    };
    Function_Cosine2D( Function_Cosine2D *f )
    {
        base       = f->base      ;
        xamplitude = f->xamplitude;
        xvacuum    = f->xvacuum   ;
        invxlength = f->invxlength;
        xphi       = f->xphi      ;
        xnumber2pi = f->xnumber2pi;
        yamplitude = f->yamplitude;
        yvacuum    = f->yvacuum   ;
        invylength = f->invylength;
        yphi       = f->yphi      ;
        ynumber2pi = f->ynumber2pi;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = "";
        info += " (base: " + std::to_string(base);
        info += ", xamplitude: " + std::to_string(xamplitude);
        info += ", xvacuum: " + std::to_string(xvacuum);
        info += ", xphi: " + std::to_string(xphi);
        info += ", xlength: " + std::to_string(1./invxlength);
        info += ", xnumber: " + std::to_string(xnumber2pi/(2.*M_PI));
        info += ")";
        return info;
    };
private:
    double base,
           xamplitude, xvacuum, invxlength, xphi, xnumber2pi,
           yamplitude, yvacuum, invylength, yphi, ynumber2pi;
};

class Function_Cosine3D : public Function
{
public:
    Function_Cosine3D( PyObject *py_profile )
    {
        double xlength, xnumber, ylength, ynumber, zlength, znumber;
        PyTools::getAttr( py_profile, "base", base );
        PyTools::getAttr( py_profile, "xamplitude", xamplitude );
        PyTools::getAttr( py_profile, "xvacuum", xvacuum );
        PyTools::getAttr( py_profile, "xlength", xlength );
        PyTools::getAttr( py_profile, "xphi", xphi );
        PyTools::getAttr( py_profile, "xnumber", xnumber );
        PyTools::getAttr( py_profile, "yamplitude", yamplitude );
        PyTools::getAttr( py_profile, "yvacuum", yvacuum );
        PyTools::getAttr( py_profile, "ylength", ylength );
        PyTools::getAttr( py_profile, "yphi", yphi );
        PyTools::getAttr( py_profile, "ynumber", ynumber );
        PyTools::getAttr( py_profile, "zamplitude", zamplitude );
        PyTools::getAttr( py_profile, "zvacuum", zvacuum );
        PyTools::getAttr( py_profile, "zlength", zlength );
        PyTools::getAttr( py_profile, "zphi", zphi );
        PyTools::getAttr( py_profile, "znumber", znumber );
        invxlength = 1./xlength;
        xnumber2pi = 2.*M_PI*xnumber;
        invylength = 1./ylength;
        ynumber2pi = 2.*M_PI*ynumber;
        invzlength = 1./zlength;
        znumber2pi = 2.*M_PI*znumber;
    };
    Function_Cosine3D( Function_Cosine3D *f )
    {
        base       = f->base      ;
        xamplitude = f->xamplitude;
        xvacuum    = f->xvacuum   ;
        invxlength = f->invxlength;
        xphi       = f->xphi      ;
        xnumber2pi = f->xnumber2pi;
        yamplitude = f->yamplitude;
        yvacuum    = f->yvacuum   ;
        invylength = f->invylength;
        yphi       = f->yphi      ;
        ynumber2pi = f->ynumber2pi;
        zamplitude = f->zamplitude;
        zvacuum    = f->zvacuum   ;
        invzlength = f->invzlength;
        zphi       = f->zphi      ;
        znumber2pi = f->znumber2pi;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = "";
        info += " (base: " + std::to_string(base);
        info += ", xamplitude: " + std::to_string(xamplitude);
        info += ", xvacuum: " + std::to_string(xvacuum);
        info += ", xphi: " + std::to_string(xphi);
        info += ", xlength: " + std::to_string(1./invxlength);
        info += ", xnumber: " + std::to_string(xnumber2pi/(2.*M_PI));
        info += ")";
        return info;
    };
private:
    double base,
           xamplitude, xvacuum, invxlength, xphi, xnumber2pi,
           yamplitude, yvacuum, invylength, yphi, ynumber2pi,
           zamplitude, zvacuum, invzlength, zphi, znumber2pi;
};


class Function_Polynomial1D : public Function
{
public:
    Function_Polynomial1D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "orders", orders );
        PyTools::getAttr( py_profile, "coeffs", coeffs );
        PyTools::getAttr( py_profile, "x0", x0 );
        n_orders = orders.size();
    };
    Function_Polynomial1D( Function_Polynomial1D *f )
    {
        orders   = f->orders;
        coeffs   = f->coeffs;
        x0       = f->x0    ;
        n_orders = f->n_orders;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (x0: " + std::to_string(x0) + ", orders: [";
        for(std::vector<unsigned int>::iterator order = orders.begin(); order != orders.end(); ++order) {
            info += " " + std::to_string(*order) ;
        }
        info += " ]";
        info += ", coeffs: [";
        for(unsigned int ix = 0 ; ix < coeffs.size() ; ix++) {
            for(unsigned int iy = 0 ; iy < coeffs[ix].size() ; iy++) {
                info += " " + std::to_string(coeffs[ix][iy]) ;
            }
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    double x0;
    unsigned int n_orders;
    std::vector<unsigned int> orders;
    std::vector<std::vector<double> > coeffs;
};


class Function_Polynomial2D : public Function
{
public:
    Function_Polynomial2D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "orders", orders );
        PyTools::getAttr( py_profile, "coeffs", coeffs );
        PyTools::getAttr( py_profile, "x0", x0 );
        PyTools::getAttr( py_profile, "y0", y0 );
        for( unsigned int i=0; i<orders.size(); i++ )
            if( coeffs[i].size() != orders[i]+1 ) {
                ERROR( "2D polynomial profile has a wrong number of coefficients for order "<<orders[i] );
            }
        n_orders = orders.size();
        n_coeffs = orders.back()+1;
    };
    Function_Polynomial2D( Function_Polynomial2D *f )
    {
        orders   = f->orders;
        coeffs   = f->coeffs;
        x0       = f->x0    ;
        y0       = f->y0    ;
        n_orders = f->n_orders;
        n_coeffs = f->n_coeffs;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (x0: " + std::to_string(x0) + ", y0: " + std::to_string(y0) + ", orders: [";
        for(std::vector<unsigned int>::iterator order = orders.begin(); order != orders.end(); ++order) {
            info += " " + std::to_string(*order) ;
        }
        info += " ]";
        info += ", coeffs: [";
        for(unsigned int ix = 0 ; ix < coeffs.size() ; ix++) {
            for(unsigned int iy = 0 ; iy < coeffs[ix].size() ; iy++) {
                info += " " + std::to_string(coeffs[ix][iy]) ;
            }
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    double x0, y0;
    unsigned int n_orders, n_coeffs;
    std::vector<unsigned int> orders;
    std::vector<std::vector<double> > coeffs;
};


class Function_Polynomial3D : public Function
{
public:
    Function_Polynomial3D( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "orders", orders );
        PyTools::getAttr( py_profile, "coeffs", coeffs );
        PyTools::getAttr( py_profile, "x0", x0 );
        PyTools::getAttr( py_profile, "y0", y0 );
        PyTools::getAttr( py_profile, "z0", z0 );
        for( unsigned int i=0; i<orders.size(); i++ )
            if( coeffs[i].size() != ( orders[i]+1 )*( orders[i]+2 )/2 ) {
                ERROR( "3D polynomial profile has a wrong number of coefficients for order "<<orders[i] );
            }
        n_coeffs = ( orders.back()+1 )*( orders.back()+2 )/2;
        n_orders = orders.size();
    };
    Function_Polynomial3D( Function_Polynomial3D *f )
    {
        orders   = f->orders;
        coeffs   = f->coeffs;
        x0       = f->x0    ;
        y0       = f->y0    ;
        z0       = f->z0    ;
        n_orders = f->n_orders;
        n_coeffs = f->n_coeffs;
    };
    double valueAt( std::vector<double> );
    std::string getInfo ()
    {
        std::string info = " (x0: " + std::to_string(x0)
                         + ", y0: " + std::to_string(y0)
                         + ", z0: " + std::to_string(z0)
                         + ", orders: [";
        for(std::vector<unsigned int>::iterator order = orders.begin(); order != orders.end(); ++order) {
            info += " " + std::to_string(*order) ;
        }
        info += " ]";
        info += ", coeffs: [";
        for(unsigned int ix = 0 ; ix < coeffs.size() ; ix++) {
            for(unsigned int iy = 0 ; iy < coeffs[ix].size() ; iy++) {
                info += " " + std::to_string(coeffs[ix][iy]) ;
            }
        }
        info += " ]";
        info += ")";
        return info;
    };
private:
    double x0, y0, z0;
    unsigned int n_orders, n_coeffs;
    std::vector<unsigned int> orders;
    std::vector<std::vector<double> > coeffs;
};


class Function_TimeConstant : public Function
{
public:
    Function_TimeConstant( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "start", start );
    };
    Function_TimeConstant( Function_TimeConstant *f )
    {
        start = f->start;
    };
    double valueAt( double );
    std::string getInfo ()
    {
        std::string info = " (" + std::to_string(start) + ")";
        return info;
    };
private:
    double start;
};


class Function_TimeTrapezoidal : public Function
{
public:
    Function_TimeTrapezoidal( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "start", start );
        PyTools::getAttr( py_profile, "plateau", plateau );
        PyTools::getAttr( py_profile, "slope1", slope1 );
        PyTools::getAttr( py_profile, "slope2", slope2 );
        invslope1 = 1./slope1;
        invslope2 = 1./slope2;
    };
    Function_TimeTrapezoidal( Function_TimeTrapezoidal *f )
    {
        start   = f->start  ;
        plateau = f->plateau;
        slope1  = f->slope1 ;
        slope2  = f->slope2 ;
        invslope1 = 1./slope1;
        invslope2 = 1./slope2;
    };
    double valueAt( double );
    std::string getInfo ()
    {
        std::string info = " (start: " + std::to_string(start)
                         + ", plateau: " + std::to_string(plateau)
                         + ", slope1: " + std::to_string(slope1)
                         + ", slope2: " + std::to_string(slope2)
                         + ")";
        return info;
    };
private:
    double start, plateau, slope1, slope2, invslope1, invslope2;
};


class Function_TimeGaussian : public Function
{
public:
    Function_TimeGaussian( PyObject *py_profile )
    {
        double duration( 0 ), sigma( 0 );
        PyTools::getAttr( py_profile, "start", start );
        PyTools::getAttr( py_profile, "duration", duration );
        PyTools::getAttr( py_profile, "sigma", sigma );
        PyTools::getAttr( py_profile, "center", center );
        PyTools::getAttr( py_profile, "order", order );
        end = start + duration;
        invsigma = 1./sigma;
    };
    Function_TimeGaussian( Function_TimeGaussian *f )
    {
        start    = f->start   ;
        end      = f->end     ;
        invsigma = f->invsigma;
        center   = f->center  ;
        order    = f->order   ;
    };
    double valueAt( double );
    std::string getInfo ()
    {
        std::string info = " (start: " + std::to_string(start)
                         + ", duration: " + std::to_string(end - start)
                         + ", sigma: " + std::to_string(1./invsigma)
                         + ", center: " + std::to_string(center)
                         + ", order: " + std::to_string(order)
                         + ")";
        return info;
    };
private:
    double start, end, invsigma, center, order;
};


class Function_TimePolygonal : public Function
{
public:
    Function_TimePolygonal( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "points", points );
        PyTools::getAttr( py_profile, "values", values );
        PyTools::getAttr( py_profile, "slopes", slopes );
        npoints = points.size();
    };
    Function_TimePolygonal( Function_TimePolygonal *f )
    {
        points = f->points;
        values = f->values;
        slopes = f->slopes;
        npoints = points.size();
    };
    double valueAt( double );
private:
    std::vector<double> points, values, slopes;
    int npoints;
};


class Function_TimeCosine : public Function
{
public:
    Function_TimeCosine( PyObject *py_profile )
    {
        double duration( 0 );
        PyTools::getAttr( py_profile, "base", base );
        PyTools::getAttr( py_profile, "amplitude", amplitude );
        PyTools::getAttr( py_profile, "start", start );
        PyTools::getAttr( py_profile, "duration", duration );
        PyTools::getAttr( py_profile, "phi", phi );
        PyTools::getAttr( py_profile, "freq", freq );
        end = start + duration;
    };
    Function_TimeCosine( Function_TimeCosine *f )
    {
        base      = f->base     ;
        amplitude = f->amplitude;
        start     = f->start    ;
        end       = f->end      ;
        phi       = f->phi      ;
        freq      = f->freq     ;
    };
    double valueAt( double );
private:
    double base, amplitude, start, end, phi, freq;
};


class Function_TimePolynomial : public Function
{
public:
    Function_TimePolynomial( PyObject *py_profile )
    {
        PyTools::getAttr( py_profile, "orders", orders );
        PyTools::getAttr( py_profile, "coeffs", coeffs );
        PyTools::getAttr( py_profile, "t0", t0 );
    };
    Function_TimePolynomial( Function_TimePolynomial *f )
    {
        orders = f->orders;
        coeffs = f->coeffs;
        t0     = f->t0    ;
    };
    double valueAt( double );
private:
    double t0;
    std::vector<int> orders;
    std::vector<double> coeffs;
};

class Function_TimeSin2Plateau : public Function
{
public:
    Function_TimeSin2Plateau( PyObject *py_profile )
    {
        //double duration;
        PyTools::getAttr( py_profile, "start", start );
        PyTools::getAttr( py_profile, "slope1", slope1 );
        PyTools::getAttr( py_profile, "plateau", plateau );
        PyTools::getAttr( py_profile, "slope2", slope2 );
        end = start + slope1 + plateau + slope2;
    };
    Function_TimeSin2Plateau( Function_TimeSin2Plateau *f )
    {
        start   = f->start;
        slope1    = f->slope1;
        plateau = f->plateau;
        slope2   = f->slope2;
        end     = f->end;
    };
    double valueAt( double );
private:
    double start, slope1, plateau, slope2, end;
};

#endif
