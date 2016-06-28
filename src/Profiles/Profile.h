#ifndef Profile_H
#define Profile_H

#include <vector>
#include <string>
#include "SmileiMPI.h"
#include "Tools.h"
#include "PyTools.h"

class Function
{
public:
    //! Default constructor
    Function(){};
    //! Default destructor
    ~Function(){};
    // spatial
    virtual double valueAt(std::vector<double>         ) {
        return 0.;
    };
    // temporal
    virtual double valueAt(double x                    ) {
        // just in case someone uses 1D space profile instead of time
        std::vector<double> v(1);
        v[0] = x;
        return valueAt(v);
    };
    // spatio-temporal
    virtual double valueAt(std::vector<double>, double ) {
        return 0.;
    };
};


//  --------------------------------------------------------------------------------------------------------------------
//! Class Profile
//  --------------------------------------------------------------------------------------------------------------------
class Profile
{
public:
    //! Default constructor
    Profile(PyObject*, unsigned int, std::string);
    //! Cloning constructor
    Profile(Profile*);
    //! Default destructor
    ~Profile();
    
    //! Get the value of the profile at some location (spatial)
    inline double valueAt(std::vector<double> coordinates) {
        return function->valueAt(coordinates);
    };
    //! Get the value of the profile at some location (temporal)
    inline double valueAt(double time) {
        return function->valueAt(time);
    };
    //! Get the value of the profile at some location (spatio-temporal)
    inline double valueAt(std::vector<double> coordinates, double time) {
        return function->valueAt(coordinates, time);
    };
    
    //! Get info on the loaded profile, to be printed later
    inline std::string getInfo() { return info; };
    
    //! Name of the profile, in the case of a built-in profile
    std::string profileName;
    
private:
    //! Object that holds the information on the profile function
    Function * function;
    
    //! String containing some info on the profile
    std::string info;
    
    //! Number of variables for the profile function
    int nvariables;
    
};//END class Profile



// Children classes for python functions

class Function_Python1D : public Function
{
public:
    Function_Python1D(PyObject *pp) : py_profile(pp) {};
    Function_Python1D(Function_Python1D *f) : py_profile(f->py_profile) {};
    double valueAt(double); // time
    double valueAt(std::vector<double>); // space
private:
    PyObject *py_profile;
};


class Function_Python2D : public Function
{
public:
    Function_Python2D(PyObject *pp) : py_profile(pp) {};
    Function_Python2D(Function_Python2D *f) : py_profile(f->py_profile) {};
    double valueAt(std::vector<double>, double); // space + time
    double valueAt(std::vector<double>); // space
private:
    PyObject *py_profile;
};


class Function_Python3D : public Function
{
public:
    Function_Python3D(PyObject *pp) : py_profile(pp) {};
    Function_Python3D(Function_Python3D *f) : py_profile(f->py_profile) {};
    double valueAt(std::vector<double>, double); // space + time
    double valueAt(std::vector<double>); // space
private:
    PyObject *py_profile;
};


class Function_Python4D : public Function
{
public:
    Function_Python4D(PyObject *pp) : py_profile(pp) {};
    Function_Python4D(Function_Python4D *f) : py_profile(f->py_profile) {};
    double valueAt(std::vector<double>, double); // space + time
private:
    PyObject *py_profile;
};



// Children classes for hard-coded functions

class Function_Constant1D : public Function
{
public:
    Function_Constant1D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
    };
    Function_Constant1D ( Function_Constant1D *f ) {
        value   = f->value  ;
        xvacuum = f->xvacuum;
    };
    double valueAt(std::vector<double>);
private:
    double value, xvacuum;
};


class Function_Constant2D : public Function
{
public:
    Function_Constant2D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
        PyTools::getAttr(py_profile, "yvacuum" , yvacuum  );
    };
    Function_Constant2D ( Function_Constant2D *f ) {
        value   = f->value  ;
        xvacuum = f->xvacuum;
        yvacuum = f->yvacuum;
    };
    double valueAt(std::vector<double>);
private:
    double value, xvacuum, yvacuum;
};


class Function_Trapezoidal1D : public Function
{
public:
    Function_Trapezoidal1D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
        PyTools::getAttr(py_profile, "xplateau", xplateau );
        PyTools::getAttr(py_profile, "xslope1" , xslope1  );
        PyTools::getAttr(py_profile, "xslope2" , xslope2  );
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
    };
    Function_Trapezoidal1D ( Function_Trapezoidal1D *f ) {
        value    = f->value   ;
        xvacuum  = f->xvacuum ;
        xplateau = f->xplateau;
        xslope1  = f->xslope1 ;
        xslope2  = f->xslope2 ;
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
    };
    double valueAt(std::vector<double>);
private:
    double value, xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2;
};


class Function_Trapezoidal2D : public Function
{
public:
    Function_Trapezoidal2D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
        PyTools::getAttr(py_profile, "xplateau", xplateau );
        PyTools::getAttr(py_profile, "xslope1" , xslope1  );
        PyTools::getAttr(py_profile, "xslope2" , xslope2  );
        PyTools::getAttr(py_profile, "yvacuum" , yvacuum  );
        PyTools::getAttr(py_profile, "yplateau", yplateau );
        PyTools::getAttr(py_profile, "yslope1" , yslope1  );
        PyTools::getAttr(py_profile, "yslope2" , yslope2  );
        invxslope1 = 1./xslope1;
        invxslope2 = 1./xslope2;
        invyslope1 = 1./yslope1;
        invyslope2 = 1./yslope2;
    };
    Function_Trapezoidal2D ( Function_Trapezoidal2D *f ) {
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
    double valueAt(std::vector<double>);
private:
    double value, 
        xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2,
        yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2;
};


class Function_Gaussian1D : public Function
{
public:
    Function_Gaussian1D ( PyObject *py_profile ) {
        double xsigma;
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
        PyTools::getAttr(py_profile, "xlength" , xlength  );
        PyTools::getAttr(py_profile, "xsigma"  , xsigma   );
        PyTools::getAttr(py_profile, "xcenter" , xcenter  );
        PyTools::getAttr(py_profile, "xorder"  , xorder   );
        invxsigma = 1./xsigma;
    };
    Function_Gaussian1D ( Function_Gaussian1D *f ) {
        value     = f->value  ;
        xvacuum   = f->xvacuum;
        xlength   = f->xlength;
        invxsigma = f->invxsigma ;
        xcenter   = f->xcenter;
        xorder    = f->xorder ;
    };
    double valueAt(std::vector<double>);
private:
    double value, xvacuum, xlength, invxsigma, xcenter;
    int xorder;
};


class Function_Gaussian2D : public Function
{
public:
    Function_Gaussian2D ( PyObject *py_profile ) {
        double xsigma, ysigma;
        PyTools::getAttr(py_profile, "value"   , value    );
        PyTools::getAttr(py_profile, "xvacuum" , xvacuum  );
        PyTools::getAttr(py_profile, "xlength" , xlength  );
        PyTools::getAttr(py_profile, "xsigma"  , xsigma   );
        PyTools::getAttr(py_profile, "xcenter" , xcenter  );
        PyTools::getAttr(py_profile, "xorder"  , xorder   );
        PyTools::getAttr(py_profile, "yvacuum" , yvacuum  );
        PyTools::getAttr(py_profile, "ylength" , ylength  );
        PyTools::getAttr(py_profile, "ysigma"  , ysigma   );
        PyTools::getAttr(py_profile, "ycenter" , ycenter  );
        PyTools::getAttr(py_profile, "yorder"  , yorder   );
        invxsigma = 1./xsigma;
        invysigma = 1./ysigma;
    };
    Function_Gaussian2D ( Function_Gaussian2D *f ) {
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
    double valueAt(std::vector<double>);
private:
    double value, 
        xvacuum, xlength, invxsigma, xcenter,
        yvacuum, ylength, invysigma, ycenter;
    int xorder, yorder;
};


class Function_Polygonal1D : public Function
{
public:
    Function_Polygonal1D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "xpoints" , xpoints );
        PyTools::getAttr(py_profile, "xvalues" , xvalues );
        PyTools::getAttr(py_profile, "xslopes" , xslopes );
        npoints = xpoints.size();
    };
    Function_Polygonal1D ( Function_Polygonal1D *f ) {
        xpoints = f->xpoints;
        xvalues = f->xvalues;
        xslopes = f->xslopes;
        npoints = xpoints.size();
    };
    double valueAt(std::vector<double>);
private:
    std::vector<double> xpoints, xvalues, xslopes;
    int npoints;
};


class Function_Polygonal2D : public Function
{
public:
    Function_Polygonal2D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "xpoints" , xpoints );
        PyTools::getAttr(py_profile, "xvalues" , xvalues );
        PyTools::getAttr(py_profile, "xslopes" , xslopes );
        npoints = xpoints.size();
    };
    Function_Polygonal2D ( Function_Polygonal2D *f ) {
        xpoints = f->xpoints;
        xvalues = f->xvalues;
        xslopes = f->xslopes;
        npoints = xpoints.size();
    };
    double valueAt(std::vector<double>);
private:
    std::vector<double> xpoints, xvalues, xslopes;
    int npoints;
};


class Function_Cosine1D : public Function
{
public:
    Function_Cosine1D ( PyObject *py_profile ) {
        double xlength, xnumber;
        PyTools::getAttr(py_profile, "base"      , base       );
        PyTools::getAttr(py_profile, "xamplitude", xamplitude );
        PyTools::getAttr(py_profile, "xvacuum"   , xvacuum    );
        PyTools::getAttr(py_profile, "xlength"   , xlength    );
        PyTools::getAttr(py_profile, "xphi"      , xphi       );
        PyTools::getAttr(py_profile, "xnumber"   , xnumber    );
        invxlength = 1./xlength;
        xfreq = 2.*M_PI*xnumber*invxlength;
    };
    Function_Cosine1D ( Function_Cosine1D *f ) {
        base       = f->base      ;
        xamplitude = f->xamplitude;
        xvacuum    = f->xvacuum   ;
        invxlength = f->invxlength;
        xphi       = f->xphi      ;
        xfreq      = f->xfreq     ;
    };
    double valueAt(std::vector<double>);
private:
    double base, xamplitude, xvacuum, invxlength, xphi, xfreq;
};


class Function_Cosine2D : public Function
{
public:
    Function_Cosine2D ( PyObject *py_profile ) {
        double xlength, xnumber, ylength, ynumber;
        PyTools::getAttr(py_profile, "base"      , base       );
        PyTools::getAttr(py_profile, "xamplitude", xamplitude );
        PyTools::getAttr(py_profile, "xvacuum"   , xvacuum    );
        PyTools::getAttr(py_profile, "xlength"   , xlength    );
        PyTools::getAttr(py_profile, "xphi"      , xphi       );
        PyTools::getAttr(py_profile, "xnumber"   , xnumber    );
        PyTools::getAttr(py_profile, "yamplitude", yamplitude );
        PyTools::getAttr(py_profile, "yvacuum"   , yvacuum    );
        PyTools::getAttr(py_profile, "ylength"   , ylength    );
        PyTools::getAttr(py_profile, "yphi"      , yphi       );
        PyTools::getAttr(py_profile, "ynumber"   , ynumber    );
        invxlength = 1./xlength;
        xfreq = 2.*M_PI*xnumber*invxlength;
        invylength = 1./ylength;
        yfreq = 2.*M_PI*ynumber*invylength;
    };
    Function_Cosine2D ( Function_Cosine2D *f ) {
        base       = f->base      ;
        xamplitude = f->xamplitude;
        xvacuum    = f->xvacuum   ;
        invxlength = f->invxlength;
        xphi       = f->xphi      ;
        xfreq      = f->xfreq     ;
        yamplitude = f->yamplitude;
        yvacuum    = f->yvacuum   ;
        invylength = f->invylength;
        yphi       = f->yphi      ;
        yfreq      = f->yfreq     ;
    };
    double valueAt(std::vector<double>);
private:
    double base, 
        xamplitude, xvacuum, invxlength, xphi, xfreq,
        yamplitude, yvacuum, invylength, yphi, yfreq;
};


class Function_Polynomial1D : public Function
{
public:
    Function_Polynomial1D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "orders", orders );
        PyTools::getAttr(py_profile, "coeffs", coeffs );
        PyTools::getAttr(py_profile, "x0"    , x0     );
    };
    Function_Polynomial1D ( Function_Polynomial1D *f ) {
        orders = f->orders;
        coeffs = f->coeffs;
        x0     = f->x0    ;
    };
    double valueAt(std::vector<double>);
private:
    double x0;
    std::vector<int> orders;
    std::vector<std::vector<double> > coeffs;
};


class Function_Polynomial2D : public Function
{
public:
    Function_Polynomial2D ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "orders", orders );
        PyTools::getAttr(py_profile, "coeffs", coeffs );
        PyTools::getAttr(py_profile, "x0"    , x0     );
        PyTools::getAttr(py_profile, "y0"    , y0     );
        for( int i=0; i<orders.size(); i++)
            if( coeffs[i].size() != orders[i]+1 )
                ERROR("2D polynomial profile has a wrong number of coefficients for order "<<orders[i]);
    };
    Function_Polynomial2D ( Function_Polynomial2D *f ) {
        orders = f->orders;
        coeffs = f->coeffs;
        x0     = f->x0    ;
        y0     = f->y0    ;
    };
    double valueAt(std::vector<double>);
private:
    double x0, y0;
    std::vector<int> orders;
    std::vector<std::vector<double> > coeffs;
};


class Function_TimeConstant : public Function
{
public:
    Function_TimeConstant ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "start", start);
    };
    Function_TimeConstant ( Function_TimeConstant *f ) {
        start = f->start;
    };
    double valueAt(double);
private:
    double start;
};


class Function_TimeTrapezoidal : public Function
{
public:
    Function_TimeTrapezoidal ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "start"  , start   );
        PyTools::getAttr(py_profile, "plateau", plateau );
        PyTools::getAttr(py_profile, "slope1" , slope1  );
        PyTools::getAttr(py_profile, "slope2" , slope2  );
        invslope1 = 1./slope1;
        invslope2 = 1./slope2;
    };
    Function_TimeTrapezoidal ( Function_TimeTrapezoidal *f ) {
        start   = f->start  ;
        plateau = f->plateau;
        slope1  = f->slope1 ;
        slope2  = f->slope2 ;
        invslope1 = 1./slope1;
        invslope2 = 1./slope2;
    };
    double valueAt(double);
private:
    double start, plateau, slope1, slope2, invslope1, invslope2;
};


class Function_TimeGaussian : public Function
{
public:
    Function_TimeGaussian ( PyObject *py_profile ) {
        double duration, sigma;
        PyTools::getAttr(py_profile, "start"   , start    );
        PyTools::getAttr(py_profile, "duration", duration );
        PyTools::getAttr(py_profile, "sigma"   , sigma    );
        PyTools::getAttr(py_profile, "center"  , center   );
        PyTools::getAttr(py_profile, "order"   , order    );
        end = start + duration;
        invsigma = 1./sigma;
    };
    Function_TimeGaussian ( Function_TimeGaussian *f ) {
        start    = f->start   ;
        end      = f->end     ;
        invsigma = f->invsigma;
        center   = f->center  ;
        order    = f->order   ;
    };
    double valueAt(double);
private:
    double start, end, invsigma, center, order;
};


class Function_TimePolygonal : public Function
{
public:
    Function_TimePolygonal ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "points" , points );
        PyTools::getAttr(py_profile, "values" , values );
        PyTools::getAttr(py_profile, "slopes" , slopes );
        npoints = points.size();
    };
    Function_TimePolygonal ( Function_TimePolygonal *f ) {
        points = f->points;
        values = f->values;
        slopes = f->slopes;
        npoints = points.size();
    };
    double valueAt(double);
private:
    std::vector<double> points, values, slopes;
    int npoints;
};


class Function_TimeCosine : public Function
{
public:
    Function_TimeCosine ( PyObject *py_profile ) {
        double duration;
        PyTools::getAttr(py_profile, "base"     , base      );
        PyTools::getAttr(py_profile, "amplitude", amplitude );
        PyTools::getAttr(py_profile, "start"    , start     );
        PyTools::getAttr(py_profile, "duration" , duration  );
        PyTools::getAttr(py_profile, "phi"      , phi       );
        PyTools::getAttr(py_profile, "freq"     , freq      );
        end = start + duration;
    };
    Function_TimeCosine ( Function_TimeCosine *f ) {
        base      = f->base     ;
        amplitude = f->amplitude;
        start     = f->start    ;
        end       = f->end      ;
        phi       = f->phi      ;
        freq      = f->freq     ;
    };
    double valueAt(double);
private:
    double base, amplitude, start, end, phi, freq;
};


class Function_TimePolynomial : public Function
{
public:
    Function_TimePolynomial ( PyObject *py_profile ) {
        PyTools::getAttr(py_profile, "orders", orders );
        PyTools::getAttr(py_profile, "coeffs", coeffs );
        PyTools::getAttr(py_profile, "t0"    , t0     );
    };
    Function_TimePolynomial ( Function_TimePolynomial *f ) {
        orders = f->orders;
        coeffs = f->coeffs;
        t0     = f->t0    ;
    };
    double valueAt(double);
private:
    double t0;
    std::vector<int> orders;
    std::vector<double> coeffs;
};



#endif
