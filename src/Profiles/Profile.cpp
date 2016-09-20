#include <cmath>

#include "ElectroMagn.h"
#include "Profile.h"
#include "PyTools.h"

using namespace std;



// Default constructor.
Profile::Profile(PyObject* py_profile, unsigned int nvariables, string name) :
    profileName(""),
    nvariables(nvariables)
{
    ostringstream info_("");
    info_ << nvariables << "D";
    
    if (!PyCallable_Check(py_profile)) {
        ERROR("Profile `"<<name<<"`: not a function");
    }
    
    // In case the function was created in "pyprofiles.py", then we transform it
    //  in a "hard-coded" function
    if( PyObject_HasAttrString(py_profile, "profileName") ) {
        
        PyTools::getAttr(py_profile, "profileName", profileName );
        
        info_ << " built-in profile `" << profileName << "`" ;
        
        if( profileName == "constant" ) {
        
            if     ( nvariables == 1 )
                function = new Function_Constant1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Constant2D(py_profile);
            else if( nvariables == 3 )
                function = new Function_Constant3D(py_profile);
            else
              ERROR("Profile `"<<name<<"`: constant() profile defined only in 1D or 2D");
        
        } else if( profileName == "trapezoidal" ){
        
            if     ( nvariables == 1 )
                function = new Function_Trapezoidal1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Trapezoidal2D(py_profile);
            else
              ERROR("Profile `"<<name<<"`: trapezoidal() profile defined only in 1D or 2D");
        
        } else if( profileName == "gaussian" ){
        
            if     ( nvariables == 1 )
                function = new Function_Gaussian1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Gaussian2D(py_profile);
            else
                ERROR("Profile `"<<name<<"`: gaussian() profile defined only in 1D or 2D");
        
        } else if( profileName == "polygonal" ){
        
            if     ( nvariables == 1 )
                function = new Function_Polygonal1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Polygonal2D(py_profile);
            else
                ERROR("Profile `"<<name<<"`: polygonal() profile defined only in 1D or 2D");
                
        } else if( profileName == "cosine" ){
        
            if     ( nvariables == 1 )
                function = new Function_Cosine1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Cosine2D(py_profile);
            else
                ERROR("Profile `"<<name<<"`: cosine() profile defined only in 1D or 2D");
                
        } else if( profileName == "polynomial" ){
        
            if     ( nvariables == 1 )
                function = new Function_Polynomial1D(py_profile);
            else if( nvariables == 2 )
                function = new Function_Polynomial2D(py_profile);
            else
                ERROR("Profile `"<<name<<"`: polynomial() profile defined only in 1D or 2D");
            
        } else if( profileName == "tconstant" ){
        
            if( nvariables == 1 )
                function = new Function_TimeConstant(py_profile);
            else
                ERROR("Profile `"<<name<<"`: tconstant() profile is only for time");
            
        } else if( profileName == "ttrapezoidal" ){
        
            if( nvariables == 1 )
                function = new Function_TimeTrapezoidal(py_profile);
            else
                ERROR("Profile `"<<name<<"`: ttrapezoidal() profile is only for time");
            
        } else if( profileName == "tgaussian" ){
        
            if( nvariables == 1 )
                function = new Function_TimeGaussian(py_profile);
            else
                ERROR("Profile `"<<name<<"`: tgaussian() profile is only for time");
            
        } else if( profileName == "tpolygonal" ){
        
            if( nvariables == 1 )
                function = new Function_TimePolygonal(py_profile);
            else
                ERROR("Profile `"<<name<<"`: tpolygonal() profile is only for time");
            
        } else if( profileName == "tcosine" ){
        
            if( nvariables == 1 )
                function = new Function_TimeCosine(py_profile);
            else
                ERROR("Profile `"<<name<<"`: tcosine() profile is only for time");
            
        } else if( profileName == "tpolynomial" ){
        
            if( nvariables == 1 )
                function = new Function_TimePolynomial(py_profile);
            else
                ERROR("Profile `"<<name<<"`: tpolynomial() profile is only for time");
            
        }
        
    }
    
    // Otherwise (if the python profile cannot be hard-coded) ....
    else {
        string message;
        
#ifdef  __DEBUG
        // Check how the profile looks like
        PyObject* repr = PyObject_Repr(py_profile);
        PyTools::convert(repr, message);
        MESSAGE(message);
        Py_XDECREF(repr);
        
        repr = PyObject_Str(py_profile);
        PyTools::convert(repr, message);
        MESSAGE(message);
        Py_XDECREF(repr);
#endif
        
        // Verify that the profile has the right number of arguments
        PyObject* inspect=PyImport_ImportModule("inspect");
        PyTools::checkPyError();
        PyObject *tuple = PyObject_CallMethod(inspect,const_cast<char *>("getargspec"),const_cast<char *>("(O)"),py_profile);
        PyObject *arglist = PyTuple_GetItem(tuple,0);
        int size = PyObject_Size(arglist);
        if (size != (int)nvariables) {
            string args("");
            for (int i=0; i<size; i++){
                PyObject *arg=PyList_GetItem(arglist,i);
                PyObject* repr = PyObject_Repr(arg);
                PyTools::convert(repr, message);
                args += message+" ";
                Py_XDECREF(repr);
            }
            WARNING ("Profile " << name << " takes "<< size <<" variables (" << args << ") but it is created with " << nvariables);
        }
        Py_XDECREF(tuple);
        Py_XDECREF(inspect);
        
        // Assign the evaluating function, which depends on the number of arguments
        if      ( nvariables == 1 ) function = new Function_Python1D(py_profile);
        else if ( nvariables == 2 ) function = new Function_Python2D(py_profile);
        else if ( nvariables == 3 ) function = new Function_Python3D(py_profile);
        else {
            ERROR("Profile `"<<name<<"`: defined with unsupported number of variables");
        }
        
        info_ << " user-defined function";
    }
    
    info = info_.str();
}


// Cloning constructor
Profile::Profile(Profile *p)
{
    profileName = p->profileName;
    nvariables  = p->nvariables ;
    info        = p->info       ;
    if( profileName != "" ) {
        if( profileName == "constant" ) {
            if     ( nvariables == 1 )
                function = new Function_Constant1D(static_cast<Function_Constant1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Constant2D(static_cast<Function_Constant2D*>(p->function));
            else if( nvariables == 3 )
                function = new Function_Constant3D(static_cast<Function_Constant3D*>(p->function));
        } else if( profileName == "trapezoidal" ){
            if     ( nvariables == 1 )
                function = new Function_Trapezoidal1D(static_cast<Function_Trapezoidal1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Trapezoidal2D(static_cast<Function_Trapezoidal2D*>(p->function));
        } else if( profileName == "gaussian" ){
            if     ( nvariables == 1 )
                function = new Function_Gaussian1D(static_cast<Function_Gaussian1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Gaussian2D(static_cast<Function_Gaussian2D*>(p->function));
        } else if( profileName == "polygonal" ){
            if     ( nvariables == 1 )
                function = new Function_Polygonal1D(static_cast<Function_Polygonal1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Polygonal2D(static_cast<Function_Polygonal2D*>(p->function));
        } else if( profileName == "cosine" ){
            if     ( nvariables == 1 )
                function = new Function_Cosine1D(static_cast<Function_Cosine1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Cosine2D(static_cast<Function_Cosine2D*>(p->function));
        } else if( profileName == "polynomial" ){
            if     ( nvariables == 1 )
                function = new Function_Polynomial1D(static_cast<Function_Polynomial1D*>(p->function));
            else if( nvariables == 2 )
                function = new Function_Polynomial2D(static_cast<Function_Polynomial2D*>(p->function));
        } else if( profileName == "tconstant" ){
            function = new Function_TimeConstant(static_cast<Function_TimeConstant*>(p->function));
        } else if( profileName == "ttrapezoidal" ){
            function = new Function_TimeTrapezoidal(static_cast<Function_TimeTrapezoidal*>(p->function));
        } else if( profileName == "tgaussian" ){
            function = new Function_TimeGaussian(static_cast<Function_TimeGaussian*>(p->function));
        } else if( profileName == "tpolygonal" ){
            function = new Function_TimePolygonal(static_cast<Function_TimePolygonal*>(p->function));
        } else if( profileName == "tcosine" ){
            function = new Function_TimeCosine(static_cast<Function_TimeCosine*>(p->function));
        } else if( profileName == "tpolynomial" ){
            function = new Function_TimePolynomial(static_cast<Function_TimePolynomial*>(p->function));
        }
    } else {
        if      ( nvariables == 1 ) function = new Function_Python1D(static_cast<Function_Python1D*>(p->function));
        else if ( nvariables == 2 ) function = new Function_Python2D(static_cast<Function_Python2D*>(p->function));
        else if ( nvariables == 3 ) function = new Function_Python3D(static_cast<Function_Python3D*>(p->function));
    }
}


Profile::~Profile()
{
    delete function;
}


// Functions to evaluate a python function with various numbers of arguments
// 1D
double Function_Python1D::valueAt(double time) {
    return PyTools::runPyFunction(py_profile, time);
}
double Function_Python1D::valueAt(vector<double> x_cell) {
    return PyTools::runPyFunction(py_profile, x_cell[0]);
}
// 2D
double Function_Python2D::valueAt(vector<double> x_cell, double time) {
    return PyTools::runPyFunction(py_profile, x_cell[0], time);
}
double Function_Python2D::valueAt(vector<double> x_cell) {
    return PyTools::runPyFunction(py_profile, x_cell[0], x_cell[1]);
}
// 3D
double Function_Python3D::valueAt(vector<double> x_cell, double time) {
    return PyTools::runPyFunction(py_profile, x_cell[0], x_cell[1], time);
}
double Function_Python3D::valueAt(vector<double> x_cell) {
    return PyTools::runPyFunction(py_profile, x_cell[0], x_cell[1], x_cell[2]);
}
// 4D
double Function_Python4D::valueAt(vector<double> x_cell, double time) {
    return PyTools::runPyFunction(py_profile, x_cell[0], x_cell[1], x_cell[2], time);
}

// Constant profiles
double Function_Constant1D::valueAt(vector<double> x_cell) {
    return (x_cell[0]>xvacuum) ? value : 0.;
}
double Function_Constant2D::valueAt(vector<double> x_cell) {
    return ((x_cell[0]>xvacuum) && (x_cell[1]>yvacuum)) ? value : 0.;
}
double Function_Constant3D::valueAt(vector<double> x_cell) {
    return ((x_cell[0]>xvacuum) && (x_cell[1]>yvacuum) && (x_cell[2]>zvacuum)) ? value : 0.;
}

// Trapezoidal profiles
inline double trapeze(double x, double plateau, double slope1, double slope2, double invslope1, double invslope2) {
    double result = 0.;
    if ( x > 0. ) {
        if ( x < slope1 ) {
            result = invslope1 * x;
        } else {
            x -= slope1 + plateau;
            if ( x < 0. ) {
                result = 1.;
            } else {
                x -= slope2;
                if ( x < 0. ) result = 1. - x * invslope2;
            }
        }
    }
    return result;
}
double Function_Trapezoidal1D::valueAt(vector<double> x_cell) {
    return value * trapeze(x_cell[0]-xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2);
}
double Function_Trapezoidal2D::valueAt(vector<double> x_cell) {
    return value
        * trapeze(x_cell[0]-xvacuum, xplateau, xslope1, xslope2, invxslope1, invxslope2)
        * trapeze(x_cell[1]-yvacuum, yplateau, yslope1, yslope2, invyslope1, invyslope2);
}

// Gaussian profiles
double Function_Gaussian1D::valueAt(vector<double> x_cell) {
    double x = x_cell[0], xfactor=0.;
    if ( x > xvacuum  && x < xvacuum+xlength )
        xfactor = exp( -pow(x-xcenter, xorder) * invxsigma );
    return value * xfactor;
}
double Function_Gaussian2D::valueAt(vector<double> x_cell) {
    double x = x_cell[0], xfactor=0.;
    double y = x_cell[1], yfactor=0.;
    if ( x > xvacuum  && x < xvacuum+xlength )
        xfactor = exp( -pow(x-xcenter, xorder) * invxsigma );
    if ( y > yvacuum  && y < yvacuum+ylength )
        yfactor = exp( -pow(y-ycenter, yorder) * invysigma );
    return value * xfactor * yfactor;
}

// Polygonal profiles
double Function_Polygonal1D::valueAt(vector<double> x_cell) {
    double x = x_cell[0];
    if( x < xpoints[0] ) return 0.;
    for( int i=1; i<npoints; i++ )
        if( x < xpoints[i] )
            return xvalues[i-1] + xslopes[i-1] * ( x - xpoints[i-1] );
    return 0.;
}
double Function_Polygonal2D::valueAt(vector<double> x_cell) {
    double x = x_cell[0];
    if( x < xpoints[0] ) return 0.;
    for( int i=1; i<npoints; i++ )
        if( x < xpoints[i] )
            return xvalues[i-1] + xslopes[i-1] * ( x - xpoints[i-1] );
    return 0.;
}

// Cosine profiles
double Function_Cosine1D::valueAt(vector<double> x_cell) {
    double x = (x_cell[0] - xvacuum) * invxlength, xfactor = 0.;
    if( x > 0. && x < 1. )
        xfactor = base + xamplitude * cos(xphi + xfreq * x);
    return xfactor;
}
double Function_Cosine2D::valueAt(vector<double> x_cell) {
    double x = (x_cell[0] - xvacuum) * invxlength, xfactor = 0.;
    double y = (x_cell[1] - yvacuum) * invylength, yfactor = 0.;
    if( x > 0. && x < 1. )
        xfactor = base + xamplitude * cos(xphi + xfreq * x);
    if( y > 0. && y < 1. )
        yfactor = base + yamplitude * cos(yphi + yfreq * y);
    return xfactor * yfactor;
}

// Polynomial profiles
double Function_Polynomial1D::valueAt(vector<double> x_cell) {
    double r = 0., xx0 = x_cell[0]-x0, xx = 1.;
    unsigned int currentOrder = 0;
    for( unsigned int i=0; i<orders.size(); i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            xx *= xx0;
        }
        r += coeffs[i][0] * xx;
    }
    return r;
}
double Function_Polynomial2D::valueAt(vector<double> x_cell) {
    double r = 0., xx0 = x_cell[0]-x0, yy0 = x_cell[1]-y0;
    vector<double> xx;
    unsigned int currentOrder = 0;
    xx.resize(orders.back()+1);
    xx[0] = 1.;
    for( unsigned int i=0; i<orders.size(); i++ ) {
        while( currentOrder<orders[i] ) {
            currentOrder += 1;
            xx[currentOrder] = xx[currentOrder-1] * yy0;
            for( unsigned int j=0; j<currentOrder; j++ ) xx[j] *= xx0;
        }
        for( unsigned int j=0; j<=orders[i]; j++ ) r += coeffs[i][j] * xx[j];
    }
    return r;
}

// Time constant profile
double Function_TimeConstant::valueAt(double time) {
    if( time > start ) return 1.;
    else               return 0.;
}

// Time trapezoidal profile
double Function_TimeTrapezoidal::valueAt(double time) {
    return trapeze(time-start, plateau, slope1, slope2, invslope1, invslope2);
}

// Time gaussian profile
double Function_TimeGaussian::valueAt(double time) {
    if( time < start || time > end) return 0.;
    else                            return exp( -pow(time-center,order) * invsigma );
}

// Time polygonal profile
double Function_TimePolygonal::valueAt(double time) {
    if( time < points[0] ) return 0.;
    for( int i=1; i<npoints; i++ )
        if( time<points[i] )
            return values[i-1] + slopes[i-1] * ( time-points[i-1] );
    return 0.;
}

// Time cosine profile
double Function_TimeCosine::valueAt(double time) {
    if( time < start || time > end ) return 0.;
    else                             return base + amplitude * cos(phi + freq*(time-start));
}

// Time polynomial profile
double Function_TimePolynomial::valueAt(double time) {
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


