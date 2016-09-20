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
              ERROR("Profile `"<<name<<"`: constant() profile defined only in 1D, 2D, or 3D");
        
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



