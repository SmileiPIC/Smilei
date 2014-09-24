#include "DensityProfile1D.h"


double DensityProfile1D::operator() (PicParams* params, unsigned int ispec, std::vector<double> x_cell) {
    
    SpeciesStructure par=params->species_param[ispec];
    
    // Constant density profile
    // ------------------------
    if (par.plasma_geometry=="constant") {
        if (   (x_cell[0]>par.vacuum_length[0])
            && (x_cell[0]<par.vacuum_length[0]+par.plasma_length[0]) ) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
    
    // Trapezoidal density profile
    // ---------------------------
    else if (par.plasma_geometry=="trap") {
        
        if(par.slope_length.size()!=0){
            
            // vacuum region
            if ( x_cell[0] < par.vacuum_length[0] ) {
                return 0.0;
            }
            // linearly increasing density
            else if ( x_cell[0] < par.vacuum_length[0]+par.slope_length[0] ) {
                return (x_cell[0]-par.vacuum_length[0]) / par.slope_length[0];
            }
            // density plateau
            else if ( x_cell[0] < par.vacuum_length[0]+par.plasma_length[0]-par.slope_length[0] ) {
                return 1.0;
            }
            // linearly decreasing density
            else if ( x_cell[0] < par.vacuum_length[0]+par.plasma_length[0] ) {
                return 1.0 - ( x_cell[0] - (par.vacuum_length[0]+par.plasma_length[0]-par.slope_length[0]) )
                /            par.slope_length[0];
            }
            // beyond the plasma
            else {
                return 0.0;
            }
        }
        else{
            // vacuum region
            if ( x_cell[0] < par.vacuum_length[0] ) {
                return 0.0;
            }
            // linearly increasing density
            else if ( x_cell[0] < par.vacuum_length[0]+par.left_slope_length[0] ) {
                return (x_cell[0]-par.vacuum_length[0]) / par.left_slope_length[0];
            }
            // density plateau
            else if ( x_cell[0] < par.vacuum_length[0]+par.plasma_length[0]-par.right_slope_length[0] ) {
                return 1.0;
            }
            // linearly decreasing density
            else if ( x_cell[0] < par.vacuum_length[0]+par.plasma_length[0] ) {
                return 1.0 - ( x_cell[0] - (par.vacuum_length[0]+par.plasma_length[0]-par.right_slope_length[0]) )
                /            par.right_slope_length[0];
            }
            
            else{
                return 0.0;
            }
        }
        
        
    }
    // Triangular density profile
    // ---------------------------
    else if (par.plasma_geometry=="triangular") {
        
        // vacuum region
        if ( x_cell[0] < par.vacuum_length[0] ) {
            return 0.0;
        }
        // linearly increasing density
        else if ( x_cell[0] < par.vacuum_length[0]+par.left_slope_length[0] ) {
            return (x_cell[0]-par.vacuum_length[0]) / par.left_slope_length[0];
        }
        // linearly decreasing density
        else if ( x_cell[0] < par.vacuum_length[0]+par.plasma_length[0] ) {
            return 1.0 - ( x_cell[0] - (par.vacuum_length[0]+par.plasma_length[0]-par.right_slope_length[0]) )
            /            par.right_slope_length[0];
        }
        
        
        else{
            return 0.0;
        }
        
    }
    // Gaussin density profile
    // ---------------------------
    else if (par.plasma_geometry=="gaussian") {
        
        double lt=par.cut[0]*par.sigma[0];
        // vacuum region
        if ( x_cell[0] < par.vacuum_length[0] ) {
            return 0.0;
        }
        //gaussian
        else if(x_cell[0] < par.vacuum_length[0]+lt ) {
            return exp(-(x_cell[0]-par.vacuum_length[0]-lt)*(x_cell[0]-par.vacuum_length[0]-lt)/(2*par.sigma[0]*par.sigma[0]));
            
        }
        else if(x_cell[0] < par.vacuum_length[0]+lt+par.plateau[0]) {
            return 1.0;
        }
        else if(x_cell[0] < par.vacuum_length[0]+ 2*lt+par.plateau[0]) {
            return exp(-(x_cell[0]-par.vacuum_length[0]-lt-par.plateau[0])*(x_cell[0]-par.vacuum_length[0]-lt-par.plateau[0])/(2*par.sigma[0]*par.sigma[0]));
        }
        
        else{
            return 0.0;
        }
    }
    // Polygonal density profile
    // ---------------------------
    else if (par.plasma_geometry=="polygonal") {
        // vacuum region
        if ( x_cell[0] < par.vacuum_length[0] ) {
            return 0.0;
        }
        else if(x_cell[0] < par.vacuum_length[0]+par.plasma_length[0]){
            for(unsigned int i=1;i<par.x_density_coor.size();++i){
                if (x_cell[0]<par.x_density_coor[i]&&x_cell[0]>=par.x_density_coor[i-1]){
                    double m = (par.density_rel_values_x[i]-par.density_rel_values_x[i-1])/
                    (par.x_density_coor[i]-par.x_density_coor[i-1]);
                    
                    return par.density_rel_values_x[i-1]+m*(x_cell[0]-par.x_density_coor[i-1]);
                }
            }
            
        }
        else{
            return 0.0;
        }
    }
    // cosine density profile
    // ---------------------------
    else if (par.plasma_geometry=="cosine") {
        // vacuum region
        if ( x_cell[0] < par.vacuum_length[0] ) {
            return 0.0;
        }
        else if (x_cell[0] < par.vacuum_length[0]+par.plasma_length[0]) {
            return 1.0+par.ampl*cos(par.mode*(x_cell[0]-par.vacuum_length[0])+par.thetax);
        }
        else{
            return 0.0;
        }
    }
    // odd  species plasma is plasma_length[0] thick and starts at vacuum_length[0]
    // even species plasma is plasma_length[0] thick and starts at vacuum_length[0] + plasma_length[0] + vacuum_length[1]
    else if (par.plasma_geometry=="separated") {
        if ((ispec%2)==0) {                
            if (x_cell[0] > par.vacuum_length[0] && x_cell[0] < par.vacuum_length[0]+ par.plasma_length[0]) {
                return 1.0;
            }
        } else {                
            if (x_cell[0] > par.vacuum_length[0]+ par.plasma_length[0]+par.vacuum_length[1] && x_cell[0] < par.vacuum_length[0]+ 2*par.plasma_length[0] + par.vacuum_length[1]) {
                return 1.0;
            }
        }
        return 0.0;
    }
    // Other density profile
    // ---------------------
    else {
        ERROR("Density profile " << par.plasma_geometry << " not yet defined in 1D");
    }
    
    return 0;
};
