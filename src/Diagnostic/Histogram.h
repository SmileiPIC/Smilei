#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include "Species.h"


// Class for each axis of the particle diags
class HistogramAxis {
public:
    HistogramAxis(std::string, double, double, int, bool, bool, std::vector<double>);
    ~HistogramAxis();
    
    //! Function that goes through the particles and find where they should go in the axis
    virtual void digitize(Species *, std::vector<double>&, unsigned int) {};
    
    //! quantity of the axis (e.g. 'x', 'px', ...)
    std::string type;
    
    //! starting/ending point for the axis binning
    double min, max;
    //! starting/ending point for the axis binning, accounting for logscale
    double actual_min, actual_max;
    //! number of bins for the axis binning
    int nbins;
    
    //! determines whether linear scale or log scale
    bool logscale;
    
    //! determines whether particles beyond min and max are counted in the first and last bin
    bool edge_inclusive;
    
    double coeff;
    
    //! List of coefficients (a,b,c) for a "composite" type of the form "ax+by+cz"
    std::vector<double> coefficients;
};

//! Children classes, corresponding to various types of axes
class HistogramAxis_x         : public HistogramAxis {
    public:
        HistogramAxis_x        (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_y         : public HistogramAxis {
    public:
        HistogramAxis_y        (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_z         : public HistogramAxis {
    public:
        HistogramAxis_z        (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_px        : public HistogramAxis {
    public:
        HistogramAxis_px       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_py        : public HistogramAxis {
    public:
        HistogramAxis_py       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_pz        : public HistogramAxis {
    public:
        HistogramAxis_pz       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_p         : public HistogramAxis {
    public:
        HistogramAxis_p        (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_gamma     : public HistogramAxis {
    public:
        HistogramAxis_gamma    (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_ekin      : public HistogramAxis {
    public:
        HistogramAxis_ekin     (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_vx        : public HistogramAxis {
    public:
        HistogramAxis_vx       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_vy        : public HistogramAxis {
    public:
        HistogramAxis_vy       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_vz        : public HistogramAxis {
    public:
        HistogramAxis_vz       (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_v         : public HistogramAxis {
    public:
        HistogramAxis_v        (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_vperp2    : public HistogramAxis {
    public:
        HistogramAxis_vperp2   (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_charge    : public HistogramAxis {
    public:
        HistogramAxis_charge   (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_chi       : public HistogramAxis {
    public:
        HistogramAxis_chi      (std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};
class HistogramAxis_composite : public HistogramAxis {
    public:
        HistogramAxis_composite(std::string type, double min, double max, int nbins, bool logscale, bool edge_inclusive, std::vector<double> coefficients)
          : HistogramAxis(type, min, max, nbins, logscale, edge_inclusive, coefficients){};
        void digitize(Species *, std::vector<double>&, unsigned int) override;
};


// Class for making a histogram of particle data
class Histogram {
public:
    Histogram(std::vector<HistogramAxis*>);
    ~Histogram();
    
    void digitize(Species *, std::vector<double>&, std::vector<int>&, std::vector<double>&);
    
    std::vector<HistogramAxis*> axes;
    
private:
    virtual void data_filling(Species*, std::vector<double>&, unsigned int) {};
};

//! Children classes, for various manners to fill the histogram
class Histogram_density         : public Histogram {
    public:
        Histogram_density        (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_charge_density  : public Histogram {
    public:
        Histogram_charge_density (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_jx_density      : public Histogram {
    public:
        Histogram_jx_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_jy_density      : public Histogram {
    public:
        Histogram_jy_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_jz_density      : public Histogram {
    public:
        Histogram_jz_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_ekin_density    : public Histogram {
    public:
        Histogram_ekin_density   (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_p_density       : public Histogram {
    public:
        Histogram_p_density      (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_px_density      : public Histogram {
    public:
        Histogram_px_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_py_density      : public Histogram {
    public:
        Histogram_py_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pz_density      : public Histogram {
    public:
        Histogram_pz_density     (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_xx     : public Histogram {
    public:
        Histogram_pressure_xx    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_yy     : public Histogram {
    public:
        Histogram_pressure_yy    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_zz     : public Histogram {
    public:
        Histogram_pressure_zz    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_xy     : public Histogram {
    public:
        Histogram_pressure_xy    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_xz     : public Histogram {
    public:
        Histogram_pressure_xz    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_pressure_yz     : public Histogram {
    public:
        Histogram_pressure_yz    (std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};
class Histogram_ekin_vx_density : public Histogram {
    public:
        Histogram_ekin_vx_density(std::vector<HistogramAxis*>axes) : Histogram(axes){};
        void data_filling(Species *, std::vector<double>&, unsigned int);
};



#endif

