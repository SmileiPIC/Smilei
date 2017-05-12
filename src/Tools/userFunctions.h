#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

class userFunctions {

    public:

        static double erfinv(double x);
        static double erfinv2(double x);

        //! Modified Bessel function of first and second kind
        static void modified_bessel_IK(double n, double x, double & I, double & dI, 
                double & K, double & dK, long maxit, double eps);

        //! Chebychev evaluation
        static double chebychev_eval(const double * c, const int m, const double x);

        //! Computation of the Gauss-Legendre abscissa and weight 
        static void gauss_legendre_coef(double xmin,double xmax, double * x, 
                double * w, int nbit, double eps);

        //! Load repartition in 1d between MPI processes
        static void distribute_load_1d(int rank,
                int nb_ranks,
                int nb_elems,
                int& imin,
                int& nb_loc_elems);

        //! Load repartition in 1d between MPI processes. 
        //! This function returns tables of indexes and length for all rank
        static void distribute_load_1d_table(
                int nb_ranks,
                int nb_elems,
                int * imin_table,
                int * length_table);

    private:


};
#endif

