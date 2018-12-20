#ifndef USERFUNCTIONS_H
#define USERFUNCTIONS_H

class userFunctions {

    public:

        static double erfinv(double x);
        static double erfinv2(double x);

        //! Modified Bessel function of first and second kind
        static void modified_bessel_IK(double n, double x, double & I, double & dI,
                double & K, double & dK, long maxit, double eps,
                bool warning = true);

        static double modified_bessel_K(double n, double x,
                long maxit, double eps,
                bool warning = true);

        //! Chebychev evaluation
        static double chebychev_eval(const double * c, const int m, const double x);

        //! Computation of the Gauss-Legendre abscissa and weight
        static void gauss_legendre_coef(double xmin,double xmax, double * x,
                double * w, int nb_iterations, double eps);

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

        //! \brief Equivalent of the function locate in the Numerical Recipes.
        //! This function uses a bijection algorithm in a monotonic double array
        //! to find the corresponding index i so that elem is between array[i]
        //! and array[i+1].
        //
        //! \param array array in which to find the value
        //! \param elem element to be found
        //! \param nb_elem number of elements
        static int search_elem_in_array(double * array,
                                                double elem,
                                                int nb_elems);

    private:


};
#endif
