#include <limits>
#include <math.h>
#include "userFunctions.h"

#include "Params.h"

//! inverse error function is taken from NIST
double userFunctions::erfinv (double x)
{
    if (x < -1. || x > 1.)
        return std::numeric_limits<double>::quiet_NaN();

    if (x == 0)
        return 0;

    int sign_x=(x > 0? 1 : -1);

    double r;
    if (x <= 0.686) {
        double x2 = x * x;
        r  = x * (((-0.140543331 * x2 + 0.914624893) * x2 + -1.645349621) * x2 + 0.886226899);
        r /= (((0.012229801 * x2 + -0.329097515) * x2 + 1.442710462) * x2 + -2.118377725) * x2 + 1.;
    } else {
        double y = sqrt (-log((1. - x) / 2.));
        r  = (((1.641345311 * y + 3.429567803) * y + -1.62490649) * y + -1.970840454);
        r /= (1.637067800 * y + 3.543889200) * y + 1.;
    }

    r *= (double)sign_x;
    x *= (double)sign_x;

    r -= (erf(r) - x) / (2. / sqrt (M_PI) * exp(-r*r));

    return r;
}

//! inverse error function is taken from M.B. Giles. 'Approximating the erfinv function'. In GPU Computing Gems, volume 2, Morgan Kaufmann, 2011.
double userFunctions::erfinv2 (double x)
{
    double w, p;
    w = -log((1.0-x)*(1.0+x));

    if ( w < 5.000000 ) {
        w = w - 2.500000;
        p = +2.81022636000e-08      ;
        p = +3.43273939000e-07 + p*w;
        p = -3.52338770000e-06 + p*w;
        p = -4.39150654000e-06 + p*w;
        p = +0.00021858087e+00 + p*w;
        p = -0.00125372503e+00 + p*w;
        p = -0.00417768164e+00 + p*w;
        p = +0.24664072700e+00 + p*w;
        p = +1.50140941000e+00 + p*w;
    } else {
        w = sqrt(w) - 3.000000;
        p = -0.000200214257      ;
        p = +0.000100950558 + p*w;
        p = +0.001349343220 + p*w;
        p = -0.003673428440 + p*w;
        p = +0.005739507730 + p*w;
        p = -0.007622461300 + p*w;
        p = +0.009438870470 + p*w;
        p = +1.001674060000 + p*w;
        p = +2.832976820000 + p*w;
    }
    return p*x;
}

// ----------------------------------------------------------------------------
//! \brief Modified Bessel function of first and second kind
//
//! \details This function computes the first and second kind Bessel function.
//! See http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html
//
//! This function is adapted from the Numerical Recipe.
//
//! \param n Bessel function order
//! \param x input Bessel parameter
//! \param K Value of the Bessel function of second kind
//! \param dK Derivative of K
//! \param I Value of the Bessel function of first kind
//! \param dI Derivative of I
//! \param maxit maximal number of iteration for convergence
//! \param esp epsilon, accuracy threhold for convergence
// ----------------------------------------------------------------------------
void userFunctions::modified_bessel_IK(double n, double x,
        double & I, double & dI,
        double & K, double & dK,
        long maxit, double eps,
        bool warning)
{

    //std::cout << "userFunctions::modified_bessel_IK("
    //          << n << ", " << x << ", " << maxit << ")" << std::endl;

    const double xmin=2.0;
    const double fpmin = 1.e-100;
    double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2;
    double gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl;
    double ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2,xx;

    long i,l,nl;

    // Coefficient for Chebychev
    const double c1[7] = {-1.142022680371168e0,6.5165112670737e-3,
        3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,3.67795e-11,
        -1.356e-13};
    const double c2[8] = {1.843740587300905e0,-7.68528408447867e-2,
        1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,2.423096e-10,
        -1.702e-13,-1.49e-15};

    // Checks
    if (x <= 0.0)
    {
        ERROR("In userFunctions::modified_bessel_IK:\n" << "Argument x is negative in modified_bessel_IK");
    }
    if (n <= 0)
    {
        ERROR("In userFunctions::modified_bessel_IK:\n" << "Argument n is negative in modified_bessel_IK");
    }
    if (maxit <= 0)
    {
        ERROR("In userFunctions::modified_bessel_IK:\n" << "Maximal number of iteration is null in modified_bessel_IK");
    }
    if (eps <= 0)
    {
        ERROR("In userFunctions::modified_bessel_IK:\n" << "Accuracy threshold, epsilon, <= 0 in modified_bessel_IK");
    }
    nl=long(n+0.5);
    xmu=n-nl;
    xmu2=xmu*xmu;
    xi=1.0/x;
    xi2=2.0*xi;
    h=n*xi;
    if (h < fpmin) h=fpmin;
    b=xi2*n;
    d=0.0;
    c=h;

    for (i=0;i<maxit;i++) {
        b += xi2;
        d=1.0/(b+d);
        c=b+1.0/c;
        del=c*d;
        h=del*h;
        if (abs(del-1.0) <= eps) break;
    }

    if (i >= maxit && warning) WARNING("In userFunctions::modified_bessel_IK:\n"
                         << "x too large in modified_bessel_IK; try asymptotic expansion");

    ril=fpmin;
    ripl=h*ril;
    ril1=ril;
    rip1=ripl;
    fact=n*xi;
    for (l=nl-1;l >= 0;l--) {
        ritemp=fact*ril+ripl;
        fact -= xi;
        ripl=fact*ritemp+ril;
        ril=ritemp;
    }

    f=ripl/ril;
    if (x < xmin) {
        x2=0.5*x;
        pimu=M_PI*xmu;
        fact = (abs(pimu) < eps ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        e=xmu*d;
        fact2 = (abs(e) < eps ? 1.0 : sinh(e)/e);
        xx=8.0*xmu*xmu-1.0;
        gam1=chebychev_eval(c1,7,xx);
        gam2=chebychev_eval(c2,8,xx);
        gampl= gam2-xmu*gam1;
        gammi= gam2+xmu*gam1;
        ff=fact*(gam1*cosh(e)+gam2*fact2*d);
        sum=ff;
        e=exp(e);
        p=0.5*e/gampl;
        q=0.5/(e*gammi);
        c=1.0;
        d=x2*x2;
        sum1=p;
        for (i=1;i<=maxit;i++) {
            ff=(i*ff+p+q)/(i*i-xmu2);
            c *= (d/i);
            p /= (i-xmu);
            q /= (i+xmu);
            del=c*ff;
            sum += del;
            del1=c*(p-i*ff);
            sum1 += del1;
            //MESSAGE(i << "/" << maxit << " " << abs(del) << " " << abs(sum) << " " << eps);
            if (abs(del) < abs(sum)*eps) break;
        }
        if (i > maxit && warning) WARNING("In userFunctions::modified_bessel_IK:\n"
                            << "Series failed to converge in modified_bessel_IK");
        rkmu=sum;
        rk1=sum1*xi2;
    } else {
        b=2.0*(1.0+x);
        d=1.0/b;
        h=delh=d;
        q1=0.0;
        q2=1.0;
        a1=0.25-xmu2;
        q=c=a1;
        a = -a1;
        s=1.0+q*delh;
        for (i=1;i<maxit;i++) {
            a -= 2*i;
            c = -a*c/(i+1.0);
            qnew=(q1-b*q2)/a;
            q1=q2;
            q2=qnew;
            q += c*qnew;
            b += 2.0;
            d=1.0/(b+a*d);
            delh=(b*d-1.0)*delh;
            h += delh;
            dels=q*delh;
            s += dels;
            if (abs(dels/s) <= eps) break;
        }
        if (i >= maxit && warning) WARNING("In userFunctions::modified_bessel_IK:\n" << "Failure to converge in cf2 in modified_bessel_IK");
        h=a1*h;
        rkmu=sqrt(M_PI/(2.0*x))*exp(-x)/s;
        rk1=rkmu*(xmu+x+0.5-h)*xi;
    }
    rkmup=xmu*xi*rkmu-rk1;
    rimu=xi/(f*rkmu-rkmup);
    I=(rimu*ril1)/ril;
    dI=(rimu*rip1)/ril;
    for (i=1;i <= nl;i++) {
        rktemp=(xmu+i)*xi2*rk1+rkmu;
        rkmu=rk1;
        rk1=rktemp;
    }
    K=rkmu;
    dK=n*xi*rkmu-rk1;
}


// ----------------------------------------------------------------------------
//! \brief Modified Bessel function of second kind K
//
//! \details This function computes the first and second kind Bessel function.
//! See http://mathworld.wolfram.com/ModifiedBesselFunctionoftheSecondKind.html
//
//! This function is adapted from the Numerical Recipe.
//
//! \param n Bessel function order
//! \param x input Bessel parameter
//! \param maxit maximal number of iteration for convergence
//! \param esp epsilon, accuracy threhold for convergence
// ----------------------------------------------------------------------------
double userFunctions::modified_bessel_K(double n, double x,
        long maxit, double eps,
        bool warning)
{

    //std::cout << "userFunctions::modified_bessel_K("
    //          << n << ", " << x << ", " << maxit << ")" << std::endl;

    const double xmin=2.0;
    const double fpmin = 1.e-100;
    double a,a1,b,c,d,del,del1,delh,dels,e,fact,fact2,ff,gam1,gam2;
    double gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ripl;
    double ritemp,rk1,rkmu,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2,xx;
    // double f,ril1,rimu,rip1,rkmup;

    long i,l,nl;

    // Final result
    double K;

    // Coefficient for Chebychev
    const double c1[7] = {-1.142022680371168e0,6.5165112670737e-3,
        3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,3.67795e-11,
        -1.356e-13};
    const double c2[8] = {1.843740587300905e0,-7.68528408447867e-2,
        1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,2.423096e-10,
        -1.702e-13,-1.49e-15};

    // Checks
    if (x <= 0.0)
    {
        ERROR("Argument x is negative in modified_bessel_K");
    }
    if (n <= 0)
    {
        ERROR("Argument n is negative in modified_bessel_K");
    }
    if (maxit <= 0)
    {
        ERROR("Maximal number of iteration is null in modified_bessel_K");
    }
    if (eps <= 0)
    {
        ERROR("Accuracy threshold, epsilon, <= 0 in modified_bessel_K");
    }
    nl=long(n+0.5);
    xmu=n-nl;
    xmu2=xmu*xmu;
    xi=1.0/x;
    xi2=2.0*xi;
    h=n*xi;
    if (h < fpmin) h=fpmin;
    b=xi2*n;
    d=0.0;
    c=h;

    for (i=0;i<maxit;i++) {
        b += xi2;
        d=1.0/(b+d);
        c=b+1.0/c;
        del=c*d;
        h=del*h;
        if (abs(del-1.0) <= eps) break;
    }

    if (i >= maxit && warning) ERROR("x too large in modified_bessel_K; try asymptotic expansion");

    ril=fpmin;
    ripl=h*ril;
    //ril1=ril;
    //rip1=ripl;
    fact=n*xi;
    for (l=nl-1;l >= 0;l--) {
        ritemp=fact*ril+ripl;
        fact -= xi;
        ripl=fact*ritemp+ril;
        ril=ritemp;
    }

    //f=ripl/ril;
    if (x < xmin) {
        x2=0.5*x;
        pimu=M_PI*xmu;
        fact = (abs(pimu) < eps ? 1.0 : pimu/sin(pimu));
        d = -log(x2);
        e=xmu*d;
        fact2 = (abs(e) < eps ? 1.0 : sinh(e)/e);
        xx=8.0*xmu*xmu-1.0;
        gam1=chebychev_eval(c1,7,xx);
        gam2=chebychev_eval(c2,8,xx);
        gampl= gam2-xmu*gam1;
        gammi= gam2+xmu*gam1;
        ff=fact*(gam1*cosh(e)+gam2*fact2*d);
        sum=ff;
        e=exp(e);
        p=0.5*e/gampl;
        q=0.5/(e*gammi);
        c=1.0;
        d=x2*x2;
        sum1=p;
        for (i=1;i<=maxit;i++) {
            ff=(i*ff+p+q)/(i*i-xmu2);
            c *= (d/i);
            p /= (i-xmu);
            q /= (i+xmu);
            del=c*ff;
            sum += del;
            del1=c*(p-i*ff);
            sum1 += del1;
            if (abs(del) < abs(sum)*eps) break;
        }
        if (i > maxit && warning) WARNING("Series failed to converge in modified_bessel_K");
        rkmu=sum;
        rk1=sum1*xi2;
    } else {
        b=2.0*(1.0+x);
        d=1.0/b;
        h=delh=d;
        q1=0.0;
        q2=1.0;
        a1=0.25-xmu2;
        q=c=a1;
        a = -a1;
        s=1.0+q*delh;
        for (i=1;i<maxit;i++) {
            a -= 2*i;
            c = -a*c/(i+1.0);
            qnew=(q1-b*q2)/a;
            q1=q2;
            q2=qnew;
            q += c*qnew;
            b += 2.0;
            d=1.0/(b+a*d);
            delh=(b*d-1.0)*delh;
            h += delh;
            dels=q*delh;
            s += dels;
            if (abs(dels/s) <= eps) break;
        }
        if (i >= maxit && warning) WARNING("Failure to converge in cf2 in modified_bessel_K");
        h=a1*h;
        rkmu=sqrt(M_PI/(2.0*x))*exp(-x)/s;
        rk1=rkmu*(xmu+x+0.5-h)*xi;
    }
    //rkmup=xmu*xi*rkmu-rk1;
    //rimu=xi/(f*rkmu-rkmup);
    //I=(rimu*ril1)/ril;
    //dI=(rimu*rip1)/ril;
    for (i=1;i <= nl;i++) {
        rktemp=(xmu+i)*xi2*rk1+rkmu;
        rkmu=rk1;
        rk1=rktemp;
    }
    K=rkmu;
    //dK=n*xi*rkmu-rk1;
    return K;
}




// ----------------------------------------------------------------------------
//! \brief Chebychev evaluation adapted from the Numerical Recipes
// ----------------------------------------------------------------------------
double userFunctions::chebychev_eval(const double * c, const int m,
        const double x)
{
    double d=0.0,dd=0.0,sv;
    int j;
    for (j=m-1;j>0;j--) {
        sv=d;
        d=2.*x*d-dd+c[j];
        dd=sv;
    }
    return x*d-dd+0.5*c[0];
}

// ----------------------------------------------------------------------------
//! \brief Computation of the abscissa and weight for the
//! Gauss-Legendre Integration between xmin and xmax.
//
//! \details This code is adpated from the Numerical Recipes.
//
//! \param xmin minimum integration boundary
//! \param xmax maximum integration boundary
//! \param x array of abscissa
//! \param w array of weight for integration
//! \param nb_iterations number of iteration for integration (array size)
//! \param eps accuracy threshold for coef computation
// ----------------------------------------------------------------------------
void userFunctions::gauss_legendre_coef(double xmin,double xmax, double * x,
        double * w, int nb_iterations, double eps)
{

    //std::cout << "userFunctions::gauss_legendre_coef" << std::endl;

    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    // Checks
    if (nb_iterations <= 0)
    {
        ERROR("Number of iteration <= 0 in gauss_legendre_coef");
    }
    if (xmax < xmin)
    {
        ERROR("xmax < xmin in gauss_legendre_coef");
    }
    if (eps <= 0)
    {
        ERROR("accuracy threshold epsilon <= 0 in gauss_legendre_coef");
    }
    // The roots are symmetric, so we only find half of them.
    m=(nb_iterations+1)/2;
    xm=0.5*(xmin+xmax);
    xl=0.5*(xmax-xmin);
    for (i=0;i<=m-1;i++) { /* Loop over the desired roots. */
        z=cos(M_PI*(i+1.0-0.25)/(nb_iterations+0.5));
        /* Starting with the above approximation to the ith root, we enter */
        /* the main loop of refinement by Newton's method.                 */
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=nb_iterations;j++) { /* Recurrence to get Legendre polynomial. */
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            // p1 is now the desired Legendre polynomial. We next compute
            // pp, its derivative, by a standard relation involving also
            // p2, the polynomial of one lower order.
            pp=nb_iterations*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp; /* Newton's method. */
        } while (fabs(z-z1) > eps);
        x[i]=xm-xl*z;      /* Scale the root to the desired interval, */
        x[nb_iterations-1-i]=xm+xl*z;  /* and put in its symmetric counterpart.   */
        w[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight             */
        w[nb_iterations-1-i]=w[i];                 /* and its symmetric counterpart. */
        //std::cout << x[i] << " " << w[i] << " " << x[nb_iterations-1-i]<< std::endl;
    }
    //std::cout << "leaving userFunctions::gauss_legendre_coef" << std::endl;
}


// ----------------------------------------------------------------------------
//! \brief Load repartition in 1d between MPI processes
//
//! \param rank MPI process rank number
//! \param nb_ranks Total number of MPI tasks
//! \param nb_elems Total number of element to be distributed
//! \param imin Index of the first element for rank
//! \param nb_loc_elems Number of element for rank
// ----------------------------------------------------------------------------
void userFunctions::distribute_load_1d(int rank,
        int nb_ranks,
        int nb_elems,
        int& imin,
        int& nb_loc_elems)
{
    // If more ranks than elements,
    // only a part of the processes will work
    if (nb_ranks >= nb_elems)
    {
        if (rank < nb_elems)
        {
            imin = rank;
            nb_loc_elems = 1;
        }
        else
        {
            imin = nb_elems;
            nb_loc_elems = 0;
        }
    }
    else
    {

        int quotient;
        int remainder;

        // Part of the load equally distributed
        quotient = nb_elems/nb_ranks;

        // Remaining load to be distributed after balanced repartition
        remainder = nb_elems%nb_ranks;

        if (rank < remainder)
        {
            imin =  rank*quotient+rank;
            nb_loc_elems = quotient + 1;
        }
        else
        {
            imin = remainder + rank*quotient;
            nb_loc_elems = quotient;
        }
    }
}

// ----------------------------------------------------------------------------
//! \brief Load repartition in 1d between MPI processes.
//! This function returns tables of indexes and length for all rank
//
//! \param nb_ranks Total number of MPI tasks
//! \param nb_elems Total number of element to be distributed
//! \param imin_table Index of the first element for rank
//! \param length_table Number of element for rank
// ----------------------------------------------------------------------------
void userFunctions::distribute_load_1d_table(
        int nb_ranks,
        int nb_elems,
        int * imin_table,
        int * length_table)
{

    // If more ranks than elements,
    // only a part of the processes will work
    if (nb_ranks >= nb_elems)
    {
        #pragma omp simd
        for(int rank = 0 ; rank < nb_elems ; rank ++)
        {
            imin_table[rank] = rank;
            length_table[rank] = 1;
        }
        #pragma omp simd
        for(int rank = nb_elems ; rank < nb_ranks ; rank ++)
        {
            imin_table[rank] = nb_elems;
            length_table[rank] = 0;
        }
    }
    else
    {

        int quotient;
        int remainder;

        // Part of the load equally distributed
        quotient = nb_elems/nb_ranks;

        // Remaining load to be distributed after balanced repartition
        remainder = nb_elems%nb_ranks;

        #pragma omp simd
        for (int rank = 0 ; rank < remainder ; rank ++)
        {
            imin_table[rank] =  rank*quotient+rank;
            length_table[rank] = quotient + 1;
        }
        #pragma omp simd
        for (int rank = remainder ; rank < nb_ranks ; rank ++)
        {
            imin_table[rank] = remainder + rank*quotient;
            length_table[rank] = quotient;
        }
    }
}


// ----------------------------------------------------------------------------
//! \brief Equivalent of the function locate in the Numerical Recipes.
//! This function uses a bijection algorithm in a monotonic double array
//! to find the corresponding index i so that elem is between array[i]
//! and array[i+1].
//
//! \param array array in which to find the value
//! \param elem element to be found
//! \param nb_elem number of elements
// ----------------------------------------------------------------------------
int userFunctions::search_elem_in_array(double * array,
                                        double elem,
                                        int nb_elems)
{
    int imin = 0; // lower bound
    int imax = nb_elems-1; // upper bound
    int imid = 0;

    //std::cerr << " userFunctions::search_elem_in_array " << std::endl;

    if (elem == array[0])
    {
        return 0;
    }
    else if (elem == array[nb_elems-1])
    {
        return nb_elems-2;
    }
    else
    {
        while(imax - imin > 1)
        {
            imid= (imin + imax)/2;
            //imid= (imin + imax)>>1;

            if (elem >= array[imid])
            {
                imin = imid;
                /*std::cerr << "elem >= array[imid]: "
                            << imin << " " << imid << " " << imax << " "
                            << elem << " > " << array[imid]  << std::endl;*/
            }
            else
            {
                imax = imid;
                //std::cerr << "elem < array[imid]: " << imin << " " << imid << " " << imax << " " << elem << " < " << array[imid] << std::endl;
            }
        }
        //std::cerr << " end userFunctions::search_elem_in_array " << std::endl;
        return imin;
    }
}
