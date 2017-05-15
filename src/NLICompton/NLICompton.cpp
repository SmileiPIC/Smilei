// ----------------------------------------------------------------------------
//! \file NLICompton.cpp
//
//! \brief Nonlinear Inverse Compton Scattering
//
//! \details This header contains the definition of the class NLICompton.
//
// ----------------------------------------------------------------------------

#include "NLICompton.h" 

#include <cstring>
#include <iostream>
#include <fstream>

#include <cmath>

#include "userFunctions.h"

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for NLICompton
// ---------------------------------------------------------------------------------------------------------------------
NLICompton::NLICompton()
{
    Integfochi.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
//
//! Initialization of the parameters for the nonlinear inverse Compton scattering
//
//! \param params Object Params for the parameters from the input script 
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::initParams(Params& params)
{

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    if( PyTools::nComponents("NLICompton") != 0 )
    {
        // Extraction of the parameter from the input file
        PyTools::extract("chipa_integfochi_min", chipa_integfochi_min, "NLICompton");
        PyTools::extract("chipa_integfochi_max", chipa_integfochi_max, "NLICompton");
        PyTools::extract("dim_integfochi", dim_integfochi, "NLICompton");
    }

    // Computation of the normalized Compton wavelength
    norm_lambda_compton = red_planck_cst*params.referenceAngularFrequency_SI/(electron_mass*c_vacuum);       

    // Computation of the factor factor_dNphdt
    factor_dNphdt = sqrt(3.)*fine_struct_cst/(2.*M_PI*norm_lambda_compton);

    // Computation of the factor for the classical radiated power
    factor_cla_rad_power = 2*fine_struct_cst/(3.*norm_lambda_compton);

    // Some additional checks
    if (chipa_integfochi_min >= chipa_integfochi_max)
    {
        ERROR("chipa_integfochi_min (" << chipa_integfochi_min 
                << ") >= chipa_integfochi_max (" << chipa_integfochi_max << ")")
    } 
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the minimum photon quantum parameter for the array xip
//
//! \details Under this value, photon energy is 
//! considered negligible
//
//! \param smpi Object of class SmileiMPI containing MPI properties 
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_chimin_table(SmileiMPI *smpi)
{

    double chipa; // Temporary particle chi value
    // For percentages
    double pct = 0.;
    double dpct = 10.;
    // table load repartition
    int * imin_table;
    int * length_table;
    // Local array
    double * buffer; 
    int err;  // error MPI
    int rank; // Rank number
    int nb_ranks; // Number of ranks
    // timers
    double t0,t1;

    // Get the MPI rank
    rank = smpi->getRank();

    // Get the number of ranks
    nb_ranks = smpi->getSize();

    // Allocation of the array Integfochi
    xip.resize(chipa_xip_dim);

    // Allocation of the table for load repartition
    imin_table = new int[nb_ranks];
    length_table = new int[nb_ranks];

    // Computation of the delta
    chipa_xip_delta = (log10(chipa_xip_max) 
            - log10(chipa_xip_min))/(chipa_xip_dim-1);

    // Load repartition
    userFunctions::distribute_load_1d_table(nb_ranks,
            chipa_xip_dim,
            imin_table,
            length_table);

} 

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the continuous quantum radiated energy during dt
//
//! \param chipa particle quantum parameter
//! \param dt time step
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::norm_rad_energy(double chipa, double dt)
{
    return g_ridgers(chipa)*dt*chipa*chipa;
}


// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the function g of Erber using the Ridgers approximation formulae
//
//! \param chipa particle quantum parameter
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::g_ridgers(double chipa)
{
    return pow(1. + 4.8*(1+chipa)*log10(1. + 1.7*chipa) + 2.44*chipa*chipa,-2./3.);
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Cross Section dNph/dt which is also 
//! the number of photons generated per time unit.
//
//! \param chipa particle quantum parameter
//! \param gfpa particle gamma factor
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_dNphdt(double chipa,double gfpa)
{

    // Log of the particle quantum parameter chipa
    double logchipa;
    double logchipam, logchipap;
    // Index
    int i_chipa;
    // final value
    double dNphdt;

    logchipa = log10(chipa);

    // Lower index for interpolation in the table integfochi
    i_chipa = int(floor(logchipa-chipa_integfochi_min)/delta_chipa_integfochi);

    // If we are not in the table...
    if (i_chipa < 0)
    {
        i_chipa = 0;
        dNphdt = Integfochi[i_chipa];
    }
    else if (i_chipa >= dim_integfochi-1)
    {
        i_chipa = dim_integfochi-2;
        dNphdt = Integfochi[i_chipa];
    }
    else
    {
        // Upper and minor values for linear interpolation
        logchipam = i_chipa*delta_chipa_integfochi + log10(chipa_integfochi_min);
        logchipap = logchipam + delta_chipa_integfochi;

        // Interpolation
        dNphdt = (Integfochi[i_chipa+1]*abs(logchipa-logchipam) + 
                Integfochi[i_chipa]*abs(logchipap - logchipa))/delta_chipa_integfochi;
    }

    return factor_dNphdt*dNphdt*chipa/gfpa;

}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the table values of integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties 
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_integfochi_table(SmileiMPI *smpi)
{

    // Test if an external table exists, we read the table...
    if (Tools::file_exists("tab_integfochi.bin"))
    {

        if (rank==0)
        {

        // Reading of the table file
        std::ifstream file;
        file.open("tab_integfochi.bin",std::ios::binary);

        if (file.is_open())
        {
 
             // Read the header
             file.read((char*)&dim_integfochi,sizeof (dim_integfochi));
             file.read((char*)&chipa_integfochi_min, sizeof (chipa_integfochi_min));
             file.read((char*)&chipa_integfochi_max, sizeof (chipa_integfochi_max));
 
             // Read the table values
             file.read((char*)&Integfochi[0], sizeof (double)*dim_integfochi);
 
             file.close();
        }

        }

        // Allgathering of the table
        MPI_Bcast(&Integfochi[0], dim_integfochi, MPI_DOUBLE, 0, 
           smpi->getGlobalComm());

    }
    // else the table is generated
    else
    {

        double chipa; // Temporary particle chi value
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int * imin_table;
        int * length_table;
        // Local array
        double * buffer; 
        int err;  // error MPI
        int rank; // Rank number
        int nb_ranks; // Number of ranks
        // timers
        double t0,t1;

        // Get the MPI rank
        rank = smpi->getRank();

        // Get the number of ranks
        nb_ranks = smpi->getSize();

        // Allocation of the array Integfochi
        Integfochi.resize(dim_integfochi);

        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];

        // Computation of the delta
        delta_chipa_integfochi = (log10(chipa_integfochi_max) 
                - log10(chipa_integfochi_min))/(dim_integfochi-1);

        // Load repartition
        userFunctions::distribute_load_1d_table(nb_ranks,
                dim_integfochi,
                imin_table,
                length_table);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]];

        MESSAGE("--- Integration F/chipa table:");

        MESSAGE("MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "Rank: " << i << " imin: " << imin_table[i] << " length: " << length_table[i] );
            }
        }

        MESSAGE("Computation:");
        dpct = std::max(dpct,100./length_table[rank]);
        t0 = MPI_Wtime();
        // Loop over the table values
        for(int i = 0 ; i < length_table[rank] ; i++)
        {
            chipa = pow(10.,(imin_table[rank] + i)*delta_chipa_integfochi + log10(chipa_integfochi_min));

            buffer[i] = NLICompton::compute_integfochi(chipa,
                    0.98e-40*chipa,0.98*chipa,400,1e-15);

            //std::cout << rank << " " << buffer[i] << std::endl;

            if (100.*i >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE(i + 1<< "/" << length_table[rank] << " - " << (int)(std::round(pct)) << "%");
            }
        }
        t1 = MPI_Wtime();

        MESSAGE("done in " << (t1 - t0) << "s");

        // Communication of the data
        err = MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
                &Integfochi[0], &length_table[0], &imin_table[0], 
                MPI_DOUBLE, smpi->getGlobalComm());    

        // Print array after communications
        /* 
           if (rank==1) {
           for (int i=0 ; i< dim_integfochi ; i++) 
           {
           std::cout << Integfochi[i] << std::endl;
           }
           }
         */

        // Free memory
        // delete buffer;
        // delete imin_table;
        delete length_table; 

    }

}


// ---------------------------------------------------------------------------------------------------------------------
// Ouput in a file of the table values of integfochi
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::output_integfochi_table(std::string format)
{

    if (format == "ascii") 
    {
        std::ofstream file;
        file.open("tab_integfochi.dat");

        if (file.is_open()) {

            file.precision(12);

            file << "Table Integration F(CHI)/CHI for Nonlinear Compton Scattering \n";

            file << "Dimension - LOG10(chi_e min) - LOG10(chi_e max) \n";

            file << dim_integfochi ;
            file << " " 
                << log10(chipa_integfochi_min) << " " 
                << log10(chipa_integfochi_max) << "\n";;

            // Loop over the table values
            for(int i = 0 ; i < dim_integfochi ; i++)
            {
                file <<  Integfochi[i] << "\n";
            }

            file.close();
        }
    }
    else if (format == "binary")
    {
        std::ofstream file;
        file.open("tab_integfochi.bin",std::ios::binary);

        if (file.is_open()) {

            file.write((char*)&dim_integfochi,sizeof (dim_integfochi));
            file.write((char*)&chipa_integfochi_min, sizeof (chipa_integfochi_min));
            file.write((char*)&chipa_integfochi_max, sizeof (chipa_integfochi_max));

            // Loop over the table values
            for(int i = 0 ; i < dim_integfochi ; i++)
            {
                file.write((char*)&Integfochi[i], sizeof (double));
            }

            file.close();
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! \brief 
//! Compute integration of F/chi between 
//! using Gauss-Legendre for a given chipa value
//
//
//! \param chipa particle (electron for instance) quantum parameter
//! \param chipmin Minimal integration value (photon quantum parameter)
//! \param chipmax Maximal integration value (photon quantum parameter)
//! \param nbit number of points for integration
//! \param eps integration accuracy
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_integfochi(double chipa,
        double chipmin,
        double chipmax,
        int nbit,
        double eps)
{

    //std::cout << "NLICompton::compute_integfochi" << std::endl;

    // Arrays for Gauss-Legendre integration
    double * gauleg_x = new double[nbit];
    double * gauleg_w = new double[nbit];
    // Photon quantum parameter
    double chiph;
    // Integration result
    double integ;
    // Synchrotron emissivity
    double sync_emi;
    // Iterator
    int i;

    // gauss Legendre coefficients
    userFunctions::gauss_legendre_coef(log10(chipmin),log10(chipmax), 
            gauleg_x, gauleg_w, nbit, eps);

    // Integration loop
    integ = 0;
#pragma omp parallel for reduction(+:integ) private(chiph,sync_emi) shared(chipa,gauleg_w,gauleg_x)
    for(i=0 ; i< nbit ; i++)
    {
        chiph = pow(10.,gauleg_x[i]);
        sync_emi = NLICompton::compute_sync_emissivity_ritus(chipa,chiph,200,1e-15);
        integ += gauleg_w[i]*sync_emi*log(10);
    }

    return integ;

}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the synchrotron emissivity following the formulae of Ritus
//
//! \param chipa particle quantum parameter
//! \param chiph photon quantum parameter
//! \param nbit number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_sync_emissivity_ritus(double chipa, 
        double chiph, int nbit, double eps)
{

    //std::cout << "NLICompton::compute_sync_emissivity_ritus" << std::endl;

    // The photon quantum parameter should be below the electron one
    if (chipa > chiph)
    {
        // Arrays for Gauss-Legendre integration
        double * gauleg_x = new double[nbit];
        double * gauleg_w = new double[nbit];
        // Values for Bessel results       
        double I,dI;
        double K,dK;
        // Parts of the formulae
        double part1,part2;
        // Iterator
        int i;

        double y = chiph/(3.*chipa*(chipa-chiph));

        // Computation of Part. 1
        // Call the modified Bessel function to get K
        userFunctions::modified_bessel_IK(2./3.,2*y,I,dI,K,dK,50000,eps);

        part1 = (2. + 3.*chiph*y)*(K);

        // Computation of Part. 2
        // Using Gauss Legendre integration 

        userFunctions::gauss_legendre_coef(log10(2*y),log10(y)+50., gauleg_x, 
                gauleg_w, nbit, eps);

        part2 = 0;
        for(i=0 ; i< nbit ; i++)
        {
            y = pow(10.,gauleg_x[i]);
            userFunctions::modified_bessel_IK(1./3.,y,I,dI,K,dK,50000,eps);
            part2 += gauleg_w[i]*K*y*log(10.);
        }

        // Factor for final result
        y = 2*chiph/(3*pow(chipa,2.));

        return (part1 - part2)*y;


    }
    else if (chipa == chiph)
    {
        return 0;
    }
    else
    {
        ERROR("In compute_sync_emissivity_ritus: chipa " << chipa << " < chiph " << chiph);
    }    
}
