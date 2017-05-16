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
        PyTools::extract("integfochi_dim", dim_integfochi, "NLICompton");
        PyTools::extract("chipa_xip_min", chipa_xip_min, "NLICompton");
        PyTools::extract("chipa_xip_max", chipa_xip_max, "NLICompton");
        PyTools::extract("xip_power", xip_power, "NLICompton");
        PyTools::extract("xip_threshold", xip_threshold, "NLICompton");
        PyTools::extract("chipa_xip_dim", chipa_xip_dim, "NLICompton");
        PyTools::extract("chiph_xip_dim", chiph_xip_dim, "NLICompton");
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
//! Computation of the minimum photon quantum parameter chiphmin for the xip array
//! and computation of the xip array.
//
//! \details Under the minimum chiph value, photon energy is 
//! considered negligible.
//
//! \param smpi Object of class SmileiMPI containing MPI properties 
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_xip_table(SmileiMPI *smpi)
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0,t1;

    // Get the MPI rank
    rank = smpi->getRank();

    t0 = MPI_Wtime();

    MESSAGE("        --- Table chiphmin for xip");

    // Test if an external table exists, we read the table...
    if (Tools::file_exists("tab_xip.bin"))
    {

        if (rank==0)
        {

            // Reading of the table file
            std::ifstream file;
            file.open("tab_xip.bin",std::ios::binary);

            if (file.is_open())
            {

                // Read the header
                file.read((char*)&chipa_xip_dim,sizeof (chipa_xip_dim));
                file.read((char*)&chiph_xip_dim,sizeof (chiph_xip_dim));
                file.read((char*)&chipa_xip_min, sizeof (chipa_xip_min));
                file.read((char*)&chipa_xip_max, sizeof (chipa_xip_max));

                // Allocation of the array xip
                xip_chiphmin_table.resize(chipa_xip_dim);
                xip_table.resize(chipa_xip_dim);

                // Read the table values
                file.read((char*)&xip_chiphmin_table[0], sizeof (double)*chipa_xip_dim);
                
                // Read the table values
                file.read((char*)&xip_table[0], sizeof (double)*chipa_xip_dim*chiph_xip_dim);

                file.close();
            }

        }

        if (rank != 0)
        { 
            xip_chiphmin_table.resize(chipa_xip_dim);
            xip_table.resize(chipa_xip_dim*chiph_xip_dim);
        }

        // Bcast of the table xip_chiphmin
        MPI_Bcast(&xip_chiphmin_table[0], chipa_xip_dim, MPI_DOUBLE, 0, 
           smpi->getGlobalComm());

        // Bcast of the table xip
        MPI_Bcast(&xip_table[0], chipa_xip_dim*chiph_xip_dim, MPI_DOUBLE, 0, 
           smpi->getGlobalComm());
    }
    // else the table is generated
    else
    {
        // Parameters:
        double chipa; // Temporary particle chi value
        double chiph; // Temporary photon chi value
        double chiph_delta; // Temporary delta for chiph
        double logchiphmin; // Temporary log10 of photon chi value
        double xip;   // Temporary xip
        double numerator;
        double denominator;
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int * imin_table;
        int * length_table;
        // Local array
        double * buffer; 
        //int err;  // error MPI
        int nb_ranks; // Number of ranks
        // Iterator
        int  k;

        // Get the number of ranks
        nb_ranks = smpi->getSize();

        // Allocation of the array xip_chiphmin_table
        xip_chiphmin_table.resize(chipa_xip_dim);

        // Allocation of the array xip_table
        xip_table.resize(chipa_xip_dim*chiph_xip_dim);

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

        // Allocation of the local buffer
        buffer = new double [length_table[rank]];

        MESSAGE("        MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "        Rank: " << i 
                                          << " imin: "   << imin_table[i] 
                                          << " length: " << length_table[i] );
            }
        }

        // 1. - Computation of xip_chiphmin_table
        MESSAGE("        Computation of chiphmin:");
        dpct = std::max(dpct,100./length_table[rank]);
        
        // Loop for chiphmin
        for(int i = 0 ; i < length_table[rank] ; i++)
        {

            xip = 1;
            chipa = pow(10.,(imin_table[rank] + i)*chipa_xip_delta + log10(chipa_xip_min));
            logchiphmin = log10(chipa);

            // Denominator of xip
            denominator = NLICompton::compute_integfochi(chipa,
                    0.98e-40*chipa,0.98*chipa,200,1e-13);

            k = 0;
            while(k < xip_power)
            {
                logchiphmin -= pow(0.1,k);
                chiph = pow(10.,logchiphmin);
                numerator = NLICompton::compute_integfochi(chipa,
                    1e-40*chiph,chiph,200,1e-13);

                if (numerator == 0)
                {
                    xip = 0;
                }
                else
                {
                    xip = numerator/denominator;
                }

                if (xip < xip_threshold)
                {
                    logchiphmin += pow(0.1,k);
                    k += 1;
                }
            }
            buffer[i] = pow(10.,logchiphmin);

            // display percentage
            if (100.*ichipa >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE( "        " << ichipa + 1 << "/" << length_table[rank] 
                                    << " - " << (int)(std::round(pct)) << "%");
            }
        }

        // Communication of the xip_chiphmin table
        MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
                &xip_chiphmin_table[0], &length_table[0], &imin_table[0], 
                MPI_DOUBLE, smpi->getGlobalComm());    

        // 2. - Computation of the xip table
        MESSAGE("        Computation of xip:");
        dpct = std::max(dpct,100./length_table[rank]);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]*chiph_xip_dim];

        // Loop for xip in the chipa dimension
        for(int ichipa = 0 ; ichipa < length_table[rank] ; ichipa++)
        {
 
            chipa = pow(10.,(imin_table[rank] + ichipa)*chipa_xip_delta + log10(chipa_xip_min));

            chiph_delta = (log10(chipa) - log10(xip_chiphmin_table[imin_table[rank] + ichipa]))
                        / (chiph_xip_dim - 1);

            // Denominator of xip
            denominator = NLICompton::compute_integfochi(chipa,
                    1e-40*chipa,chipa,250,1e-15);

            // Loop in the chiph dimension
            for (int ichiph = 0 ; ichiph < chiph_xip_dim ; ichiph ++)
            {
                // Local chiph value
               chiph = pow(10.,ichiph*chiph_delta + log10(xip_chiphmin_table[imin_table[rank] + ichipa])); 

               /* std::cout << "rank: " << rank 
                         << " " << chipa 
                         << " " << chiph 
                         << " " << xip_chiphmin_table[imin_table[rank] + ichipa] << std::endl; */

               // Denominator of xip
               numerator = NLICompton::compute_integfochi(chipa,
                       1e-40*chiph,chiph,250,1e-15);

               // Update local buffer value
               buffer[ichipa*chiph_xip_dim + ichiph] = std::min(1.,numerator / denominator);

            }


            // display percentage
            if (100.*ichipa >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE( "        " << ichipa + 1 << "/" << length_table[rank] 
                                    << " - " << (int)(std::round(pct)) << "%");
            }
        }

        // Update length_table and imin_table
        for (int i = 0 ; i < nb_ranks ; i++)
        {
            length_table[i] *= chiph_xip_dim;
            imin_table[i] *= chiph_xip_dim;
        }

        // Communication of the xip_chiphmin table
        MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
                &xip_table[0], &length_table[0], &imin_table[0], 
                MPI_DOUBLE, smpi->getGlobalComm());    

        delete buffer;
        delete length_table;
        delete imin_table;

    }

    t1 = MPI_Wtime();
    MESSAGE("        done in " << (t1 - t0) << "s");

} 

// ---------------------------------------------------------------------------------------------------------------------
//! File output of xip_chiphmin_table and xip_table
//
//! \param format output format (ascii, binary)
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::output_xip_table(std::string format)
{

    if (format == "ascii") 
    {
        std::ofstream file;
        file.open("tab_xip.dat");

        if (file.is_open()) {

            file.precision(12);

            file << "Table xip_chiphmin and xip for Nonlinear Compton Scattering \n";

            file << "Dimension chipa - Dimension chiph - LOG10(chipa min) - LOG10(chipa max) \n";

            file << chipa_xip_dim << " "
                 << chiph_xip_dim << " "
                 << log10(chipa_xip_min) << " " 
                 << log10(chipa_xip_max) << "\n";;

            // Loop over the xip_chiphmin values
            for(int i = 0 ; i < chipa_xip_dim ; i++)
            {
                file <<  xip_chiphmin_table[i] << "\n";
            }

            // Loop over the xip values
            for(int ichipa = 0 ; ichipa < chipa_xip_dim ; ichipa++)
            {
                for(int ichiph = 0 ; ichiph < chiph_xip_dim ; ichiph++)
                {
                    file <<  xip_table[ichipa*chiph_xip_dim+ichiph] << " ";
                }
                file << "\n";
            }

            file.close();
        }
    }
    else if (format == "binary")
    {
        std::ofstream file;
        file.open("tab_xip.bin",std::ios::binary);

        if (file.is_open()) {

            double temp0, temp1;

            temp0 = log10(chipa_xip_min);
            temp1 = log10(chipa_xip_max);

            file.write((char*)&chipa_xip_dim,sizeof (int));
            file.write((char*)&chiph_xip_dim,sizeof (int));
            file.write((char*)&temp0, sizeof (double));
            file.write((char*)&temp1, sizeof (double));

            // Write all table values of xip_chimin_table
            file.write((char*)&xip_chiphmin_table[0], sizeof (double)*chipa_xip_dim);

            // Write all table values of xip_table
            file.write((char*)&xip_table[0], sizeof (double)*chipa_xip_dim*chiph_xip_dim);

            file.close();
        }
    }
    else
    {
        MESSAGE("The output format " << format << " is not recognized");
    }
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
    return pow(1. + 4.8*(1.+chipa)*log10(1. + 1.7*chipa) + 2.44*chipa*chipa,-2./3.);
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

    // Parameters
    int rank; // Rank number 
    // timers
    double t0,t1;

    // Get the MPI rank
    rank = smpi->getRank();

    MESSAGE("        --- Integration F/chipa table:");

    t0 = MPI_Wtime();

    // Test if an external table exists, we read the table...
    if (Tools::file_exists("tab_integfochi.bin"))
    {

        MESSAGE("        Reading of the table:");
        
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

                // Resize of the array Integfochi before reading
                Integfochi.resize(dim_integfochi);

                // Read the table values
                file.read((char*)&Integfochi[0], sizeof (double)*dim_integfochi);

                file.close();
            }

        }

        // Resize of the array Integfochi before communication
        if (rank != 0)
        {
            Integfochi.resize(dim_integfochi);
        }

        // Allgathering of the table
        MPI_Bcast(&Integfochi[0], dim_integfochi, MPI_DOUBLE, 0, 
           smpi->getGlobalComm());


    }
    // else the table is generated
    else
    {
        // Temporary particle chi value
        double chipa; 
        // For percentages
        double pct = 0.;
        double dpct = 10.;
        // table load repartition
        int * imin_table;
        int * length_table;
        // Local array
        double * buffer; 
        int nb_ranks; // Number of ranks
        //int err;  // error MPI

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


        MESSAGE("        MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "        Rank: " << i 
                                          << " imin: " << imin_table[i] 
                                          << " length: " << length_table[i] );
            }
        }

        MESSAGE("        Computation:");
        dpct = std::max(dpct,100./length_table[rank]);
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
                MESSAGE("       " << i + 1<< "/" << length_table[rank] << " - " << (int)(std::round(pct)) << "%");
            }
        }

        // Communication of the data
        MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
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

    t1 = MPI_Wtime();
    MESSAGE("        done in " << (t1 - t0) << "s");
}


// ---------------------------------------------------------------------------------------------------------------------
//! Ouput in a file of the table values of integfochi
//
//! \param format output format (ascii, binary)
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

            file << "Dimension chipa - LOG10(chipa min) - LOG10(chipa max) \n";

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

            double temp0, temp1;

            temp0 = log10(chipa_integfochi_min);
            temp1 = log10(chipa_integfochi_max);

            file.write((char*)&dim_integfochi,sizeof (dim_integfochi));
            file.write((char*)&temp0, sizeof (double));
            file.write((char*)&temp1, sizeof (double));

            // Loop over the table values
            for(int i = 0 ; i < dim_integfochi ; i++)
            {
                file.write((char*)&Integfochi[i], sizeof (double));
            }

            file.close();
        }
    }
    else
    {
        MESSAGE("The output format " << format << " is not recognized");
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
