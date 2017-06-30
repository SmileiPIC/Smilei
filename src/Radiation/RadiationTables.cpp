// ----------------------------------------------------------------------------
//! \file RadiationTables.cpp
//
//! \brief This class contains the tables and the functions to generate them
//! for the Nonlinear Inverse Compton Scattering
//
//! \details The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#include "RadiationTables.h"

#include <cstring>
#include <iostream>
#include <fstream>

#include <cmath>

#include "userFunctions.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for RadiationTables
// -----------------------------------------------------------------------------
RadiationTables::RadiationTables()
{
    Integfochi.resize(0);
    xip_chiphmin_table.resize(0);
    xip_table.resize(0);

    integfochi_computed = false;
    xip_computed = false;

}

// -----------------------------------------------------------------------------
// Destructor for RadiationTables
// -----------------------------------------------------------------------------
RadiationTables::~RadiationTables()
{
}

// -----------------------------------------------------------------------------
//
//! Initialization of the parameters for the nonlinear inverse Compton
//! scattering
//
//! \param params Object Params for the parameters from the input script
// -----------------------------------------------------------------------------
void RadiationTables::initParams(Params& params)
{

    if (params.hasMCRadiation || params.hasContinuousRadiation)
    {
        TITLE("Initializing Radiation loss")

        // Preliminary checks
        if (params.referenceAngularFrequency_SI <= 0.)
            ERROR("The parameter `referenceAngularFrequency_SI` needs "
                  << "to be defined and positive to compute radiation losses");

    }

    if (params.hasMCRadiation)
    {
        MESSAGE("        The Monte-Carlo Compton radiation module is requested by some species.\n");
    }

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    // We read the properties
    if( PyTools::nComponents("RadiationLoss") != 0 )
    {

        // If Monte-Carlo radiation loss is requested
        if (params.hasMCRadiation)
        {

            // Extraction of the parameter from the input file
            PyTools::extract("chipa_h_min", chipa_h_min, "RadiationLoss");
            PyTools::extract("chipa_h_max", chipa_h_max, "RadiationLoss");
            PyTools::extract("h_dim", chipa_h_dim, "RadiationLoss");
            PyTools::extract("chipa_integfochi_min", chipa_integfochi_min, "RadiationLoss");
            PyTools::extract("chipa_integfochi_max", chipa_integfochi_max, "RadiationLoss");
            PyTools::extract("integfochi_dim", dim_integfochi, "RadiationLoss");
            PyTools::extract("chipa_xip_min", chipa_xip_min, "RadiationLoss");
            PyTools::extract("chipa_xip_max", chipa_xip_max, "RadiationLoss");
            PyTools::extract("xip_power", xip_power, "RadiationLoss");
            PyTools::extract("xip_threshold", xip_threshold, "RadiationLoss");
            PyTools::extract("chipa_xip_dim", chipa_xip_dim, "RadiationLoss");
            PyTools::extract("chiph_xip_dim", chiph_xip_dim, "RadiationLoss");
            PyTools::extract("output_format", output_format, "RadiationLoss");

            // Path to the databases
            PyTools::extract("table_path", table_path, "RadiationLoss");

            // Discontinuous minimum threshold
            PyTools::extract("chipa_disc_min_threshold",
                             chipa_disc_min_threshold,"RadiationLoss");

            // Additional regularly used parameters
            log10_chipa_xip_min = log10(chipa_xip_min);
            log10_chipa_integfochi_min = log10(chipa_integfochi_min);
            inv_chiph_xip_dim_minus_one = 1./(chiph_xip_dim - 1.);
        }
    }

    // Computation of some parameters
    if (params.hasMCRadiation || params.hasContinuousRadiation)
    {

        // Computation of the normalized Compton wavelength
        norm_lambda_compton = red_planck_cst*params.referenceAngularFrequency_SI
                            / (electron_mass*c_vacuum*c_vacuum);

        // Computation of the factor factor_dNphdt
        factor_dNphdt = sqrt(3.)*fine_struct_cst/(2.*M_PI*norm_lambda_compton);

        // Computation of the factor for the classical radiated power
        factor_cla_rad_power = 2.*fine_struct_cst/(3.*norm_lambda_compton);

        MESSAGE( "        factor_cla_rad_power: " << factor_cla_rad_power)

    }

    if (params.hasMCRadiation)
    {
        MESSAGE( "        chipa_disc_min_threshold: " << chipa_disc_min_threshold);
        MESSAGE( "        table path: " << table_path);
    }

    MESSAGE("")

    // Some additional checks
    if (params.hasMCRadiation)
    {
        if (chipa_integfochi_min >= chipa_integfochi_max)
        {
            ERROR("chipa_integfochi_min (" << chipa_integfochi_min
                    << ") >= chipa_integfochi_max (" << chipa_integfochi_max << ")")
        }
        if (chipa_xip_min >= chipa_xip_max)
        {
            ERROR("chipa_xip_min (" << chipa_xip_min
                    << ") >= chipa_xip_max (" << chipa_xip_max << ")")
        }
    }
}

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the photon quantum parameter chiph for emission
//! ramdomly and using the tables xip and chiphmin
//
//! \param chipa particle quantum parameter
// ---------------------------------------------------------------------------------------------------------------------
double RadiationTables::compute_chiph_emission(double chipa)
{
    // Log10 of chipa
    double logchipa;
    double chiph;
    double chiph_xip_delta;
    // Random xip
    double xip;
    int ichipa;
    int ichiph;
    // For the interpolation
    double log10_chiphm;
    double log10_chiphp;
    double d;
    int ixip;

    logchipa = log10(chipa);

    // -------------------------------
    // index of chipa in xip_table
    // -------------------------------
    // Use floor so that chipa corresponding to ichipa is <= given chipa
    ichipa = int(floor((logchipa-log10_chipa_xip_min)*(inv_chipa_xip_delta)));

    // Checking that ichipa is in the range of the tables
    // Else we use the values at the boundaries
    if (ichipa < 0)
    {
        ichipa = 0;
    }
    else if (ichipa > chipa_xip_dim-1)
    {
        ichipa = chipa_xip_dim-1;
    }

    // ---------------------------------------
    // Search of the index ichiph for chiph
    // ---------------------------------------

    // First, we compute a random xip in [0,1]
    xip = Rand::uniform();

    // If the randomly computed xip if below the first one of the row,
    // we take the first one which corresponds to the minimal photon chiph
    if (xip <= xip_table[ichipa*chiph_xip_dim])
    {
        ichiph = 0;
        xip = xip_table[ichipa*chiph_xip_dim];
    }
    // Above the last xip of the row, the last one corresponds
    // to the maximal photon chiph
    else if (xip > xip_table[(ichipa+1)*chiph_xip_dim-2])
    {
        ichiph = chiph_xip_dim-2;
        xip = xip_table[(ichipa+1)*chiph_xip_dim-1];
        // If nearest point: ichiph = chiph_xip_dim-1
    }
    else
    {
        // Search for the corresponding index ichiph for xip
        ichiph = userFunctions::search_elem_in_array(
            &xip_table[ichipa*chiph_xip_dim],xip,chiph_xip_dim);
    }

    // Corresponding chipa for ichipa
    logchipa = ichipa*chipa_xip_delta+log10_chipa_xip_min;

    // Delta for the corresponding chipa
    chiph_xip_delta = (logchipa - xip_chiphmin_table[ichipa])
                    *inv_chiph_xip_dim_minus_one;

    // --------------------------------------------------------------------
    // Compute chiph
    // This method is slow but more accurate than taking the nearest point
    // --------------------------------------------------------------------

    // Computation of the final chiph by interpolation
    log10_chiphm = ichiph*chiph_xip_delta
           + xip_chiphmin_table[ichipa];
    log10_chiphp = log10_chiphm + chiph_xip_delta;

    ixip = ichipa*chiph_xip_dim + ichiph;

    d = (xip - xip_table[ixip]) / (xip_table[ixip+1] - xip_table[ixip]);

    // Chiph after linear interpolation in the logarithmic scale
    chiph = pow(10.,log10_chiphm*(1.0-d) + log10_chiphp*(d));

    // ------------------------------------------------------------
    // Compute chiph
    // Fastest method using the nearest point but less accurate
    // ------------------------------------------------------------

    /*chiph = pow(10.,ichiph*chiph_xip_delta
           + xip_chiphmin_table[ichipa]);*/

    // Debugging
    /*std::cerr << "ichiph: " << ichiph << " "
              << "ichipa: " << ichipa << " "
              << "" << xip_table[ichipa*chiph_xip_dim + ichiph] << " < "
              << "xip: " << xip << " "
              << " < " << xip_table[ichipa*chiph_xip_dim + ichiph+1] << " "
              << "logchipa: " << logchipa << " "
              << "" << pow(10,(ichipa)*chipa_xip_delta + log10(chipa_xip_min)) << " < "
              << "chipa: " << chipa << " "
              << pow(10,(ichipa+1)*chipa_xip_delta + log10(chipa_xip_min)) << " "
              << "chiph: " << chiph << " "
              << std::endl;*/

    return chiph;
}

// -----------------------------------------------------------------------------
// TABLE COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the table h that is a discetization of the h function
//! in the stochastic model of Niel.
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::compute_h_table(SmileiMPI *smpi)
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0,t1;
    // tabel_exists
    bool table_exists;

    // Get the MPI rank
    rank = smpi->getRank();

    MESSAGE("        --- h table:");

    // Initial timer
    t0 = MPI_Wtime();

    // If external tables are available, we read them
    table_exists = RadiationTables::read_h_table(smpi);

    // Else we compute them
    if (!table_exists)
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
        h_table.resize(dim_integfochi);

        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];

        // Computation of the delta
        chipa_h_delta = (log10(chipa_h_max)
                - log10_chipa_h_min)/(chipa_h_dim-1);

        // Inverse delta
        chipa_h_inv_delta = 1./chipa_h_delta;

        // Load repartition
        userFunctions::distribute_load_1d_table(nb_ranks,
                chipa_h_dim,
                imin_table,
                length_table);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]];


        MESSAGE("            MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "            Rank: " << i
                                          << " imin: " << imin_table[i]
                                          << " length: " << length_table[i] );
            }
        }

        MESSAGE("            Computation:");
        dpct = std::max(dpct,100./length_table[rank]);
        // Loop over the table values
        for(int i = 0 ; i < length_table[rank] ; i++)
        {
            chipa = pow(10.,(imin_table[rank] + i)*chipa_h_delta
                  + log10_chipa_h_min);

            // 100 iterations is in theory sufficient to get the convergence
            // at 1e-15
            buffer[i] = RadiationTables::get_h_Niel(chipa,400,1e-15);

            //std::cout << rank << " " << buffer[i] << std::endl;

            if (100.*i >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE("            " << i + 1<< "/" << length_table[rank]
                                       << " - " << (int)(std::round(pct))
                                       << "%");
            }
        }

        // Communication of the data
        MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
                &h_table[0], &length_table[0], &imin_table[0],
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

        // flag computed at true
        h_computed = true;

        // Free memory
        // delete buffer;
        // delete imin_table;
        delete length_table;

    }

    // Final timer
    t1 = MPI_Wtime();
    MESSAGE("        done in " << (t1 - t0) << "s");
}

// -----------------------------------------------------------------------------
//! Computation of the table values of integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::compute_integfochi_table(SmileiMPI *smpi)
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0,t1;
    // tabel_exists
    bool table_exists;

    // Get the MPI rank
    rank = smpi->getRank();

    MESSAGE("        --- Integration F/chipa table:");

    // Initial timer
    t0 = MPI_Wtime();

    // If external tables are available, we read them
    table_exists = RadiationTables::read_integfochi_table(smpi);

    // Else we compute them
    if (!table_exists)
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
                - log10_chipa_integfochi_min)/(dim_integfochi-1);

        // Inverse delta
        inv_delta_chipa_integfochi = 1./delta_chipa_integfochi;

        // Load repartition
        userFunctions::distribute_load_1d_table(nb_ranks,
                dim_integfochi,
                imin_table,
                length_table);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]];


        MESSAGE("            MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "            Rank: " << i
                                          << " imin: " << imin_table[i]
                                          << " length: " << length_table[i] );
            }
        }

        MESSAGE("            Computation:");
        dpct = std::max(dpct,100./length_table[rank]);
        // Loop over the table values
        for(int i = 0 ; i < length_table[rank] ; i++)
        {
            chipa = pow(10.,(imin_table[rank] + i)*delta_chipa_integfochi
                  + log10_chipa_integfochi_min);

            buffer[i] = RadiationTables::compute_integfochi(chipa,
                    0.98e-40*chipa,0.98*chipa,400,1e-15);

            //std::cout << rank << " " << buffer[i] << std::endl;

            if (100.*i >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE("            " << i + 1<< "/" << length_table[rank] << " - " << (int)(std::round(pct)) << "%");
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

        // flag computed at true
        integfochi_computed = true;

        // Free memory
        // delete buffer;
        // delete imin_table;
        delete length_table;

    }

    // Final timer
    t1 = MPI_Wtime();
    MESSAGE("        done in " << (t1 - t0) << "s");
}

// -----------------------------------------------------------------------------
//! Computation of the minimum photon quantum parameter chiphmin
//! for the xip array and computation of the xip array.
//
//! \details Under the minimum chiph value, photon energy is
//! considered negligible.
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// ---------------------------------------------------------------------------------------------------------------------
void RadiationTables::compute_xip_table(SmileiMPI *smpi)
{

    // Parameters
    int rank; // Rank number
    // timers
    double t0,t1;
    // Flag table exists
    bool table_exists;

    // Get the MPI rank
    rank = smpi->getRank();

    t0 = MPI_Wtime();

    MESSAGE("        --- Table chiphmin and xip:");

    // If external tables are available, we read them
    table_exists = RadiationTables::read_xip_table(smpi);

    // Else we compute them
    if (!table_exists)
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
                - log10_chipa_xip_min)/(chipa_xip_dim-1);

        // Inverse of delta
        inv_chipa_xip_delta = 1./chipa_xip_delta;

        // Load repartition
        userFunctions::distribute_load_1d_table(nb_ranks,
                chipa_xip_dim,
                imin_table,
                length_table);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]];

        MESSAGE("            MPI repartition:");
        // Print repartition
        if (rank==0)
        {
            for(int i =0 ; i < nb_ranks ; i++)
            {
                MESSAGE( "            Rank: " << i
                                          << " imin: "   << imin_table[i]
                                          << " length: " << length_table[i] );
            }
        }

        // 1. - Computation of xip_chiphmin_table
        MESSAGE("            Computation of log10(chiphmin):");
        dpct = std::max(dpct,100./length_table[rank]);

        // Loop for chiphmin
        for(int ichipa = 0 ; ichipa < length_table[rank] ; ichipa++)
        {

            xip = 1;
            logchiphmin = (imin_table[rank] + ichipa)*chipa_xip_delta
                  + log10_chipa_xip_min;
            chipa = pow(10.,logchiphmin);

            // Denominator of xip
            denominator = RadiationTables::compute_integfochi(chipa,
                    0.98e-40*chipa,0.98*chipa,200,1e-13);

            k = 0;
            while(k < xip_power)
            {
                logchiphmin -= pow(0.1,k);
                chiph = pow(10.,logchiphmin);
                numerator = RadiationTables::compute_integfochi(chipa,
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
            buffer[ichipa] = logchiphmin;

            // display percentage
            if (100.*ichipa >= length_table[rank]*pct)
            {
                pct += dpct;
                MESSAGE( "            " << ichipa + 1 << "/" << length_table[rank]
                                    << " - " << (int)(std::round(pct)) << "%");
            }
        }

        // Communication of the xip_chiphmin table
        MPI_Allgatherv(&buffer[0], length_table[rank], MPI_DOUBLE,
                &xip_chiphmin_table[0], &length_table[0], &imin_table[0],
                MPI_DOUBLE, smpi->getGlobalComm());

        // 2. - Computation of the xip table
        MESSAGE("            Computation of xip:");
        dpct = std::max(dpct,100./length_table[rank]);

        // Allocation of the local buffer
        buffer = new double [length_table[rank]*chiph_xip_dim];

        // Loop for xip in the chipa dimension
        pct = 0;
        dpct = std::max(dpct,100./(length_table[rank]*chiph_xip_dim));
        for(int ichipa = 0 ; ichipa < length_table[rank] ; ichipa++)
        {

            chipa = pow(10.,(imin_table[rank] + ichipa)*chipa_xip_delta + log10_chipa_xip_min);

            chiph_delta = (log10(chipa) - xip_chiphmin_table[imin_table[rank] + ichipa])
                        / (chiph_xip_dim - 1);

            // Denominator of xip
            denominator = RadiationTables::compute_integfochi(chipa,
                    1e-40*chipa,chipa,250,1e-15);

            // Loop in the chiph dimension
            for (int ichiph = 0 ; ichiph < chiph_xip_dim ; ichiph ++)
            {
                // Local chiph value
               chiph = pow(10.,ichiph*chiph_delta +
                   xip_chiphmin_table[imin_table[rank] + ichipa]);

               /* std::cout << "rank: " << rank
                         << " " << chipa
                         << " " << chiph
                         << " " << xip_chiphmin_table[imin_table[rank] + ichipa] << std::endl; */

               // Denominator of xip
               numerator = RadiationTables::compute_integfochi(chipa,
                       1e-40*chiph,chiph,250,1e-15);

               // Update local buffer value
               buffer[ichipa*chiph_xip_dim + ichiph] = std::min(1.,numerator / denominator);

               // If == 1, end of the loop
               if (buffer[ichipa*chiph_xip_dim + ichiph] == 1)
               {
                   for (int i = ichiph+1 ; i < chiph_xip_dim ; i ++)
                   {
                       buffer[ichipa*chiph_xip_dim + i] = 1.;
                   }
                   ichiph = chiph_xip_dim;
               }

               // display percentage
               if (100.*(ichipa*chiph_xip_dim+ichiph)
                   >= length_table[rank]*chiph_xip_dim*pct)
               {
                   pct += dpct;
                   MESSAGE( "            " << ichipa*chiph_xip_dim+ichiph + 1
                                           << "/"
                                           << length_table[rank]*chiph_xip_dim
                                           << " - " << (int)(std::round(pct))
                                           << "%");
               }

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

        // flag computed at true
        xip_computed = true;

        // clean temporary arrays
        delete buffer;
        delete length_table;
        delete imin_table;

    }

    t1 = MPI_Wtime();
    MESSAGE("        done in " << (t1 - t0) << "s");

}

// ---------------------------------------------------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//
//! \param params list of simulation parameters
//! \param smpi MPI parameters
// ---------------------------------------------------------------------------------------------------------------------
void RadiationTables::compute_tables(Params& params, SmileiMPI *smpi)
{
    // These tables are loaded only if if one species has Monte-Carlo Compton radiation
    if (params.hasMCRadiation)
    {
        RadiationTables::compute_integfochi_table(smpi);
        RadiationTables::compute_xip_table(smpi);
    }
}

// -----------------------------------------------------------------------------
// TABLE OUTPUTS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Ouput in a file of the table values of integfochi
//
// -----------------------------------------------------------------------------
void RadiationTables::output_integfochi_table()
{

    if (output_format == "ascii")
    {
        std::ofstream file;
        file.open(table_path + "/tab_integfochi.dat");

        if (file.is_open()) {

            file.precision(12);

            file << "Table Integration F(CHI)/CHI for Nonlinear Compton Scattering \n";

            file << "Dimension chipa - chipa min - chipa max \n";

            file << dim_integfochi ;
            file << " "
                << chipa_integfochi_min << " "
                << chipa_integfochi_max << "\n";;

            // Loop over the table values
            for(int i = 0 ; i < dim_integfochi ; i++)
            {
                file <<  Integfochi[i] << "\n";
            }

            file.close();
        }
    }
    else if (output_format == "binary")
    {
        std::ofstream file;
        file.open(table_path + "/tab_integfochi.bin",std::ios::binary);

        if (file.is_open()) {

            double temp0, temp1;

            temp0 = chipa_integfochi_min;
            temp1 = chipa_integfochi_max;

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
    // HDF5
    // The table is written as a dataset
    else if (output_format == "hdf5")
    {

        hid_t       fileId;
        hid_t       datasetId;
        hid_t       dataspaceId;
        hsize_t     dims;
        std::string buffer;

        buffer = table_path + "/radiation_tables.h5";

        // We first check whether the file already exists
        // If yes, we simply open the file
        if (Tools::file_exists(buffer))
        {
            fileId = H5Fopen(buffer.c_str(),
                              H5F_ACC_RDWR,
                              H5P_DEFAULT);

        }
        // Else, we create the file
        else
        {
            fileId  = H5Fcreate(buffer.c_str(),
                                 H5F_ACC_TRUNC,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
        }

        // Create the data space for the dataset.
        dims = dim_integfochi;
        dataspaceId = H5Screate_simple(1, &dims, NULL);

        // Creation of the dataset
        datasetId = H5Dcreate(fileId,
                              "integfochi",
                              H5T_NATIVE_DOUBLE,
                              dataspaceId,
                              H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

        // Fill the dataset
        H5Dwrite(datasetId, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 &Integfochi[0]);

        // Create attributes
        H5::attr(datasetId, "chipa_min", chipa_integfochi_min);
        H5::attr(datasetId, "chipa_max", chipa_integfochi_max);
        H5::attr(datasetId, "chipa_dim", dim_integfochi);

        // Close everything
        H5Dclose(datasetId);
        H5Sclose(dataspaceId);
        H5Fclose(fileId);

    }
    else
    {
        MESSAGE("The table output format " << output_format
             << " is not recognized");
    }
}

// -----------------------------------------------------------------------------
//! File output of xip_chiphmin_table and xip_table
//
// -----------------------------------------------------------------------------
void RadiationTables::output_xip_table()
{

    if (output_format == "ascii")
    {
        std::ofstream file;
        file.open(table_path + "/tab_xip.dat");

        if (file.is_open()) {

            file.precision(12);

            file << "Table xip_chiphmin and xip for Nonlinear Compton Scattering \n";

            file << "Dimension chipa - Dimension chiph - chipa min - chipa max \n";

            file << chipa_xip_dim << " "
                 << chiph_xip_dim << " "
                 << chipa_xip_min << " "
                 << chipa_xip_max << "\n";;

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
    else if (output_format == "binary")
    {
        std::ofstream file;
        file.open(table_path + "/tab_xip.bin",std::ios::binary);

        if (file.is_open()) {

            double temp0, temp1;

            temp0 = chipa_xip_min;
            temp1 = chipa_xip_max;

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
    // HDF5
    // The table is written as a dataset
    else if (output_format == "hdf5")
    {

        hid_t       fileId;
        hid_t       datasetId;
        hid_t       dataspaceId;
        hsize_t     dims[2];
        std::string buffer;

        buffer = table_path + "/radiation_tables.h5";

        // We first check whether the file already exists
        // If yes, we simply open the file
        if (Tools::file_exists(buffer))
        {
            fileId = H5Fopen(buffer.c_str(),
                              H5F_ACC_RDWR,
                              H5P_DEFAULT);

        }
        // Else, we create the file
        else
        {
            fileId  = H5Fcreate(buffer.c_str(),
                                 H5F_ACC_TRUNC,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT);
        }

        // Creation of the datasat chiphmin
        dims[0] = chipa_xip_dim;
        dataspaceId = H5Screate_simple(1, dims, NULL);
        datasetId = H5Dcreate(fileId,
                              "xip_chiphmin",
                              H5T_NATIVE_DOUBLE,
                              dataspaceId,
                              H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

        // Fill the dataset chiphmin
        H5Dwrite(datasetId, H5T_NATIVE_DOUBLE,
               H5S_ALL, H5S_ALL, H5P_DEFAULT,
               &xip_chiphmin_table[0]);

        // Attribute creation
        H5::attr(datasetId, "chipa_min", chipa_xip_min);
        H5::attr(datasetId, "chipa_max", chipa_xip_max);
        H5::attr(datasetId, "chipa_dim", chipa_xip_dim);
        H5::attr(datasetId, "power", xip_power);
        H5::attr(datasetId, "threshold", xip_threshold);

        // Creation of the datasat chiphmin xip
        dims[0] = chipa_xip_dim;
        dims[1] = chiph_xip_dim;
        dataspaceId = H5Screate_simple(2, dims, NULL);
        datasetId = H5Dcreate(fileId,
                              "xip",
                              H5T_NATIVE_DOUBLE,
                              dataspaceId,
                              H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

        // Fill the dataset
        H5Dwrite(datasetId, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 &xip_table[0]);

        // Attribute creation
        H5::attr(datasetId, "chipa_min", chipa_xip_min);
        H5::attr(datasetId, "chipa_max", chipa_xip_max);
        H5::attr(datasetId, "chipa_dim", chipa_xip_dim);
        H5::attr(datasetId, "chiph_dim", chiph_xip_dim);
        H5::attr(datasetId, "power", xip_power);
        H5::attr(datasetId, "threshold", xip_threshold);

        // Close everything
        H5Dclose(datasetId);
        H5Sclose(dataspaceId);
        H5Fclose(fileId);
    }
    else
    {
        MESSAGE("The output format " << output_format << " is not recognized");
    }
}

// -----------------------------------------------------------------------------
//! Output the computed tables so that thay can be read at the next run.
//! Table output by the master MPI rank
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::output_tables(SmileiMPI *smpi)
{
    // Sequential output
    if (smpi->isMaster())
    {
        // If tables have been computed, they are output on the disk
        // to be used for the next run
        if (integfochi_computed)
        {
            RadiationTables::output_integfochi_table();
        }
        if (xip_computed)
        {
            RadiationTables::output_xip_table();
        }
    }
}

// -----------------------------------------------------------------------------
// PHYSICAL COMPUTATION
// -----------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------------------------------------
//! Computation of the Cross Section dNph/dt which is also
//! the number of photons generated per time unit.
//
//! \param chipa particle quantum parameter
//! \param gfpa particle gamma factor
// ---------------------------------------------------------------------------------------------------------------------
double RadiationTables::compute_dNphdt(double chipa,double gfpa)
{

    // Log of the particle quantum parameter chipa
    double logchipa;
    double logchipam;
    double logchipap;
    // Index
    int ichipa;
    // final value
    double dNphdt;

    logchipa = log10(chipa);

    // Lower index for interpolation in the table integfochi
    ichipa = int(floor((logchipa-log10_chipa_integfochi_min)
                 *inv_delta_chipa_integfochi));

    // If we are not in the table...
    if (ichipa < 0)
    {
        ichipa = 0;
        dNphdt = Integfochi[ichipa];
    }
    else if (ichipa >= dim_integfochi-1)
    {
        ichipa = dim_integfochi-2;
        dNphdt = Integfochi[ichipa];
    }
    else
    {
        // Upper and lower values for linear interpolation
        logchipam = ichipa*delta_chipa_integfochi + log10_chipa_integfochi_min;
        logchipap = logchipam + delta_chipa_integfochi;

        // Interpolation
        dNphdt = (Integfochi[ichipa+1]*fabs(logchipa-logchipam) +
                Integfochi[ichipa]*fabs(logchipap - logchipa))*inv_delta_chipa_integfochi;
    }

    /*std::cerr << "factor_dNphdt: " << factor_dNphdt << " "
              << "dNphdt: " << dNphdt << " "
              << "chipa: " << chipa << " "
              << "ichipa: " << ichipa << " "
              << "" << logchipam << " < logchipa: " << logchipa << " < " << logchipap << " "
              << "delta_chipa_integfochi: " << delta_chipa_integfochi << " "
              << "Integfochi[ichipa]: " << Integfochi[ichipa] << " "
              << "Integfochi[ichipa+1]: " << Integfochi[ichipa+1] << " "
              << std::endl;*/

    return factor_dNphdt*dNphdt*chipa/gfpa;

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
double RadiationTables::compute_integfochi(double chipa,
        double chipmin,
        double chipmax,
        int nbit,
        double eps)
{

    //std::cout << "RadiationTables::compute_integfochi" << std::endl;

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
        sync_emi = RadiationTables::compute_sync_emissivity_ritus(chipa,chiph,200,1e-15);
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
double RadiationTables::compute_sync_emissivity_ritus(double chipa,
        double chiph, int nbit, double eps)
{

    //std::cout << "RadiationTables::compute_sync_emissivity_ritus" << std::endl;

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
        ERROR("In compute_sync_emissivity_ritus: chipa " << chipa
                                                         << " < chiph "
                                                         << chiph);
    }
}

// -----------------------------------------------------------------------------
//! Return the value of the function h(chipa) of Niel et al.
//! Use an integration of Gauss-Legendre
//
//! \param chipa particle quantum parameter
//! \param nbit number of iterations for the Gauss-Legendre integration
//! \param eps epsilon for the modified bessel function
// -----------------------------------------------------------------------------
double RadiationTables::get_h_Niel(double chipa,
        int nbit, double eps)
{
    // Arrays for Gauss-Legendre integration
    double * gauleg_x = new double[nbit];
    double * gauleg_w = new double[nbit];

    return 0;

}

// -----------------------------------------------------------------------------
// TABLE READING
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Read the external table h
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_h_table(SmileiMPI *smpi)
{
    // Flag database available
    bool table_exists = false;

    // Test if an external table exists, if yes we read the table...
    // Binary table
    if (Tools::file_exists(table_path + "/tab_h.bin"))
    {

        table_exists = true;

        if (smpi->getRank()==0)
        {

            // Reading of the table file
            std::ifstream file;
            file.open(table_path + "/tab_h.bin",std::ios::binary);

            if (file.is_open())
            {

                // Read the header
                file.read((char*)&chipa_h_dim,
                          sizeof (chipa_h_dim));
                file.read((char*)&chipa_h_min,
                          sizeof (chipa_h_min));
                file.read((char*)&chipa_h_max,
                          sizeof (chipa_h_max));

                // Resize of the array Integfochi before reading
                h_table.resize(chipa_h_dim);

                // Read the table values
                file.read((char*)&h_table[0], sizeof (double)*chipa_h_dim);

                file.close();
            }

        }

    }
    // HDF5 format
    else if (Tools::file_exists(table_path + "/radiation_tables.h5"))
    {
         hid_t       fileId;
         hid_t       datasetId;
         std::string buffer;

         if (smpi->getRank()==0)
         {

             buffer = table_path + "/radiation_tables.h5";

             fileId = H5Fopen(buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

             datasetId = H5Dopen2(fileId, "h", H5P_DEFAULT);

             // If this dataset exists, we read it
             if (datasetId > 0)
             {

                 table_exists = true;

                 // First, we read attributes
                 H5::getAttr(datasetId, "chipa_dim", chipa_h_dim);
                 H5::getAttr(datasetId, "chipa_min", chipa_h_min);
                 H5::getAttr(datasetId, "chipa_max", chipa_h_max);

                 // Resize of the array Integfochi before reading
                 h_table.resize(chipa_h_dim);

                 // then the dataset
                 H5Dread(datasetId,
                        H5T_NATIVE_DOUBLE, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT,
                        &h_table[0]);

                 H5Dclose(datasetId);
                 H5Fclose(fileId);
             }
             // Else, we will have to compute it
             else
             {
                 table_exists = false;
             }
        }

        // Bcast table_exists
        MPI_Bcast(&table_exists, 1, MPI_INTEGER, 0,smpi->getGlobalComm());

    }

    // If the table exists, they have been read...
    if (table_exists)
    {

        MESSAGE("            Reading of the external database");
        MESSAGE("            Dimension quantum parameter: " << dim_integfochi);
        MESSAGE("            Minimum particle quantum parameter chi: " << chipa_integfochi_min);
        MESSAGE("            Maximum particle quantum parameter chi: " << chipa_integfochi_max);

        // Bcast the table to all MPI ranks
        RadiationTables::bcast_integfochi_table(smpi);
    }

    return table_exists;

}

// -----------------------------------------------------------------------------
//! Read the external table integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_integfochi_table(SmileiMPI *smpi)
{
    // Flag database available
    bool table_exists = false;

    // Test if an external table exists, if yes we read the table...
    // Binary table
    if (Tools::file_exists(table_path + "/tab_integfochi.bin"))
    {

        table_exists = true;

        if (smpi->getRank()==0)
        {

            // Reading of the table file
            std::ifstream file;
            file.open(table_path + "/tab_integfochi.bin",std::ios::binary);

            if (file.is_open())
            {

                // Read the header
                file.read((char*)&dim_integfochi,sizeof (dim_integfochi));
                file.read((char*)&chipa_integfochi_min,
                          sizeof (chipa_integfochi_min));
                file.read((char*)&chipa_integfochi_max,
                          sizeof (chipa_integfochi_max));

                // Resize of the array Integfochi before reading
                Integfochi.resize(dim_integfochi);

                // Read the table values
                file.read((char*)&Integfochi[0], sizeof (double)*dim_integfochi);

                file.close();
            }

        }

    }
    // HDF5 format
    else if (Tools::file_exists(table_path + "/radiation_tables.h5"))
    {
         hid_t       fileId;
         hid_t       datasetId;
         std::string buffer;

         if (smpi->getRank()==0)
         {

             buffer = table_path + "/radiation_tables.h5";

             fileId = H5Fopen(buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

             datasetId = H5Dopen2(fileId, "integfochi", H5P_DEFAULT);

             // If this dataset exists, we read it
             if (datasetId > 0)
             {

                 table_exists = true;

                 // First, we read attributes
                 H5::getAttr(datasetId, "chipa_dim", dim_integfochi);
                 H5::getAttr(datasetId, "chipa_min", chipa_integfochi_min);
                 H5::getAttr(datasetId, "chipa_max", chipa_integfochi_max);

                 // Resize of the array Integfochi before reading
                 Integfochi.resize(dim_integfochi);

                 // then the dataset
                 H5Dread(datasetId,
                        H5T_NATIVE_DOUBLE, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT,
                        &Integfochi[0]);

                 H5Dclose(datasetId);
                 H5Fclose(fileId);
             }
             // Else, we will have to compute it
             else
             {
                 table_exists = false;
             }
        }

        // Bcast table_exists
        MPI_Bcast(&table_exists, 1, MPI_INTEGER, 0,smpi->getGlobalComm());

    }

    // If the table exists, they have been read...
    if (table_exists)
    {

        MESSAGE("            Reading of the external database");
        MESSAGE("            Dimension quantum parameter: " << dim_integfochi);
        MESSAGE("            Minimum particle quantum parameter chi: " << chipa_integfochi_min);
        MESSAGE("            Maximum particle quantum parameter chi: " << chipa_integfochi_max);

        // Bcast the table to all MPI ranks
        RadiationTables::bcast_integfochi_table(smpi);
    }

    return table_exists;

}

// -----------------------------------------------------------------------------
//! Read the external table xip_chiphmin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool RadiationTables::read_xip_table(SmileiMPI *smpi)
{
    // Flag database available
    bool table_exists = false;

    // Test if an external table exists, we read the table...
    if (Tools::file_exists(table_path + "/tab_xip.bin"))
    {

        table_exists = true;

        if (smpi->getRank()==0)
        {

            // Reading of the table file
            std::ifstream file;
            file.open(table_path + "/tab_xip.bin",std::ios::binary);

            if (file.is_open())
            {

                // Read the header
                file.read((char*)&chipa_xip_dim,sizeof (chipa_xip_dim));
                file.read((char*)&chiph_xip_dim,sizeof (chiph_xip_dim));
                file.read((char*)&chipa_xip_min, sizeof (chipa_xip_min));
                file.read((char*)&chipa_xip_max, sizeof (chipa_xip_max));

                // Allocation of the array xip
                xip_chiphmin_table.resize(chipa_xip_dim);
                xip_table.resize(chipa_xip_dim*chiph_xip_dim);

                // Read the table values
                file.read((char*)&xip_chiphmin_table[0], sizeof (double)*chipa_xip_dim);

                // Read the table values
                file.read((char*)&xip_table[0], sizeof (double)*chipa_xip_dim*chiph_xip_dim);

                file.close();
            }

        }
    }
    // HDF5 format
    else if (Tools::file_exists(table_path + "/radiation_tables.h5"))
    {
         if (smpi->getRank()==0)
         {

             hid_t       fileId;
             hid_t       datasetId_chiphmin;
             hid_t       datasetId_xip;
             std::string buffer;

             buffer = table_path + "/radiation_tables.h5";

             fileId = H5Fopen(buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

             datasetId_chiphmin = H5Dopen2(fileId, "xip_chiphmin", H5P_DEFAULT);
             datasetId_xip = H5Dopen2(fileId, "xip", H5P_DEFAULT);

             // If this dataset exists, we read it
             if (datasetId_chiphmin > 0 && datasetId_xip > 0)
             {

                 table_exists = true;

                 // First, we read attributes
                 H5::getAttr(datasetId_xip, "chipa_dim", chipa_xip_dim);
                 H5::getAttr(datasetId_xip, "chiph_dim", chiph_xip_dim);
                 H5::getAttr(datasetId_xip, "chipa_min", chipa_xip_min);
                 H5::getAttr(datasetId_xip, "chipa_max", chipa_xip_max);

                 // Allocation of the array xip
                 xip_chiphmin_table.resize(chipa_xip_dim);
                 xip_table.resize(chipa_xip_dim*chiph_xip_dim);

                 // then the dataset for chiphmin
                 H5Dread(datasetId_chiphmin,
                        H5T_NATIVE_DOUBLE, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT,
                        &xip_chiphmin_table[0]);

                 // then the dataset for xip
                 H5Dread(datasetId_xip,
                        H5T_NATIVE_DOUBLE, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT,
                        &xip_table[0]);

                H5Dclose(datasetId_xip);
                H5Dclose(datasetId_chiphmin);
                H5Fclose(fileId);
             }
             // Else, we will have to compute it
             else
             {
                 table_exists = false;
             }
        }

        // Bcast table_exists
        MPI_Bcast(&table_exists, 1, MPI_INTEGER, 0,smpi->getGlobalComm());

    }

    // If the table exists, they have been read...
    if (table_exists)
    {

        MESSAGE("            Reading of the external database");
        MESSAGE("            Dimension particle chi: " << chipa_xip_dim);
        MESSAGE("            Dimension photon chi: " << chiph_xip_dim);
        MESSAGE("            Minimum particle chi: " << chipa_xip_min);
        MESSAGE("            Maximum particle chi: " << chipa_xip_max);

        // Bcast the table to all MPI ranks
        RadiationTables::bcast_xip_table(smpi);
    }

    return table_exists;
}

// -----------------------------------------------------------------------------
// TABLE COMMUNICATIONS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Bcast of the external table integfochi
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::bcast_integfochi_table(SmileiMPI *smpi)
{
    // Position for MPI pack and unack
    int position;
    // buffer size for MPI pack and unpack
    int buf_size;

    // -------------------------------------------------------
    // Bcast of all the parameters
    // We pack everything in a buffer
    // --------------------------------------------------------

    // buffer size
    if (smpi->getRank() == 0)
    {
        MPI_Pack_size(1, MPI_INTEGER, smpi->getGlobalComm(), &position);
        buf_size = position;
        MPI_Pack_size(2, MPI_DOUBLE, smpi->getGlobalComm(), &position);
        buf_size += position;
        MPI_Pack_size(dim_integfochi, MPI_DOUBLE, smpi->getGlobalComm(),
                      &position);
        buf_size += position;
    }

    MESSAGE("            Buffer size: " << buf_size);

    // Exchange buf_size with all ranks
    MPI_Bcast(&buf_size, 1, MPI_INTEGER, 0,smpi->getGlobalComm());

    // Packet that will contain all parameters
    char * buffer = new char[buf_size];

    // Proc 0 packs
    if (smpi->getRank() == 0)
    {
        position = 0;
        MPI_Pack(&dim_integfochi,
             1,MPI_INTEGER,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&chipa_integfochi_min,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&chipa_integfochi_max,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

        MPI_Pack(&Integfochi[0],dim_integfochi,
            MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

    }

    // Bcast all parameters
    MPI_Bcast(&buffer[0], buf_size, MPI_PACKED, 0,smpi->getGlobalComm());

    // Other ranks unpack
    if (smpi->getRank() != 0)
    {
        position = 0;
        MPI_Unpack(buffer, buf_size, &position,
                   &dim_integfochi, 1, MPI_INTEGER,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &chipa_integfochi_min, 1, MPI_DOUBLE,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &chipa_integfochi_max, 1, MPI_DOUBLE,smpi->getGlobalComm());

        // Resize table before unpacking values
        Integfochi.resize(dim_integfochi);

        MPI_Unpack(buffer, buf_size, &position,&Integfochi[0],
                    dim_integfochi, MPI_DOUBLE,smpi->getGlobalComm());

    }

    log10_chipa_integfochi_min = log10(chipa_integfochi_min);

    // Computation of the delta
    delta_chipa_integfochi = (log10(chipa_integfochi_max)
            - log10_chipa_integfochi_min)/(dim_integfochi-1);

    // Inverse delta
    inv_delta_chipa_integfochi = 1./delta_chipa_integfochi;
}

// -----------------------------------------------------------------------------
//! Bcast of the external table xip_chiphmin and xip
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void RadiationTables::bcast_xip_table(SmileiMPI *smpi)
{
    // Position for MPI pack and unack
    int position = 0;
    // buffer size for MPI pack and unpack
    int buf_size = 0;

    // -------------------------------------------
    // Bcast of all the parameters
    // We pack everything in a buffer
    // -------------------------------------------

    // Compute the buffer size
    if (smpi->getRank() == 0)
    {
        MPI_Pack_size(2, MPI_INTEGER, smpi->getGlobalComm(), &position);
        buf_size = position;
        MPI_Pack_size(2, MPI_DOUBLE, smpi->getGlobalComm(), &position);
        buf_size += position;
        MPI_Pack_size(chipa_xip_dim, MPI_DOUBLE, smpi->getGlobalComm(),
                      &position);
        buf_size += position;
        MPI_Pack_size(chipa_xip_dim*chiph_xip_dim, MPI_DOUBLE,
                      smpi->getGlobalComm(), &position);
        buf_size += position;
    }

    MESSAGE("            Buffer size for MPI exchange: " << buf_size);

    // Exchange buf_size with all ranks
    MPI_Bcast(&buf_size, 1, MPI_INTEGER, 0,smpi->getGlobalComm());

    // Packet that will contain all parameters
    char * buffer = new char[buf_size];

    // Proc 0 packs
    if (smpi->getRank() == 0)
    {
        position = 0;
        MPI_Pack(&chipa_xip_dim,
             1,MPI_INTEGER,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&chiph_xip_dim,
             1,MPI_INTEGER,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&chipa_xip_min,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&chipa_xip_max,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

        MPI_Pack(&xip_chiphmin_table[0],chipa_xip_dim,
            MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

        MPI_Pack(&xip_table[0],chipa_xip_dim*chiph_xip_dim,
            MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());
    }

    // Bcast all parameters
    MPI_Bcast(&buffer[0], buf_size, MPI_PACKED, 0,smpi->getGlobalComm());

    // Other ranks unpack
    if (smpi->getRank() != 0)
    {
        position = 0;
        MPI_Unpack(buffer, buf_size, &position,
                   &chipa_xip_dim, 1, MPI_INTEGER,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &chiph_xip_dim, 1, MPI_INTEGER,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &chipa_xip_min, 1, MPI_DOUBLE,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &chipa_xip_max, 1, MPI_DOUBLE,smpi->getGlobalComm());

        // Resize tables before unpacking values
        xip_chiphmin_table.resize(chipa_xip_dim);
        xip_table.resize(chipa_xip_dim*chiph_xip_dim);

        MPI_Unpack(buffer, buf_size, &position,&xip_chiphmin_table[0],
                    chipa_xip_dim, MPI_DOUBLE,smpi->getGlobalComm());

        MPI_Unpack(buffer, buf_size, &position,&xip_table[0],
            chipa_xip_dim*chiph_xip_dim, MPI_DOUBLE,smpi->getGlobalComm());
    }

    // Log10 of chipa_xip_min for efficiency
    log10_chipa_xip_min = log10(chipa_xip_min);

    // Computation of the delta
    chipa_xip_delta = (log10(chipa_xip_max)
           - log10_chipa_xip_min)/(chipa_xip_dim-1);

    // Inverse of delta
    inv_chipa_xip_delta = 1./chipa_xip_delta;

    // Inverse chiph discetization (regularly used)
    inv_chiph_xip_dim_minus_one = 1./(chiph_xip_dim - 1.);

}
