// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheelerTables.cpp
//
//! \brief This file is for the methods of the class
//! MultiphotonBreitWheelerTables.
//! This class contains the methods and tools to generate and manage
//! the physical tables for the multiphoton Breit-wheeler process.
//
//! \details
//! The implementation is adapted from the thesis results of M. Lobet
//! See http://www.theses.fr/2015BORD0361
// ----------------------------------------------------------------------------

#include "MultiphotonBreitWheelerTables.h"

// -----------------------------------------------------------------------------
// INITILIZATION AND DESTRUCTION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// Constructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::MultiphotonBreitWheelerTables()
{
    // T table
    T_table.resize(0);
    T_computed = false;
}

// -----------------------------------------------------------------------------
// Destructor for MutliphotonBreitWheelerTables
// -----------------------------------------------------------------------------
MultiphotonBreitWheelerTables::~MultiphotonBreitWheelerTables()
{
}

// -----------------------------------------------------------------------------
//
//! Initialization of the parameters for the multiphoton Breit-Wheeler pair
//! creation
//
//! \param params Object Params for the parameters from the input script
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::initialization(Params& params)
{
    if (params.hasMultiphotonBreitWheeler)
    {
        TITLE("Initializing mutliphoton Breit-Wheeler")

        // Preliminary checks
        if (params.referenceAngularFrequency_SI <= 0.)
            ERROR("The parameter `referenceAngularFrequency_SI` needs "
                  << "to be defined and positive to compute radiation losses");

    }

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    // We read the properties
    if(PyTools::nComponents("MultiphotonBreitWheeler"))
    {

        // Extraction of the parameter from the input file
        PyTools::extract("T_chiph_min", T_chiph_min, "MultiphotonBreitWheeler");
        PyTools::extract("T_chiph_max", T_chiph_max, "MultiphotonBreitWheeler");
        PyTools::extract("T_dim", T_dim, "MultiphotonBreitWheeler");

        T_log10_chiph_min = log10(T_chiph_min);

        // Format of the tables
        PyTools::extract("output_format", output_format, "MultiphotonBreitWheeler");

        // Path to the databases
        PyTools::extract("table_path", table_path, "MultiphotonBreitWheeler");
    }

    // Computation of some parameters
    if (params.hasMultiphotonBreitWheeler)
    {
        // Computation of the normalized Compton wavelength
        norm_lambda_compton = params.red_planck_cst*params.referenceAngularFrequency_SI
                            / (params.electron_mass*params.c_vacuum*params.c_vacuum);

        // Computation of the factor factor_dNBWdt
        factor_dNBWdt = params.fine_struct_cst/(norm_lambda_compton);
    }

    // Messages and checks
    if (params.hasMultiphotonBreitWheeler)
    {
        MESSAGE( "        table path: " << table_path);

        if (T_chiph_min >= T_chiph_max)
        {
            ERROR("T_chiph_min (" << T_chiph_min
                    << ") >= T_chiph_max (" << T_chiph_max << ")")
        }
    }

    MESSAGE("")
}

// -----------------------------------------------------------------------------
// PHYSICAL COMPUTATION
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Computation of the production rate of pairs per photon
//! \param chiph photon quantum parameter
//! \param gamma photon normalized energy
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_dNBWdt(double chiph, double gamma)
{
    // ________________________________________
    // Parameters

    // Log of the photon quantum parameter chipa
    double logchiph;
    double logchiphm;
    double logchiphp;
    // Index
    int ichiph;
    // final value
    double dNBWdt;

    // ________________________________________
    // Computation

    logchiph = log10(chiph);

    // Lower index for interpolation in the table integfochi
    ichiph = int(floor((logchiph-T_log10_chiph_min)
                 *T_chiph_inv_delta));

     // If chiph is below the lower bound of the table
     // An asymptotic approximation is used
     if (ichiph < 0)
     {
         ichiph = 0;
         dNBWdt = 0.46*exp(-8./(3.*chiph));
     }
     // If chiph is above the upper bound of the table
     // An asymptotic approximation is used
     else if (ichiph >= T_dim-1)
     {
        ichiph = T_dim-2;
        dNBWdt = 0.38*pow(chiph,-1./3.);
     }
     else
     {
         // Upper and lower values for linear interpolation
         logchiphm = ichiph*T_chiph_delta + T_log10_chiph_min;
         logchiphp = logchiphm + T_chiph_delta;

         // Interpolation
         dNBWdt = (T_table[ichiph+1]*fabs(logchiph-logchiphm) +
                 T_table[ichiph]*fabs(logchiphp - logchiph))*T_chiph_inv_delta;
     }
     return factor_dNBWdt*dNBWdt*chiph/gamma;
}

// -----------------------------------------------------------------------------
//! Computation of the value T(chiph) using the approximated
//! formula of Erber
//! \param chiph photon quantum parameter
//! \param nbit number of iteration for the Bessel evaluation
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_Erber_T(double chiph,int nbit,
                                                    double eps)
{
    // Values for Bessel results
    //double I,dI;
    double K;

    //userFunctions::modified_bessel_IK(1./3.,4./(3.*chiph),I,dI,K,dK,nbit,eps);
    K = userFunctions::modified_bessel_K(1./3.,4./(3.*chiph),nbit,eps);

    return 0.16*K*K/chiph;
}

// -----------------------------------------------------------------------------
//! Computation of the value T(chiph) using the formula of Ritus
//! \param chiph photon quantum parameter
//! \param nbit number of iteration for the Gauss-Legendre integration
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_Ritus_T(double chiph,
                      int nbit,double eps)
{

    // Parameters:
    // Arrays for Gauss-Legendre integration
    double * gauleg_x = new double[nbit];
    double * gauleg_w = new double[nbit];
    int i;
    double u;
    double T;

    // gauss Legendre coefficients
    userFunctions::gauss_legendre_coef(log10(0.5*chiph)-50.,log10(0.5*chiph),
            gauleg_x, gauleg_w, nbit, eps);

    T = 0;
    for(i=0 ; i< nbit ; i++)
    {
        u = pow(10.,gauleg_x[i]);
        T += gauleg_w[i]*compute_Ritus_dTdchi(chiph,u,nbit,eps)*u*log(10);
    }
    return T*2;

}

// -----------------------------------------------------------------------------
//! Computation of the value dT/dchipa(chiph) using the formula of Ritus
//! \param chiph photon quantum parameter
//! \param chipa particle quantum parameter
//! \param nbit number of iteration for the Gauss-Legendre integration
//! \param eps epsilon for the Bessel evaluation
// -----------------------------------------------------------------------------
double MultiphotonBreitWheelerTables::compute_Ritus_dTdchi(double chiph,
                      double chipa,int nbit,double eps)
{

    // Parameters:
    // Arrays for Gauss-Legendre integration
    double * gauleg_x = new double[nbit];
    double * gauleg_w = new double[nbit];
    double y,u;
    double p1,p2;
    int i;
    //double I,dI;
    //double K,dK;

    y = (chiph/(3.*chipa*(chiph-chipa)));

    //userFunctions::modified_bessel_IK(2./3.,2.*y,I,dI,K,dK,500,1e-16);
    //p1 = (2. - 3.*chiph*y)*K;

    p1 = (2. - 3.*chiph*y)*userFunctions::modified_bessel_K(2./3.,2*y,500,1e-16);

    // gauss Legendre coefficients
    userFunctions::gauss_legendre_coef(log10(2*y),log10(2*y)+50.,
            gauleg_x, gauleg_w, nbit, eps);

    // Integration loop
    p2 = 0;
    //#pragma omp parallel for reduction(+:p2) private(u) shared(gauleg_w,gauleg_x)
    for(i=0 ; i< nbit ; i++)
    {
        u = pow(10.,gauleg_x[i]);
        p2 += gauleg_w[i]*userFunctions::modified_bessel_K(1./3.,u,500,1e-16)*u*log(10.);
        //userFunctions::modified_bessel_IK(1./3.,u,I,dI,K,dK,500,1e-16);
        //p2 += gauleg_w[i]*K*u*log(10.);
    }

    return (p2 - p1)/(M_PI*sqrt(3.)*pow(chiph,2));
}

// -----------------------------------------------------------------------------
// TABLE COMPUTATION
// -----------------------------------------------------------------------------

//! Computation of the table T_table that is a discetization of the T function
//! for the multiphoton Breit-Wheeler process
//! \param smpi Object of class SmileiMPI containing MPI properties
void MultiphotonBreitWheelerTables::compute_T_table(SmileiMPI *smpi)
{
    // Parameters
    int rank; // Rank number
    // timers
    double t0,t1;
    // tabel_exists
    bool table_exists;

    // Get the MPI rank
    rank = smpi->getRank();

    MESSAGE("        --- T table:");

    // Initial timer
    t0 = MPI_Wtime();

    // If external tables are available, we read them
    table_exists = MultiphotonBreitWheelerTables::read_T_table(smpi);

    // Else we compute them
    if (!table_exists)
    {

        // Temporary photon chi value
        double chiph;
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

        // Allocation of the array T_table
        T_table.resize(T_dim);

        // Allocation of the table for load repartition
        imin_table = new int[nb_ranks];
        length_table = new int[nb_ranks];

        // Computation of the delta
        T_chiph_delta = (log10(T_chiph_max)
                - T_log10_chiph_min)/(T_dim-1);

        // Inverse delta
        T_chiph_inv_delta = 1./T_chiph_delta;

        // Load repartition
        userFunctions::distribute_load_1d_table(nb_ranks,
                T_dim,
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
            chiph = pow(10.,(imin_table[rank] + i)*T_chiph_delta
                  + T_log10_chiph_min);

            buffer[i] = MultiphotonBreitWheelerTables::compute_Ritus_T(chiph,200,1e-15);
            //buffer[i] = MultiphotonBreitWheelerTables::compute_Erber_T(chiph,500000,1e-10);

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
                &T_table[0], &length_table[0], &imin_table[0],
                MPI_DOUBLE, smpi->getGlobalComm());

        // flag computed at true
        T_computed = true;

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
//! Output the computed tables so that thay can be read at the next run.
//
//! \param params list of simulation parameters
//! \param smpi MPI parameters
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::compute_tables(Params& params, SmileiMPI *smpi)
{
    // These tables are loaded only if if one species has Monte-Carlo Compton radiation
    if (params.hasMultiphotonBreitWheeler)
    {
        MultiphotonBreitWheelerTables::compute_T_table(smpi);
    }
}

// -----------------------------------------------------------------------------
// TABLE OUTPUTS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Ouput in a file of the table values of T for the
//! mutliphoton Breit-Wheeler process
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::output_T_table()
{

    if (output_format == "ascii")
    {
        std::ofstream file;
        file.open(table_path + "/tab_T.dat");

        if (file.is_open()) {

            file.precision(12);

            file << "Multiphoton Breit-Wheeler T table\n";

            file << "Dimension chiph - chiph min - chiph max \n";

            file << T_dim;
            file << " "
                << T_chiph_min << " "
                << T_chiph_max << "\n";;

            // Loop over the table values
            for(int i = 0 ; i < T_dim ; i++)
            {
                file <<  T_table[i] << "\n";
            }

            file.close();
        }
    }
    else if (output_format == "binary")
    {
        std::ofstream file;
        file.open(table_path + "/tab_T.bin",std::ios::binary);

        if (file.is_open()) {

            file.write((char*)&T_dim,sizeof (T_dim));
            file.write((char*)&T_chiph_min, sizeof (double));
            file.write((char*)&T_chiph_max, sizeof (double));

            // Loop over the table values
            for(int i = 0 ; i < T_dim ; i++)
            {
                file.write((char*)&T_table[i], sizeof (double));
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

        buffer = table_path + "/multiphoton_Breit_Wheeler_tables.h5";

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
        dims = T_dim;
        dataspaceId = H5Screate_simple(1, &dims, NULL);

        // Creation of the dataset
        datasetId = H5Dcreate(fileId,
                              "h",
                              H5T_NATIVE_DOUBLE,
                              dataspaceId,
                              H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);

        // Fill the dataset
        H5Dwrite(datasetId, H5T_NATIVE_DOUBLE,
                 H5S_ALL, H5S_ALL, H5P_DEFAULT,
                 &T_table[0]);

        // Create attributes
        H5::attr(datasetId, "chiph_min", T_chiph_min);
        H5::attr(datasetId, "chiph_max", T_chiph_max);
        H5::attr(datasetId, "chiph_dim", T_dim);

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
//! Output the computed tables so that thay can be read at the next run.
//! Table output by the master MPI rank
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::output_tables(SmileiMPI *smpi)
{
    // Sequential output
    if (smpi->isMaster())
    {
        // If tables have been computed, they are output on the disk
        // to be used for the next run
        if (T_computed)
        {
            MultiphotonBreitWheelerTables::output_T_table();
        }
    }
}

// -----------------------------------------------------------------------------
// TABLE READING
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Read the external table T
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
bool MultiphotonBreitWheelerTables::read_T_table(SmileiMPI *smpi)
{
    // Flag database available
    bool table_exists = false;

    // Test if an external table exists, if yes we read the table...
    // Binary table
    if (Tools::file_exists(table_path + "/tab_T.bin"))
    {

        table_exists = true;

        if (smpi->getRank()==0)
        {

            // Reading of the table file
            std::ifstream file;
            file.open(table_path + "/tab_T.bin",std::ios::binary);

            if (file.is_open())
            {

                // Read the header
                file.read((char*)&T_dim,
                          sizeof (T_dim));
                file.read((char*)&T_chiph_min,
                          sizeof (T_chiph_min));
                file.read((char*)&T_chiph_max,
                          sizeof (T_chiph_max));

                // Resize of the array integfochi_table before reading
                T_table.resize(T_dim);

                // Read the table values
                file.read((char*)&T_table[0], sizeof (double)*T_dim);

                file.close();
            }

        }

    }
    // HDF5 format
    else if (Tools::file_exists(table_path + "/multiphoton_Breit_Wheeler_tables.h5"))
    {
         hid_t       fileId;
         hid_t       datasetId;
         std::string buffer;

         if (smpi->getRank()==0)
         {

             buffer = table_path + "/multiphoton_Breit_Wheeler_tables.h5";

             fileId = H5Fopen(buffer.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

             datasetId = H5Dopen2(fileId, "h", H5P_DEFAULT);

             // If this dataset exists, we read it
             if (datasetId > 0)
             {

                 table_exists = true;

                 // First, we read attributes
                 H5::getAttr(datasetId, "chiph_dim", T_dim);
                 H5::getAttr(datasetId, "chiph_min", T_chiph_min);
                 H5::getAttr(datasetId, "chiph_max", T_chiph_max);

                 // Resize of the array integfochi_table before reading
                 T_table.resize(T_dim);

                 // then the dataset
                 H5Dread(datasetId,
                        H5T_NATIVE_DOUBLE, H5S_ALL,
                        H5S_ALL, H5P_DEFAULT,
                        &T_table[0]);

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
        MESSAGE("            Dimension quantum parameter: "
                             << T_dim);
        MESSAGE("            Minimum photon quantum parameter chi: "
                             << T_chiph_min);
        MESSAGE("            Maximum photon quantum parameter chi: "
                             << T_chiph_max);

        // Bcast the table to all MPI ranks
        MultiphotonBreitWheelerTables::bcast_T_table(smpi);
    }

    return table_exists;

}

// -----------------------------------------------------------------------------
// TABLE COMMUNICATIONS
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
//! Bcast of the external table T
//
//! \param smpi Object of class SmileiMPI containing MPI properties
// -----------------------------------------------------------------------------
void MultiphotonBreitWheelerTables::bcast_T_table(SmileiMPI *smpi)
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
        MPI_Pack_size(T_dim, MPI_DOUBLE, smpi->getGlobalComm(),
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
        MPI_Pack(&T_dim,
             1,MPI_INTEGER,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&T_chiph_min,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());
        MPI_Pack(&T_chiph_max,
             1,MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

        MPI_Pack(&T_table[0],T_dim,
            MPI_DOUBLE,buffer,buf_size,&position,smpi->getGlobalComm());

    }

    // Bcast all parameters
    MPI_Bcast(&buffer[0], buf_size, MPI_PACKED, 0,smpi->getGlobalComm());

    // Other ranks unpack
    if (smpi->getRank() != 0)
    {
        position = 0;
        MPI_Unpack(buffer, buf_size, &position,
                   &T_dim, 1, MPI_INTEGER,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &T_chiph_min, 1, MPI_DOUBLE,smpi->getGlobalComm());
        MPI_Unpack(buffer, buf_size, &position,
                   &T_chiph_max, 1, MPI_DOUBLE,smpi->getGlobalComm());

        // Resize table before unpacking values
        T_table.resize(T_dim);

        MPI_Unpack(buffer, buf_size, &position,&T_table[0],
                    T_dim, MPI_DOUBLE,smpi->getGlobalComm());

    }

    T_log10_chiph_min = log10(T_chiph_min);

    // Computation of the delta
    T_chiph_delta = (log10(T_chiph_max)
            - T_log10_chiph_min)/(T_dim-1);

    // Inverse delta
    T_chiph_inv_delta = 1./T_chiph_delta;
}
