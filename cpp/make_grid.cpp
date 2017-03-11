/*! \file

Contains the main function for the simuated photometry by integrating code (make-grid)

Integrates filters against template SEDs at specified redshifts and with specified reddenings.

Usage: 

make-grid <list of SEDs> <filters> <model grid spec> <output file>

list of SEDs: A list of input template SEDs to use when building the grid
filters: A list of filters to apply to the template SEDs
model grid spec: specifies redshift, reddening, redenning law.
output: the name of the output file containing the model grid. If the filename ends with .photoz, the
grid is stored in a compact binary format. If the filename ends with .template, only the modified
templates (reddened and redshifted with scaled emission lines) are saved. In all other cases, the
grid is stored as a human-readable (ASCII) file.

\todo rlaw should be an enumeration

something like the following should be introduced to replace the sed_integral arrays

    struct Sed
    {
        struct MetaData
         {
             real_t wavelength;
             real_t frequency;
             real_t frequency_step;
         };

         //Sed(const uint n_seds, const uint n_bins)

         VECTOR<MetaData> info;
         VECTOR<VECTOR<real_t>> sed;
    };
**/

#include "grid_tools.hpp"

#include <boost/algorithm/string.hpp>

#include <omp.h>

#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>


int main(int argc, char* argv[])
{
    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " SED_list Filter_list "
                     "input_catalog output_file\n";
        return -1;
    }

    const uint n_columns_in_catalog = 10; // the expected number of columns in the catalog

    const char* templatelist_filename = argv[1];
    const char* filterlist_filename = argv[2];
    const char* input_catalog_file = argv[3];
    const char* output_file = argv[4];

    std::string outfile(output_file);

    if (outfile.size() < 10)
    {
        std::cerr << "The output filename must consist of at least 10 characters.\n";

        return -1;
    }

    // does the user want to output a binary grid file?

    const bool output_binary = (outfile.substr(outfile.size() - 7 , 7) == ".photoz") ? true : false;

    // does the user just want to output the templates?

    const bool output_template = (outfile.substr(outfile.size() - 9 , 9) == ".template") ? true : false;

    ResampledData output;

    // resample the filters and templates. Results are stored in output object.

    ResampleFiltersAndTemplates(filterlist_filename, templatelist_filename, output);

    // read the catalog contents

    StringList catalog = ReadLines(input_catalog_file);

    // allocate memory for the input data

    real_t** input_data = new real_t*[catalog.size()];

    for (uint i = 0; i < catalog.size(); ++i)
    {
        input_data[i] = new real_t [n_columns_in_catalog];
    }

    // convert catalog entries to double precision floating point and store in input_data[][]

    for (uint i = 0; i < catalog.size(); ++i)
    {
        std::string line = catalog[i];

        StringList column;

        boost::split(column, line,
                     boost::algorithm::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        // we expect each row of the catalog to have n_columns_in_catalog entries. Refuse to continue
        // if this expectation is not met.

        if (column.size() != n_columns_in_catalog)
        {
            std::cerr << "The catalog is corrupt. Each line is supposed to have " << n_columns_in_catalog << " columns.";
            exit(2);
        }

        for (uint j = 0; j < n_columns_in_catalog; ++j)
        {
            input_data[i][j] = std::stod(column[j]);
        }
    }

    // if we just want to output the template corresponding to the input lines:

    if (output_template)
    {
        WriteTemplates(outfile, output, catalog.size(), input_data);

        return 0;
    }


    // Do integrals and get photometry

    std::cerr << "Calculating photometry.\n";

    // allocate memory for fluxes

    real_t** flux = new real_t* [catalog.size()];

    for (int i = 0; i < catalog.size(); ++i)
    {
        flux[i] = new real_t[output.n_filters];
    }

    const uint number_of_threads = omp_get_max_threads();

    // A cache array is created here so that integrate_sed() doesn't have to repeatedly allocate and
    // free memory. We make sure that each sub-vector is allocted by the thread that will be using
    // (this is important for efficient use of cache and NUMA machines).

    VECTOR<VECTOR<real_t>> cache(number_of_threads);

    #pragma omp parallel for schedule(static, 1)

    for (uint i = 0; i < number_of_threads; ++i)
    {
        const uint thread_id = omp_get_thread_num();

        VECTOR<real_t> x(output.n_wave, 0);
        cache[thread_id].swap(x);
    }

    #pragma omp parallel for schedule(dynamic)

    for (int i = 0; i < catalog.size(); ++i) // for each input line
    {
        const real_t redshift = input_data[i][4];
        const int sed_type = input_data[i][5];
        const real_t Ebv = input_data[i][6];
        const int rlaw = input_data[i][7];   // extinction law
        const float line_correction = input_data[i][9]; // emission line scaling factor

        if ((sed_type > 0) && (redshift >= 0) && (rlaw > 0) && (Ebv >= 0))
        {
            const uint thread_id = omp_get_thread_num();

            double* sed_fluxes = output.sed_integral_flux[sed_type - 1];

            // scale the emission lines if the line correction factor is not equal to 1.0.
            /// \todo: make this clean up after itself, using RAII or create a special cache for
            /// the sed_fluxes array.

            if (line_correction != 1.0) // the for loop is outside of the function in order to prevent
            {                           // function calls when the function is not needed.
                sed_fluxes = ScaleEmissionLines(output.n_wave,
                                                output.sed_integral_flux[sed_type - 1],
                                                output.sed_integral_wave,
                                                line_correction);
            }

            integrate_sed(output.dlogwave,             // log wavelength spacing
                          output.sed_integral_wave,    // SED wavelengths (A)
                          output.dnu_nu.data(),        // sed_integral_dnu[i] / sed_integral_nu[i];
                          sed_fluxes,  // SED fluxes (Jy)
                          output.n_wave,               // number of points?
                          output.filter_n.data(),      // number of elements per filter
                          output.filter_starting_index.data(), // starting index per filter
                          output.filter_transmissions, // filter transmissions
                          output.n_filters,            // number of filters
                          output.red_Data[rlaw - 1],   // the associated reddening vector
                          redshift,             // redshift
                          Ebv,                  // E(B-V)
                          cache[thread_id].data(),
                          flux[i]);             // output fluxes (Jy)

            // clean up the array allocated by ScaleEmissionLines()
            if (sed_fluxes != output.sed_integral_flux[sed_type - 1])
            {
                delete [] sed_fluxes;
            }
        }
    }

    std::cerr << "Writing grid model to " << output_file << std::endl;

    // now that we have calculated the output photometry, we're ready to write it into the output file!

    if (output_binary)
    {
        SavePhotoZBinary(output_file, output.n_filters, catalog.size(), input_data, flux);
    }
    else // output ASCII file.
    {
        std::ofstream output_catalog(output_file);

        if (!output_catalog.is_open())
        {
            std::cerr << "Error: Cannot create file, " << output_file << "\n";
            exit(2);
        }

        for (int i = 0; i < catalog.size(); ++i) // for each object (each entry in the catalog)
        {
            // insert the contents of the input catalog into the output catalog---except for the
            // line scale factor

            StringList entry;

            boost::split(entry, catalog[i],
                         boost::algorithm::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            for (uint col = 0; col < 9; ++col) { output_catalog << entry[col] << " "; }

            const float scale_factor = input_data[i][8];

            real_t * const object_flux = flux[i];

            // save convolved photometry

            for (int j = 0; j < output.n_filters; ++j)
            {
                output_catalog << std::scientific << std::setprecision(6) << std::setw(14) << scale_factor * object_flux[j];
            }

            output_catalog << "\n";
        }
    }
}
