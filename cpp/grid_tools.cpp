/*! \file

  Contains the bulk of the code for the simulated photometry / model grid generation utilities.
*/


#include "grid_tools.hpp"
#include "binary_search.hpp"

#include <omp.h>

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>
#include <algorithm>

#define ebv_step 0.055

const uint NFLT = 1953261166; // "nflt" ASCII in little endian encoding, used in the binary file format
const uint ZBIN = 1852400250; // "zbin" ASCII in little endian encoding, used in the binary grid file format


/*! Stores the model grid data from fluxes[][] in a binary file format that can be quickly read by
the photo-z executable. */
void SavePhotoZBinary(const char* filename,
                      const uint n_filters,
                      const uint n_models,
                      double** metadata,
                      double** fluxes)
{
    std::ofstream model_grid_file(filename, std::ofstream::binary);

    if (!model_grid_file.is_open())
    {
        std::cerr << "Error: could not open file " << filename << "\n";
        exit(2);
    }

    // Now we will normalize the models so that the total flux in each model is 1 (because the
    // photo-z code works on the assumption that the models are normalized).

    VECTOR<float> flux(n_filters * n_models, 0);

#pragma omp parallel for schedule(dynamic)

    for (uint model = 0; model < n_models; ++model)
    {
        float total_flux = 0;

        double* input_model_fluxes = fluxes[model];

        for (uint filter = 0; filter < n_filters; ++filter)
        {
            total_flux += input_model_fluxes[filter];
        }

        const float scale = (total_flux > 0) ? 1.0f / total_flux : 0.0f;

        float* model_fluxes = &flux.data()[model * n_filters];

        for (uint filter = 0; filter < n_filters; ++filter)
        {
            model_fluxes[filter] = scale * input_model_fluxes[filter];
        }
    }

    // Insert magic number to indicate that (1) this is a photo-z model grid file and (2), the next
    // four bytes contain an integer specifying number of filters that we are using.

    model_grid_file.write((char*) &NFLT, sizeof(NFLT));

    model_grid_file.write((char*) &n_filters, sizeof(n_filters));

    float prev_redshift_bin = -1;

    const std::streamsize chunk_size_in_bytes = sizeof(float) * n_filters;

    for (uint model = 0; model < n_models; ++model)
    {
        const float model_redshift = metadata[model][4];

        if (model_redshift != prev_redshift_bin)
        {
            // Insert number ZBIN to indicate that the next 4 bytes contain a floating point number
            // which specifies the redshift bin of the following models.

            model_grid_file.write((char*) &ZBIN, sizeof(ZBIN));

            model_grid_file.write((char*) &model_redshift, sizeof(model_redshift));

            prev_redshift_bin = model_redshift;
        }

        // write the template number used by this model

        const int template_id = metadata[model][gridspec_column::template_id];

        model_grid_file.write((char*) &template_id, sizeof(template_id));

        float* model_fluxes = &flux.data()[model * n_filters];

        // write all of the fluxes for the current model:

        model_grid_file.write((char*) model_fluxes, chunk_size_in_bytes);
    }

    model_grid_file.close();
}


VECTOR<float> UniqueRedshiftsInCatalog(const char* catalog_file)
{
    StringList catalog_contents = ReadLines(catalog_file);

    VECTOR<float> all_redshifts;

    for (const auto& line : catalog_contents)
    {
        StringList strings;

        boost::split(strings, line,
                     boost::algorithm::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        const float redshift = std::stof(strings[catalog_column::redshift]);

        if (redshift < 6.5 && redshift >= 0)
        {
            all_redshifts.push_back(redshift);
        }
    }

    std::sort(all_redshifts.begin(), all_redshifts.end());

    VECTOR<float> unique_redshifts;

    unique_redshifts.push_back(all_redshifts.at(0)); // there is at least one uniqe redshift.

    float previous_z = 0;

    for (const auto& z : all_redshifts)
    {
        if (z - previous_z > 0.0005)
        {
            unique_redshifts.push_back(z);
            previous_z = z;
        }
    }

    return unique_redshifts;
}


inline const ChiSquaredResult ChiSquared(const VECTOR<Flux>& observed_fluxes,
                                         const float total_weighted_observed_flux,
                                         const VECTOR<float> &template_fluxes)
{
    // compute the normalization

    float template_sum = 0;

    const uint n_filters = observed_fluxes.size();

    for (uint i = 0; i < n_filters; ++i)
    {
        const Flux ofluxi = observed_fluxes[i];
        template_sum += (ofluxi.value > 0.0f) ? template_fluxes[i] * ofluxi.weight: 0.0f; // only base the sum on the non-zero elements
    }

    const float scale_factor = total_weighted_observed_flux / template_sum;

    float chi2 = 0;

    for (uint i = 0; i < n_filters; ++i)
    {
        const Flux ofluxi = observed_fluxes[i];

        if (ofluxi.value <= 0) { continue; } // skip this iteration; the observed flux is zero

        const float diff = ofluxi.value - scale_factor * template_fluxes[i];

        chi2 += diff * diff * ofluxi.weight;
    }

    return {chi2, scale_factor};
}


void MakePhotometryFittingGrid(const char* catalog_filename,
                               const char* filterlist_filename,
                               const char* templatelist_filename,
                               const char* output_filename)
{
    VECTOR<float> redshifts = UniqueRedshiftsInCatalog(catalog_filename);

    ResampledData output;

    // resample the filters and templates. Results are stored in the 'output' object.
    ResampleFiltersAndTemplates(filterlist_filename,
                                templatelist_filename,
                                output,
                                142, 0.01);

    std::cerr << "There are " << redshifts.size() << " unique redshifts in the catalog.\n";

    std::cerr << "\nIntegrating filters against the templates...\n";

    const int ebv_max = (1 + ebv_step) / ebv_step;

    output.n_ebv = ebv_max;

    const uint grid_size = output.n_template * output.n_rlaw * ebv_max;

    const uint grid_size_per_template = output.n_rlaw * ebv_max;

    VECTOR<VECTOR<real_t>> template_filter_fluxes(grid_size);

    for (auto& fbvec : template_filter_fluxes)
    {
        VECTOR<real_t> blank_vec(output.n_filters);
        fbvec.swap(blank_vec);
    }

    // validate output filename and determine the starting and ending indices of the slice.

    SliceLimits limits(output_filename, redshifts.size());

    std::ofstream fitting_grid(output_filename);

    if (!fitting_grid.is_open())
    {
        std::cerr << "Error: Cannot create file, " << output_filename << "\n";
        exit(2);
    }

    for (uint z_idx = limits.first(); z_idx < limits.last(); ++ z_idx)
    {
        const float redshift = redshifts[z_idx];

        const uint number_of_threads = omp_get_max_threads();

        // cache array so that integrate_sed() doesn't have to repeatedly allocate and free memory.
        // We make sure that each sub-vector is allocted by the thread that will be using it.
        // (this is important for efficient use of CPU cache and NUMA)

        VECTOR<VECTOR<real_t>> cache(number_of_threads);

        #pragma omp parallel for schedule(static, 1)

        for (uint i = 0; i < number_of_threads; ++i)
        {
            const uint thread_id = omp_get_thread_num();

            VECTOR<real_t> c(output.n_wave, 0);
            cache[thread_id].swap(c);
        }

        // loop over, templates, reddening law, ebv:

        #pragma omp parallel for schedule(dynamic)

        for (int template_id = 0; template_id < output.n_template; ++template_id)
        {
            const uint thread_id = omp_get_thread_num();

            for (int rlaw = 0; rlaw < output.n_rlaw; ++rlaw)
            {
                for (int ebv_idx = 0; ebv_idx < ebv_max; ++ebv_idx)
                {
                    real_t ebv = ebv_step * ebv_idx;

                    // run the integrate_sed code and get the observed filter fluxes

                    integrate_sed(output.dlogwave,             // log wavelength spacing
                                  output.sed_integral_wave,    // SED wavelengths (A)
                                  output.dnu_nu.data(),        // sed_integral_dnu[i] / sed_integral_nu[i] : all i;
                                  output.sed_integral_flux[template_id],  // SED fluxes (Jy)
                                  output.n_wave,                // number of points?
                                  output.filter_n.data(),     // number of elements per filter
                                  output.filter_starting_index.data(), // starting index per filter
                                  output.filter_transmissions, // filter transmissions
                                  output.n_filters,            // number of filters
                                  output.red_Data[rlaw],   // the associated reddening vector
                                  redshift,             // redshift
                                  ebv,                  // E(B-V)
                                  cache[thread_id].data(),
                                  template_filter_fluxes[grid_size_per_template * template_id + ebv_max * rlaw + ebv_idx].data()); // output fluxes (Jy)
                }
            }
        }

        // write the fluxes

        for(const auto& output_fluxes : template_filter_fluxes)
        {
            for(const auto flux_value : output_fluxes) { fitting_grid << flux_value << " "; }

            fitting_grid << "\n";
        }
    }

    // write the footer

    if (limits.last() == redshifts.size())
    {
        fitting_grid << output.n_filters << " "
                     << output.n_rlaw << " "
                     << output.n_ebv << " "
                     << output.n_template
                     << std::endl;
    }

    std::cerr << "Finished writing the model grid for fitting.\n";

    fitting_grid.close();
}


void FindBestFits(const char* catalog_filename,
                  const char* model_grid_filename,
                  const char* output_filename)
{
    // The ID, logmass, nbands, RA, DEC, I-band magnitude, and redshift  columns at the beginning
    // of each line. All remaining entries are fluxes through filters.

    const uint n_metadata = 7;

    // identify the unique redshift values in the catalog

    std::cerr << "reading catalog to get unique redshifts\n";
    VECTOR<float> unique_redshifts = UniqueRedshiftsInCatalog(catalog_filename);

    // create the search tree that will be used for finding the index of the nearest redshift bin for
    // each object.

    BinaryTree<float> index_search(unique_redshifts);

    // read the contents of the catalog

    StringList catalog = ReadLines(catalog_filename);

    StringList model_grid = ReadLines(model_grid_filename);

    // read model grid footer. The final line is supposed to contain only:
    //       (0) number of filters
    //       (1) number of reddening laws
    //       (2) number of ebv steps
    //       (3) number of templates

    StringList footer;

    boost::split(footer, model_grid.back(),
                 boost::algorithm::is_any_of("\t "),
                 boost::algorithm::token_compress_on);

    if (footer.size() != 4)
    {
        std::cerr << "The model grid file, " <<  model_grid_filename  << ", has a corrupt footer.\n";
        exit(1);
    }

    const uint n_filters = std::stoi(footer[0]);
    const uint n_rlaw = std::stoi(footer[1]);
    const uint n_ebv = std::stoi(footer[2]);
    const uint n_templates = std::stoi(footer[3]);

    // remove footer from the model grid

    model_grid.pop_back();

    std::cerr << "\nNow we search for the optimal set of parameters for each object in the catalog.\n";

    // identify the bounding indices for the slice of interest

    SliceLimits slice_index(output_filename, catalog.size());

    std::ofstream best_fits(output_filename);

    if (!best_fits.is_open())
    {
        std::cerr << "Error: Cannot create file, " << output_filename << "\n";
        exit(2);
    }

    for (uint line = slice_index.first(); line < slice_index.last(); ++line)
    {
        std::string& line_contents = catalog[line];
        StringList col;

        boost::split(col, line_contents,
                     boost::algorithm::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        const uint apparent_number_of_filters_in_catalog = (col.size() - n_metadata) / 2;

        if ( apparent_number_of_filters_in_catalog != n_filters)
        {
            std::cerr << "Error: The catalog uses " <<  apparent_number_of_filters_in_catalog
                      << " filters, but there are " << n_filters << " filters in "
                      << model_grid_filename << ".\n";

            exit(2);
        }

        // the redshift of the current object:

        const real_t redshift = std::stod(col[catalog_column::redshift]);

        // the starting major index in the model grid

        const uint redshift_index = index_search.index_of(redshift);

        const uint models_per_template = n_rlaw * n_ebv;

        const uint models_per_redshift = n_templates * models_per_template;

        // The fist entry in the model_grid that corresponds to the appropriate redshift

        const uint first_index = models_per_redshift * redshift_index;

        // an array for storing the current object's fluxes from the catalog:

        VECTOR<Flux> observed_filter_fluxes(n_filters, {0.0f, 0.0f});

        float weighted_flux = 0;

        for (uint i = n_metadata; i < col.size() - n_filters; ++i)
        {
            const float flux = std::stod(col[i]);
            const float uncertainty = std::stod(col[i + n_filters]);
            const float weight = 1.0f / (uncertainty * uncertainty);

            weighted_flux += flux * weight;

            observed_filter_fluxes[i - n_metadata] = { flux , weight };
        }

        VECTOR<OptimalParams> template_chi2(n_templates);

        const uint number_of_threads = omp_get_max_threads();

        VECTOR<VECTOR<float>> template_filter_fluxes(number_of_threads);

        #pragma omp parallel for schedule(static, 1)

        for (uint i = 0; i < number_of_threads; ++i)
        {
            const uint thread_id = omp_get_thread_num();

            VECTOR<float> tff(n_filters, 0);
            template_filter_fluxes[thread_id].swap(tff);
        }

        // loop over templates, reddening law, ebv:

        #pragma omp parallel for schedule(dynamic)

        for (int template_id = 0; template_id < n_templates; ++template_id)
        {
            const uint thread_id = omp_get_thread_num();

            // set up quantities for identifying the optimal parameters

            ChiSquaredResult min_rlaw_chi2 = { 2e30, 1 };
            int optimal_rlaw = -9;

            ChiSquaredResult min_ebv_chi2 = { 1e30, 1 };
            float optimal_ebv = -9;

            for (int rlaw = 0; rlaw < n_rlaw; ++rlaw)
            {
                for (int ebv_idx = 0; ebv_idx < n_ebv; ++ebv_idx)
                {
                    real_t ebv = ebv_step * ebv_idx;

                    const uint offset = models_per_template * template_id + n_ebv * rlaw + ebv_idx;

                    const uint grid_idx = first_index + offset;

                    StringList fluxes;

                    boost::split(fluxes, model_grid[grid_idx],
                                 boost::algorithm::is_any_of("\t "),
                                 boost::algorithm::token_compress_on);

                    VECTOR<float>& tff = template_filter_fluxes[thread_id];

                    for (uint fidx = 0; fidx < n_filters; ++fidx)
                    {
                        tff[fidx] = std::stod(fluxes[fidx]);
                    }

                    // do chi2 comparison between input fluxes and template fluxes

                    ChiSquaredResult chi2 = ChiSquared(observed_filter_fluxes, weighted_flux, tff);

                    // identify the value of ebv which minimized chi2

                    if (chi2.value < min_ebv_chi2.value)
                    {
                        optimal_ebv = ebv;
                        min_ebv_chi2 = chi2;
                    }
                }

                // identify value of rlaw which minimized chi2

                if (min_ebv_chi2.value < min_rlaw_chi2.value)
                {
                    optimal_rlaw = rlaw;
                    min_rlaw_chi2 = min_ebv_chi2;
                }
            }

            // store min_rlaw_chi2 for this template

            OptimalParams optimal_params = {min_rlaw_chi2, optimal_ebv, optimal_rlaw};

            template_chi2[template_id] = optimal_params;
        }

        // identify the optimal template

        ChiSquaredResult min_template_chi2 = { 3e30, 1 };
        int best_template = -9;

        for (int template_id = 0; template_id < n_templates; ++template_id)
        {
            const OptimalParams& best = template_chi2[template_id];

            if (best.chi2.value < min_template_chi2.value)
            {
                best_template = template_id;
                min_template_chi2 = best.chi2;
            }
        }

        const OptimalParams& optimal = template_chi2[best_template];

        // the first 5 entries of the output file are: ID, RA, DEC, I-band magnitude, and redshift

        best_fits << std::setw(9) << std::right << std::stoul(col[catalog_column::id])
                  << std::setprecision(5) << std::fixed
                  << std::setw(11) << std::stod(col[catalog_column::RA]) << " "
                  << std::setw(9) << std::stod(col[catalog_column::DEC])
                  << std::setw(10) << std::stod(col[catalog_column::I_magnitude])
                  << std::setprecision(4) << std::fixed
                  << std::setw(8) << std::stod(col[catalog_column::redshift]);

        // now add the template number, reddening, reddening law, and scale factor

        best_fits << std::setw(5) << best_template + 1
                  << std::setw(8) << optimal.ebv
                  << std::setw(3) << optimal.rlaw + 1
                  << std::setprecision(6) << std::scientific
                  << std::setw(13) << optimal.chi2.scale_factor << "\n";
    }

    best_fits.close();
}


void ResampleFiltersAndTemplates(const char* filterlist,
                                 const char* templatelist,
                                 ResampledData& output,
                                 const uint Rmax,
                                 const float dzmin)
{
    VECTOR<VECTOR<double>> template_fluxes;
    VECTOR<VECTOR<double>> template_wavelengths;
    VECTOR<VECTOR<double>> filter_transmissions;
    VECTOR<VECTOR<double>> filter_wavelengths;

    // read all of the filters

    ReadTwoColumnFiles(filterlist, filter_wavelengths, filter_transmissions);

    // read all of the templates

    ReadTwoColumnFiles(templatelist, template_wavelengths, template_fluxes);

    const int n_filters = filter_wavelengths.size();
    const int n_template = template_wavelengths.size();

    std::cerr << "There are " << n_filters << " filters and " << n_template << " templates\n";
    std::cerr << "Resampling the filters and templates so that they share common "
                 "wavelength bins\n";

    ////////////// transform data into the format that is needed by set_sampling() /////////////////

    VECTOR<double*> wavedata(n_template, 0);
    VECTOR<double*> fluxdata(n_template, 0);
    VECTOR<int> n_elements(n_template, 0);

    for (uint i = 0; i < n_template; ++i)
    {
        wavedata[i] = template_wavelengths[i].data();
        fluxdata[i] = template_fluxes[i].data();
        n_elements[i] = template_wavelengths[i].size();
    }

    double** sed_waveData = wavedata.data();

    double** sed_fluxData = fluxdata.data();

    int* sed_Nelem = n_elements.data();

    int Nwave;

    double lmin;

    double dlogwave;

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // prepare to resample the SEDs

    // sed_waveData (double**) : SED wavelengths
    // sed_Nelem (int*) number of elements per SED
    // n_template (int) number of seds
    // dzmin (float) the redshift interval in the grid (because this was designed for making a grid)
    // Rmax (float) the maximum resolution of the filters
    // output: Nwave (int*) the number of wavelengths that should be sampled
    // output: lmin (real_t*) minimum wavelength
    // output: dlogwave (real_t*) log wavelength spacing

    set_sampling(sed_waveData, sed_Nelem, n_template, dzmin, Rmax, &Nwave, &lmin, &dlogwave);

    output.n_wave = Nwave;
    output.n_filters = n_filters;
    output.n_template = n_template;
    output.dlogwave = dlogwave;

    ///////////////////////// allocate quantities needed by resample_seds() ////////////////////////

    output.sed_integral_wave = new real_t[Nwave];

    real_t* sed_integral_nu = new real_t[Nwave];
    real_t* sed_integral_dnu = new real_t[Nwave];

    output.sed_integral_flux = new real_t*[n_template];

    for (int i = 0; i < n_template; ++i)
    {
        output.sed_integral_flux[i] = new real_t[Nwave];
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // sed_waveData ''
    // sed_fluxData (double**) sed fluxes
    // sed_Nelem ''
    // n_template ''
    // output sed_integral_wave (real_t*) output wavelengths
    // output sed_integral_nu (real_t*) output frequencies
    // output sed_integral_dnu (real_t*) output frequency spacings
    // output sed_integral_flux (real_t*) output fluxes
    // Nwave ''
    // lmin ''
    // dlogwave ''

    resample_seds(sed_waveData,
                  sed_fluxData,
                  sed_Nelem,
                  n_template,
                  output.sed_integral_wave, // these should probably be replaced with a single Sed object
                  sed_integral_nu,   // ''
                  sed_integral_dnu,  // ''
                  output.sed_integral_flux, // ''
                  Nwave,
                  lmin,
                  dlogwave);

    // create a vector to store the ratio, dnu / nu in order to speed up the integration process.
    VECTOR<float> dnu_nu(Nwave, 0.0f);

    for (uint i = 0; i < Nwave; ++i)
    {
        dnu_nu[i] = sed_integral_dnu[i] / sed_integral_nu[i];
        if (std::isnan(dnu_nu[i]))
        {
            std::cerr << "sed_integral_nu[i] = " << sed_integral_nu[i] << "\n";
            std::cerr << "sed_integral_dnu[i] = " << sed_integral_dnu[i] << "\n";
        }
    }

    output.dnu_nu.swap(dnu_nu);

    // allocate memory for reddening vectors

    const int Nrlaw = 5;

    output.n_rlaw = Nrlaw;

    output.red_Data = new real_t*[Nrlaw];

    for (int i = 0; i < Nrlaw; ++i)
    {
        output.red_Data[i] = new real_t[Nwave];
    }

    // make the appropriate reddening vectors that can be subsequently applied to the data.

    make_reddening_vec(output.sed_integral_wave, output.red_Data, Nrlaw, Nwave);

    /////////////// setup the quantities needed by setup_resample_filters() ////////////////////////

    VECTOR<int> filter_starting_index(n_filters, 0);
    VECTOR<int> filter_n(n_filters, 0);

    VECTOR<double*> fwavedata(n_filters, 0);
    VECTOR<double*> ffluxdata(n_filters, 0);
    VECTOR<int> fn_elements(n_filters, 0);

    for (uint i = 0; i < n_filters; ++i)
    {
        fwavedata[i] = filter_wavelengths[i].data();
        ffluxdata[i] = filter_transmissions[i].data();
        fn_elements[i] = filter_wavelengths[i].size();
    }

    double** filter_waveData = fwavedata.data();

    double** filter_fluxData = ffluxdata.data();

    int* filter_Nelem = fn_elements.data();

    ////////////////////////////////////////////////////////////////////////////////////////////////

    // sed_integral_wave (real_t*) ''
    // lmin (real_t) ''
    // dlogwave (real_t) ''
    // filter_waveData (real_t**) filter wavelengths
    // filter_fluxData (real_t**) filter transmission
    // n_filters (int) ''
    // filter_Nelem (int*) number of wavelengths per filter
    // output filter_starting_index.data() (int*) the starting index of the filter
    // output filter_n.data() number of elements in the output / resampled filter

    setup_resample_filters(output.sed_integral_wave,
                           lmin,
                           dlogwave,
                           filter_waveData,
                           filter_fluxData,
                           n_filters,
                           filter_Nelem,
                           filter_starting_index.data(),
                           filter_n.data());

    output.filter_transmissions = new real_t* [n_filters];

    for (int i = 0; i < n_filters; i++)
    {
        output.filter_transmissions[i] = new real_t [filter_n.data()[i]];
    }

    // resample the filters
    // inputs are output SED wavelengths (A), output filter wavelengths (double array; A), output
    // filter transmission (double array), number of filters, number of wavelengths per filter,
    // number of wavelengths down to start sampling each filter from, output number of elements per
    // filter, output filter transmissions (double array)

    resample_filters(output.sed_integral_wave,
                     filter_waveData,
                     filter_fluxData,
                     n_filters,
                     filter_Nelem,
                     filter_starting_index.data(),
                     filter_n.data(),
                     output.filter_transmissions);

    output.filter_starting_index.swap(filter_starting_index);
    output.filter_n.swap(filter_n);
}

void WriteTemplates(const std::string output_file,
                    const ResampledData &resampled_templates,
                    const int n_entries,
                    real_t** catalog_data)
{

    VECTOR<real_t> template_fluxes(resampled_templates.n_wave, 0);

    std::ofstream outfile(output_file);

    if (!outfile.is_open())
    {
        std::cerr << "Error: Cannot create file, " << output_file << "\n";
        exit(2);
    }

    for (int i = 0; i < n_entries; ++i) // for each input line of the catalog
    {
        const real_t redshift = catalog_data[i][4];
        const int sed_type = catalog_data[i][5];
        const real_t Ebv = catalog_data[i][6];
        const int rlaw = catalog_data[i][7];
        const float scale_factor = catalog_data[i][8];
        const float line_correction = catalog_data[i][9]; // emission line scaling factor

        if ((sed_type > 0) && (redshift >= 0) && (rlaw > 0) && (Ebv >= 0))
        {
            double* sed_fluxes = resampled_templates.sed_integral_flux[sed_type - 1];

            // scale the emission lines if the line correction factor is not equal to 1.0.
            /// \todo: make this clean up after itself, using RAII or create a special cache for
            /// the sed_fluxes array.

            if (line_correction != 1.0) // the for loop is outside of the function in order to prevent
            {                           // function call overhead when the function is not needed.
                sed_fluxes = ScaleEmissionLines(resampled_templates.n_wave,
                                                resampled_templates.sed_integral_flux[sed_type - 1],
                                                resampled_templates.sed_integral_wave,
                                                line_correction);
            }

            ComputeSed(resampled_templates.sed_integral_wave,    // SED wavelengths (A)
                       sed_fluxes,  // SED fluxes (Jy)
                       resampled_templates.n_wave,               // number of points?
                       resampled_templates.red_Data[rlaw - 1],   // the associated reddening vector
                       redshift,              // redshift
                       Ebv,                   // E(B-V)
                       template_fluxes.data());

            // write wavelengths

            const float zplusone = 1.0 + redshift;

            for (uint i = 0; i < resampled_templates.n_wave; ++i)
            {
                outfile << zplusone * resampled_templates.sed_integral_wave[i] << " ";
            }

            outfile << "\n";

            for (const real_t flux : template_fluxes)
            {
                outfile << flux * scale_factor << " "; // why is the -1 needed?
            }

            outfile << "\n";

            // clean up the array allocated by ScaleEmissionLines()
            if (sed_fluxes != resampled_templates.sed_integral_flux[sed_type - 1])
            {
                delete [] sed_fluxes;
            }
        }
        else
        {
            std::cerr << "some entries in the input catalog are out of bounds.\n";
        }
    }

    outfile.close();
}


double* ScaleEmissionLines(const uint n_wavelengths,
                           const double *fluxes,
                           const double *wavelengths,
                           const float line_correction)
{
    double* corrected_fluxes = new double[n_wavelengths];

    for (uint i = 0; i < n_wavelengths; ++i)
    {
        const float lambda = wavelengths[i];

        const bool is_line1 = (lambda > 3715.0f && lambda < 3742.0f); // [OII]3727
        const bool is_line2 = (lambda > 4850.0f && lambda < 4873.0f); // H-beta
        const bool is_line3 = (lambda > 4950.0f && lambda < 4970.0f); // [OIII]4959
        const bool is_line4 = (lambda > 4995.0f && lambda < 5020.0f); // [OIII]5007
        const bool is_line5 = (lambda > 6540.0f && lambda < 6595.0f); // [OIII]5007

        if (is_line1 || is_line2 || is_line3 || is_line4 || is_line5)
        {
            corrected_fluxes[i] = line_correction * fluxes[i];
            //std::cout << "scaling by " << line_correction << "\n old = " << fluxes[i] << ", new = " << corrected_fluxes[i] << "\n";
        }
        else
        {
            corrected_fluxes[i] = fluxes[i];
        }
    }

    return corrected_fluxes;
}


