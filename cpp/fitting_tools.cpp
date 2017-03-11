/*! \file
 *
 * Contains the bulk of the actual redshift fitting code.
 */

#include "fitting_tools.hpp"
#include "basic_utils.hpp"

#include <boost/algorithm/string.hpp>

#include <omp.h>

#include <fstream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <cmath>

using namespace std;

// constants used by the .photoz model grid format:

const uint NFLT = 1953261166; // "nflt" ASCII in little endian encoding
const uint ZBIN = 1852400250; // "zbin" ASCII in little endian encoding


uint RedshiftPDF::ReadModelGrid(const char *model_grid_file_name)
{
    std::cerr << "Reading " << model_grid_file_name << std::endl;

    std::ifstream model_grid_file(model_grid_file_name, std::ifstream::binary);

    if (!model_grid_file.is_open())
    {
        std::cerr << "Error: could not open file " << model_grid_file_name << "\n";
        exit(3); /// \todo In general, error reporting should be handled more consistently / systematically
                 /// with well-documented, specific error codes.
    }

    // get the size of the file, in bytes

    model_grid_file.seekg(0, model_grid_file.end);
    const unsigned long file_size = model_grid_file.tellg();
    model_grid_file.seekg (0, model_grid_file.beg);

    // If the file is smaller than 16 bytes, it can't possibly be a valid file. Realistically, it
    // it will be more than 1 GB in size, but the format can't possibly correct with less than 16 bytes.

    if (file_size < 16)
    {
        std::cerr << "The file " << model_grid_file_name << " is not large enough to be a real "
                     "model grid. The file may be corrupt. Size = " << file_size << " bytes.\n";
        exit(4);
    }

    uint magic_number = 0;

    // the first entry in the file must be NFLT, in order for the file to be valid.

    model_grid_file.read((char*) &magic_number, sizeof magic_number);

    if (magic_number != NFLT)
    {
        std::cerr << "Error: " << model_grid_file_name << " is not a valid photo-z model grid file.\n";
        exit(4);
    }

    // get the number of filters

    uint n_filters = 0;

    model_grid_file.read((char*) &n_filters, sizeof n_filters);

    // perform a basic check of the file size. It has to be at least large enough to hold the magic
    // numbers, the fluxes, the redshift bin, the model metadata, and the n_filters specification

    if (file_size < sizeof(float) * (n_filters + 2) + kn_metadata_ * sizeof(int) + 2 * sizeof(uint))
    {
        std::cerr << "The file " << model_grid_file_name << " is corrupt.\n";
        exit(4);
    }

    float* buffer = new float[n_filters + kn_metadata_];

    ModelSeds* current_redshift_vector;

    float current_redshift = -2.0f;

    while (model_grid_file.tellg() < file_size - 5)
    {
        // the position in the stream, as we enter this iteration of the while loop
        const auto position = model_grid_file.tellg();

        uint next_number;

        model_grid_file.read((char*) &next_number, sizeof next_number);

        if (next_number == ZBIN)
        {
            // the next four bytes are supposed to contain a redshift value

            float redshift;

            model_grid_file.read((char*) &redshift, sizeof redshift);

            if (redshift > 100.0f || redshift < 0.0f)
            {
                cerr << "Input file error: Redshift out of range. z = " << redshift << " in file, "
                     << model_grid_file_name << ", at byte = "
                     << model_grid_file.tellg() - (std::ifstream::pos_type) (sizeof (float)) << "\n";

                exit(4);
            }

            if (current_redshift != redshift)
            {
                // we're working with a new redshift bin.
                current_redshift = redshift;
                current_redshift_vector = &model_grid_[current_redshift];
            }
        }
        else // set position back to where it was before we checked for ZBIN
        {
            model_grid_file.seekg(position);
        }

        // read fluxes for one model SED

        model_grid_file.read((char*) buffer, sizeof(float) * (n_filters + kn_metadata_));

        // copy the fluxes from the buffer to the end of the vector

        current_redshift_vector->insert(current_redshift_vector->end(),
                                        buffer,
                                        buffer + n_filters + kn_metadata_);
    }

    model_grid_file.close();

    delete [] buffer;

    return n_filters;
}


uint RedshiftPDF::ReadSourceSeds(const char* object_sed_filename)
{
    uint n_filters = 0;

    ifstream object_sed_ifs(object_sed_filename);

    if (!object_sed_ifs.is_open())
    {
       std::cerr << "Error: " << object_sed_filename << "could not be opened\n";
       exit(1);
    }
    else
    {
       std::cerr << "Reading " << object_sed_filename << endl;

       string line;

       while (std::getline(object_sed_ifs, line))
       {
          // trim whitespace at the beginning and end of the line then split the line into a
          // vector of strings, delimited by whitespace.

          boost::trim(line);

          VECTOR<std::string> strs;

          boost::split(strs, line,
                       boost::algorithm::is_any_of("\t "),
                       boost::algorithm::token_compress_on);

          SourceInfo current_object;

          // the first 3 entries are the object ID, x pixel coordinate, y pixel coordinate

          current_object.id = std::stoi(strs[0]);
          current_object.xpos = std::stod(strs[1]);
          current_object.ypos = std::stod(strs[2]);

          // the remaining entries are the sed: flux_0, uncertainty_0, flux_1, uncertainty_1, ...

          auto& fluxes = current_object.sed;

          float total_sed_flux = 0;

          float total_noise = 0;

          float snr = 0;

          float max_flux = 0;

          for (uint i = 3; i < strs.size(); i += 2)
          {
              const float fluxi = std::stod(strs[i]);
              const float uncertaintyi = std::stod(strs[i + 1]);

              // ith weight, weighti = 1.0 / uncertaintyi ^ 2

              const float w1 = 1.0f / uncertaintyi;
              const float weighti = w1 * w1;

              // compute contribution to the weighted_flux for the current object

              total_sed_flux += fluxi;

              total_noise += uncertaintyi * uncertaintyi;

              snr += fluxi * w1;

              max_flux = (fluxi < max_flux) ? max_flux : fluxi;

              // add the ith Flux element to current_object.sed

              fluxes.push_back( {fluxi, weighti} );
          }

          // compute the total weighted flux

          float weighted_sed_flux = 0;
          float total_weight = 0;
          float total_squared_flux = 0;
          float total_weighted_squared_flux = 0;

          for (auto& flux : fluxes)
          {
              const float squared_flux = flux.value * flux.value;

              weighted_sed_flux += flux.value * flux.weight;
              total_weight += flux.weight;
              total_squared_flux += squared_flux;
              total_weighted_squared_flux += squared_flux * flux.weight;
          }

          // store the total weighted flux and the total raw flux

          current_object.weighted_flux = weighted_sed_flux;

          current_object.total_flux = total_sed_flux;

          current_object.total_squared_flux = total_squared_flux;

          current_object.total_weighted_squared_flux = total_weighted_squared_flux;

          current_object.total_weight = total_weight;

          current_object.total_noise = total_noise;

          current_object.snr = snr;

          object_list_.push_back(current_object);

          if (n_filters == 0) { n_filters = strs.size() - 3; }
       }

       object_sed_ifs.close();
    }

    return n_filters / 2; // divide by two because the line contains (flux, error) pairs.
}


VECTOR<float> RedshiftPDF::ReadTemplateProbabilities(const char *filename)
{
    StringList lines = ReadLines(filename);

    const uint n_templates = lines.size();

    VECTOR<float> template_probabilities(n_templates + 1, 0.0);

    for (auto& line : lines)
    {
        StringList line_contents;

        boost::split(line_contents, line,
                     boost::algorithm::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        cerr << "template_id = " << line_contents[0] << " p = " << line_contents[1] << "\n";
        const uint template_id = std::stoi(line_contents[0]);
        const float probability = std::stod(line_contents[1]);

        template_probabilities[template_id] = probability;
    }

    return template_probabilities;
}


void RedshiftPDF::RedshiftProbDensity(const float redshift,
                                      SourceInfo& object,
                                      const ModelSeds& model_sed)
{
    // for each model in the model grid for this particular redshift, perform the chi^2 comparison

    const uint kn_models = model_sed.size() / (n_filters_ + kn_metadata_);  // number of model SEDs in this redshift bin

    double prob_density = 0;  // probability density for this object at this redshift

    double min_chisq = 9e100;

    int best_template = 0;

    for (uint model = 0; model < kn_models; ++model)
    {
        // fetch the address of the first entry of the current model's SED and treat it as
        // the first address in a regular C array.

        const float* mod_sed = &model_sed.data()[(n_filters_ + kn_metadata_) * model];

        // This entry needs to be interpreted as an int, but reinterpret_cast only works with pointers
        // and references, so we get the address of the float, reinterpret as the address of an int
        // then de-reference
        const int ktemplate_id = *reinterpret_cast<const int*>(&mod_sed[0]); // the first entry is the template id.

        if (ktemplate_id == 0 || ktemplate_id >= template_probability_.size())
        {
            cerr << "\nERROR: There is either an error in the model grid file OR the model grid \n"
                    "is incompatible with the template probability file that you provided. \n"
                    "template_id = " << ktemplate_id << " does not exist in the template probability "
                    "file.\n\n";
            exit(-1);
        }

        // if the probability of this template is zero, move on the the next model

        if (template_probability_[ktemplate_id] == 0.0f) { continue; }

        float total_weighted_template_flux = 0;
        float total_weighted_observed_flux = 0;

        for (uint i = 0; i < n_filters_; ++i)
        {
            const float kmodel_fluxi = mod_sed[kn_metadata_ + i];
            const SourceInfo::Flux& object_sed = object.sed[i];

            const float shared_product = object_sed.weight * kmodel_fluxi;

            total_weighted_observed_flux += shared_product * object_sed.value;

            total_weighted_template_flux += shared_product * kmodel_fluxi;
        }

        // compute the scale factor:
        // \frac{ \sum_i F_i * f_i * w_i } { \sum_i f_i^2 w_i }
        // where F_i is the measured flux associated with ith filter, f_i is the template flux
        // associated with the ith filter, and the weight, w_i = 1 / dF_i ^ 2, where dF_i is the
        // error in the measurement of the flux associated with the ith filter.

        const float inverse_template_weight = 1.0f / total_weighted_template_flux;

        const float kscale_factor = total_weighted_observed_flux * inverse_template_weight;

        // compute chi^2

        double chisq = 0;

        for (int i = 0; i < n_filters_; ++i)
        {
            const SourceInfo::Flux& object_sed = object.sed[i];

            const float diff = object_sed.value - mod_sed[kn_metadata_ + i] * kscale_factor;

            chisq += diff * diff * object_sed.weight;
        }

        // sum up the individual contributions, taking the template probability into account when
        // computing the probability.

        prob_density += template_probability_[ktemplate_id] * sqrt(inverse_template_weight) * exp(-0.5 * chisq);

        min_chisq = min(min_chisq, chisq);

        if (min_chisq == chisq)
        {
            best_template = ktemplate_id;
        }
    }

    // prob_density is divided by the number of contributions, kn_models, in order to account for
    // the possibility that redshift bins do not contain the same number of models.

    object.pdf.push_back( { redshift, prob_density / kn_models, min_chisq, best_template } );
}


void RedshiftPDF::ComputePDF(const RedshiftGrid &model_grid, SourceList &object_list)
{
    cerr << "\nn_filters = " << n_filters_ << "\n"
            "n_objects = " << object_list.size() << "\n"
            "n_redshift_bins = " << model_grid.size() << "\nBeginning to compute probability"
            " distributions\n";

    # pragma omp parallel for schedule(dynamic)

    for (int i = 0; i < object_list.size(); ++i) // parallel loop over all objects
    {
        for (const auto& redshift_bin : model_grid) // loop over redshift bins
        {
            // compute the contribution to the probability distribution for object at the given
            // redshift

            RedshiftProbDensity(redshift_bin.first,  // the redshift value.
                                object_list[i],      // all information about the current object.
                                redshift_bin.second);// flattened vector of model SEDs at this redshift.
        }
    }
}


void RedshiftPDF::SavePDFs(const char *pdf_filename) const
{
    ofstream pdfs(pdf_filename);

    if (pdfs.is_open())
    {
        // write a header for the file:

        pdfs << "# id xpos ypos redshifts: ";

        auto& zeroth_object = object_list_[0];

        // the redshift values

        for (auto& entry : zeroth_object.pdf) { pdfs << " " << entry.redshift; }

        pdfs << "\n";

        // write the PDFs

        for (const auto& object : object_list_ )
        {
            pdfs << object.id << " " << object.xpos << " " << object.ypos;

            for (auto& entry : object.pdf) { pdfs << " " << entry.density; }

            pdfs << "\n";
        }

        pdfs.close();
    }
    else { cerr << "Error: Could not open file, \"" << pdf_filename << "\".\n"; }
}


void RedshiftPDF::SaveRedshifts(const char *filename) const
{
    ofstream redshift_catalog(filename);

    redshift_catalog << setw(20) << "id"
                     << setw(20) << "ra"
                     << setw(20) << "dec"
                     << setw(20) << "z"
                     << setw(20) << "std-dev"
                     << setw(20) << "skewness"
                     << setw(20) << "kurtosis"
                     << setw(20) << "mlz"
                     << setw(20) << "min-chisq@mlz"
                     << setw(20) << "min-chisq-global"
                     << setw(20) << "mean-chisq-global"
                     << setw(20) << "best-template@mlz"
                     << setw(20) << "best-template"
                     << setw(20) << "mean-SNR\n";

    if (redshift_catalog.is_open())
    {
        for (const auto& object : object_list_)
        {
            const SourceInfo::ProbDist& pdf = object.pdf;

            // compute estimated redshift: z = z_i * dp_i / sum_i dp_i
            // and uncertainty: dz = sqrt( sum_i (dp_i * (z_i - z) ^ 2) / sum_i dp_i )

            double pdf_sum = 0; // sum_i dp_i
            double zdp_sum = 0; // z_i * dp_i
            double chisq_dp_sum = 0; // min_chisq_i * dp_i

            double min_chisq = 1e100; // global minimum chisq
            int best_template = 0;

            SourceInfo::ProbDensity most_probable = {0, 0, 0, 0};

            for (const SourceInfo::ProbDensity& redshift_density : pdf)
            {
                zdp_sum += redshift_density.redshift * redshift_density.density;
                pdf_sum += redshift_density.density;
                chisq_dp_sum += redshift_density.min_chisq;

                min_chisq = min(min_chisq, redshift_density.min_chisq);

                if (min_chisq == redshift_density.min_chisq)
                {
                    best_template = redshift_density.best_template;
                }

                if (redshift_density.density > most_probable.density)
                {
                    most_probable = redshift_density;
                }
            }

            double redshift_estimate, std_dev, skewness, kurtosis, mean_chisq;

            redshift_estimate = std_dev = skewness = kurtosis = mean_chisq -9.0;

            if (pdf_sum > 0)
            {
                const double inv_pdf_sum = 1.0 / pdf_sum;

                redshift_estimate = zdp_sum * inv_pdf_sum;
                mean_chisq = chisq_dp_sum * inv_pdf_sum;

                double second_moment = 0;
                double third_moment = 0;
                double fourth_moment = 0;

                for (const SourceInfo::ProbDensity& redshift_density : pdf)
                {
                    const float deviation = redshift_density.redshift - redshift_estimate;
                    const float deviation2 = deviation * deviation;

                    second_moment += deviation2 * redshift_density.density;
                    third_moment += deviation2 * deviation * redshift_density.density;
                    fourth_moment += deviation2 * deviation2 * redshift_density.density;
                }

                // standard deviation is the square root of the second moment (the variance)
                // sqrt (sum diff^2 * dp / sum dp)

                std_dev = sqrt(second_moment * inv_pdf_sum);

                if (std_dev > 0)
                {
                    const double inv_std_dev = 1.0 / std_dev;
                    const double inv_std_dev2 = inv_std_dev * inv_std_dev;

                    // skewness is the third moment divided by std_dev^3
                    // (sum diff^3 * dp / sum dp) / std_dev^3

                    skewness = third_moment * inv_pdf_sum * inv_std_dev2 * inv_std_dev;

                    // kurtosis is the fourth moment divided by std_dev^4
                    // (sum diff^4) * dp / sum dp / std_dev^4

                    kurtosis = fourth_moment * inv_pdf_sum * inv_std_dev2 * inv_std_dev2;
                }

                // make sure that the outputs are valid numbers. If one is not finite, then
                // set all values to -9 (filter out NaNs and Infs)

                bool valid_output = (std::isfinite(std_dev) &&
                                     std::isfinite(redshift_estimate) &&
                                     std::isfinite(kurtosis) &&
                                     std::isfinite(skewness));

                if ( !valid_output )
                {
                    redshift_estimate = std_dev = skewness = kurtosis = -9.0;
                }
            }
            else { pdf_sum = -9.0; }

            if (!std::isfinite(pdf_sum)) { pdf_sum = -9.0; }

            redshift_catalog << setw(20) << object.id
                             //<< std::setprecision(8) << std::setiosflags(ios::scientific | ios::internal)
                             << std::setprecision(8)
                             << setw(20) << object.xpos
                             << setw(20) << object.ypos
                             << setw(20) << redshift_estimate
                             << setw(20) << std_dev
                             << setw(20) << skewness
                             << setw(20) << kurtosis
                             << setw(20) << most_probable.redshift
                             << setw(20) << most_probable.min_chisq
                             << setw(20) << min_chisq
                             << setw(20) << mean_chisq
                             << setw(20) << most_probable.best_template
                             << setw(20) << best_template
                             << setw(20) << object.snr / n_filters_ << "\n";
        }

        redshift_catalog.close();
    }
    else { cerr << "could not open file, \"" << filename << "\"\n"; }
}

