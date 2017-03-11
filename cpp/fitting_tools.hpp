/*! \file
 *
 * Contains classes used for the actual redshfift fitting routine.
*/

#ifndef FITTING_TOOLS_H
#define FITTING_TOOLS_H

#include <glob.h>

#ifdef HAS_FOLLY
    #include <folly/FBVector.h>
    #define VECTOR folly::fbvector
#else
    #include <vector>
    #define VECTOR std::vector
#endif

#include <map>
#include <iostream>

/*!
 * @brief ModelSeds is used for storing fluxes of the model SEDs.
 */
typedef VECTOR<float> ModelSeds;

/*!
 * @brief RedshiftGrid stores the seds of each redshift bin in an efficient manner.
 */
typedef std::map<float, ModelSeds> RedshiftGrid;

/*!
 * @brief The SourceInfo struct stores all information about an observed source / SED.
 */
struct SourceInfo
{
    /*!
    * \brief The Flux struct stores flux data for SourceInfo objects.
    */
   struct Flux
   {
      float value;
      float weight;
   };

   /*!
    * \brief The ProbDensity struct stores redshift probability density information.
    */
   struct ProbDensity
   {
       float redshift;
       double density;
       double min_chisq;
       int best_template;
   };

   typedef VECTOR<Flux> FluxVector;
   typedef VECTOR<ProbDensity> ProbDist;

   uint id;
   float xpos;
   float ypos;
   float weighted_flux = 0; // sum of flux_i * weight_i
   float total_flux = 0; // sum of flux_i
   float total_squared_flux = 0; // sum of flux_i * flux_i
   float total_weighted_squared_flux = 0; // sum of flux_i * flux_i * weight_i
   float total_weight = 0; // sum of weight_i
   float total_noise = 0; // sum of \sigma_i^2
   float snr = 0;
   FluxVector sed;
   ProbDist pdf;
};

/*!
 * @brief SourceList stores a collection of SourceInfo objects (i.e., information for a collection
 * of observed sources).
 */
typedef VECTOR<SourceInfo> SourceList;

/*!
 * @brief The RedshiftPDF class does all of the work needed to compute redshift probability
 * distribution functions and moments of the distributions.
 */
class RedshiftPDF final
{
public:

    /*!
     * \brief RedshiftPDF constructor reads the model grid file(s), the file containing observed SEDs,
     * and the template probability file. When all data is loaded and parsed, it begins the process
     * of computing the redshift probability distributions (PDFs).
     *
     * \param model_grid_filename is the filename (or filename pattern, containing wildcards), specifying
     * the model grid file(s), produced by the make-grid program.
     *
     * \param object_sed_filename is the name of a text file containing the observed fluxes of an object,
     * with one object per line. The order of the filter fluxes in the file must match the order of the
     * filters, used to generate the model grid.
     *
     * \param template_probability_filename is a two-column text file with the first column representing
     * the ID number of a template and the second column representing the probability that the particular
     * template will match an observed galaxy.
     */
    RedshiftPDF(const char* model_grid_filename,
                const char* object_sed_filename,
                const char* template_probability_filename)
        : kmodel_grid_file_(model_grid_filename),
          kobject_sed_file_(object_sed_filename)

    {
        uint filters_in_grid = -1;
        uint filters_in_seds = -1;

        // if multiple threads are available and the code is compiled with OpenMP support, the following
        // two sections will execute simultaneously (i.e., the model grid and SEDs will be read and parsed
        // simultaneosly.)

        #pragma omp parallel sections
        {
             #pragma omp section
             {
                  glob_t glob_results;

                  int glob_status = glob(kmodel_grid_file_, GLOB_ERR, nullptr, &glob_results);

                  if (glob_status != 0)
                  {
                      std::cerr << "glob() failed to find files matching pattern, " << kmodel_grid_file_
                           << ". Exited with code: " << glob_status << "\n";

                      exit(glob_status);
                  }

                  for (uint part = 0; part < glob_results.gl_pathc; ++part)
                  {
                      filters_in_grid = ReadModelGrid(glob_results.gl_pathv[part]);
                  }

                  globfree(&glob_results);
             }

             #pragma omp section
             {
                   filters_in_seds = ReadSourceSeds(object_sed_filename);
             }
         }

        // The number of filters (flux bins) present in the model grid must math the number present
        // in the observed SED list.

        if (filters_in_grid != filters_in_seds || (filters_in_grid + filters_in_seds) < 2 )
        {
            std::cerr << "Error: filters_in_grid = " << filters_in_grid << ", "
                         "filters_in_seds = " << filters_in_seds << "\n";
            exit(11);
        }

        n_filters_ = filters_in_grid;

        template_probability_ = ReadTemplateProbabilities(template_probability_filename);

        ComputePDF(model_grid_, object_list_);
    }

    ~RedshiftPDF() {}

    VECTOR<float> ReadTemplateProbabilities(const char* filename);

    void SavePDFs(const char* pdf_filename) const;

    void SaveRedshifts(const char* filename) const;

private:

    /// Reads the model grid from model_grid_filename and populates the model_grid data structure.
    /// The model grid is stored in a custom (somewhat fragile) binary format, defined in the source
    /// for the make-grid executable (see make_grid.cpp).
    /// Returns the apparent number of filters in the model grid file.
    uint ReadModelGrid(const char* model_grid_file_name);

    /// Reads the source SEDs from object_sed_filename and populates the object_list data structure.
    /// Returns the apparent number of filters in the model grid file.
    uint ReadSourceSeds(const char* object_sed_filename);

    /*! Compares all of the model SEDs at a given redshift with a single object SED and updates
        the object's redshift PDF by adding the entry for this particular redshift. In order to
        compare the measured fluxes with the templates, the templates are first scaled to match the
        imput data by multiplying it by a scale factor, \f$ N \f$, such that

        \f[
             f_i' = N f_i
        \f]

        where \f$f_i'\f$ is the scaled template flux associated with the \f$i\f$th filter and \f$f_i\f$
        is the original template flux associated with the \f$i\f$th filter.

        The scale factor, \f$N\f$, is computed as

        \f[
            N = \frac{ \sum_i F_i f_i  w_i } { \sum_i f_i^2 w_i }
        \f]

        where \f$F_i\f$ is the measured flux associated with the \f$i\f$th filter and the weight, \f$w_i\f$,
        associated with \f$F_i\f$ is defined as

        \f[
            w_i \equiv \frac{1}{dF_i}
        \f]

        where \f$dF_i\f$ is the error in the measurement of \f$F_i\f$.

        Once the scaling is performed, \f$\chi^2\f$ is computed:

        \f[
             \chi^2 = \sum_i (F_i - f_i')^2 w_i
        \f]

        Finally, \f$P(z)\f$ is computed by summing over all \f$\ch^2\f$ values for the redshift bin:

        \f[
            P(z) = \sum_i P(T_i) e^{-\chi_i^2 / 2}
        \f]

        Where \f$P(T_i)\f$ is the probability that the template associated with \f$\chi_i^2\f$ is the
        correct template. Note that $P(T_i)$ is provided as an input file, template_probability_filename.
        The probabilities are computed by matching the template library to the COSMOS photometry.
        For consistency, only objects that are detectable at this stage of the pipeline (i.e., using
        the detection limits of PanSTARRS and WISE) are included in the template probability determination.
     **/
    void RedshiftProbDensity(const float redshift,
                             SourceInfo& object,
                             const ModelSeds& model_sed) __attribute__ ((hot));

    /// Computes the redshift probability distribution function for each object in object_list
    /// by comparing each object's SED to each model SED in the model_grid;
    void ComputePDF(const RedshiftGrid& model_grid, SourceList& object_list);

    const char* kmodel_grid_file_;
    const char* kobject_sed_file_;
    const uint kn_metadata_ = 1; // number of 4-byte words needed to store model metadata 
    VECTOR<float> template_probability_;
    uint n_filters_;
    RedshiftGrid model_grid_;
    SourceList object_list_;
};

#endif // FITTING_TOOLS_H
