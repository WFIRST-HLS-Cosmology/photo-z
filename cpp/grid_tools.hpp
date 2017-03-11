/*! \file
 *
 * The primary header for the simulated photometry / model grid generation utilities.
 */

#ifndef GRID_TOOLS_H
#define GRID_TOOLS_H

#include "sed.hpp"
#include "basic_utils.hpp"

#ifdef HAS_FOLLY
    #include <folly/FBVector.h>
    #define VECTOR folly::fbvector
#else
    #include <vector>
    #define VECTOR std::vector
#endif

#include <string>
#include <iostream>
#include <iomanip>



/*! The catalog_column namespace serves a purpose similar to that of an enum class. Specifically,
 * the code becomes more readable when the enum is prefixed with a descriptive name. This allows us
 * to use the enumeration values as array indices, which isn't possible with enum classes unless we
 * static_cast<uint> the enumerator.
 */
namespace catalog_column
{

/*! The raw_catalog enum specifies the column indices for the metadata columns of the
 * input photometry catalog. The remaining columns of the catalog are assumed to be fluxes. Note
 * that the column indices in this catalog are different from those of the grid specification catalog.
 */
enum raw_catalog : uint
{
    id = 0,
    logmass = 1,
    nbands = 2,
    RA = 3,
    DEC = 4,
    I_magnitude = 5,
    redshift = 6
};
}

/*! The gridspec_column namespace wraps the grid_specification column enumeration to make the column
 * indices more readable.
 *
 */
namespace gridspec_column
{
enum gridspec_file : uint
{
    id = 0,
    RA = 1,
    DEC = 2,
    I_magnitude = 3,
    redshift = 4,
    template_id = 5,
    reddening = 6,
    reddening_law = 7,
    normalization = 8,
    line_scale = 9
};
}



/*! The ChiSquaredResult struct is a convenience struct for storing the results of the chi squared
 * comparison operation.
 */
struct ChiSquaredResult
{
    float value;
    float scale_factor;
};

/*!
 * \brief The Flux struct stores flux data (a flux-weight pair).
 */
struct Flux
{
    float value;
    float weight;
};

//

/*!
 * \brief The OptimalParams struct is a convenience struct for storing the optimal parameters of the
 * fit.
 */
struct OptimalParams
{
    ChiSquaredResult chi2;
    float ebv;
    int rlaw;
};

/*!
 * @brief The SliceLimits class parses an output filename and parses sice information, if present.
 */
class SliceLimits
{
public:

    /*!
     * @brief Computes the indices of the slice specified in the filename.
     * @param filename is the name of the output file. If the name ends with .part_NUM1_of_NUM2, the
     * NUM1 is interpreted as a numerator and NUM2 is interpreted as a denominator (a slice number
     * and the total number of slices).
     * @param array_size is the total size of the array that will be sliced into NUM2 slices.
     *
     * When supplied with a filename that contains slice information, in the format
     * .part_NUM1_of_NUM2, this class comptutes the indices in an array of length array_size that
     * bound the slice specified in the filename. The slice is [ first(), last() ).
     */
    SliceLimits(const char* filename, uint array_size)
    {
        std::string name(filename);

        size_t index_of_dotpart = name.rfind(".part_");

        if(std::string::npos == index_of_dotpart) // name does not contain .part_
        {
            first_ = 0;
            last_ = array_size;
        }
        else
        {
            uint index_of_fraction = index_of_dotpart + 6; // 6 = length of ".part_"

            std::string fraction = name.substr(index_of_fraction);

            uint sep_index = fraction.find("_of_");

            uint part = std::stoi(fraction.substr(0, sep_index));

            uint n_slices = std::stoi(fraction.substr(sep_index + 4)); // 4 = length of "_of_"

            if (part ==  0 || part > n_slices || n_slices == 0)
            {
                std::cerr << "ERROR: " << name << " is an invalid slice name\n";
                exit(2);
            }

            first_ = (part - 1) * array_size / n_slices;

            last_ = (part < n_slices) ? part * array_size / n_slices : array_size;
        }
    }

    /*!
     * @brief first returns the first index in the array belonging to the specified slice.
     */
    const uint first() const { return first_; }

    /*!
     * @brief last returns the index in the array, just beyond the specified slice.
     */
    const uint last() const { return last_; }

private:

    uint first_;
    uint last_;
};


/*! Stores the model grid data from fluxes[][] in a binary file format that can be quickly read by
    the photo-z executable. The first four bytes of a photoz file are "nflt" in ASCII. This is followed
    by four bytes that should be interpreted as an unsigned integer containing the number of filters
    used to construct the SEDs. The third group of four bytes contains "zbin". The symbol "zbin" is
    always immediately followed by a single precision floating point number (4 bytes) and then a list
    of SEDs, which each contain the same number of flux values, corresponding to the nflt filters.
    Subsequent occurances of "zbin" indicate that the following seds are in a new redshift bin.

    No attempt is made to check the byte order because it is assumed that the same machine that stored
    the photoz file will later read it; there should be no movement between Little Endian and Big Endian
    machines during the course of one simulation. */
void SavePhotoZBinary(const char* filename,
                      const uint n_filters,
                      const uint n_models,
                      double** metadata,
                      double** fluxes);


/*!
 * \brief The ResampledData struct stores all of the output of the ResampleFiltersAndTemplates function.
 */
struct ResampledData
{
    int n_wave;       // should be implicit because it's the length of sed_integral_wave
    int n_filters;    // should be implicit because it's the length of filter_transmissions
    int n_template;
    int n_rlaw;
    int n_ebv;
    double dlogwave;
    VECTOR<int> filter_n;    // should be implicit because the filter_transmissions should be returned as a vector
    VECTOR<int> filter_starting_index;
    VECTOR<float> dnu_nu; // really should be combined with the wavelenghts, below
    double* sed_integral_wave;     // SED wavelengths (A)
    double** sed_integral_flux;
    double** filter_transmissions; // filter_integral_flux
    double** red_Data;
};

/*!
 * @brief ChiSquared measures the difference between two vectors of photometry measurements
 * @param observed_fluxes contains the integrated fluxes from observations
 * @param template_fluxes contains the integrated fluxes from a template SED
 * @return returns a measure of the difference between observed_fluxes and template_fluxes
 */
inline const ChiSquaredResult ChiSquared(const VECTOR<Flux>& observed_fluxes,
                                         const float total_weighted_observed_flux,
                                         const VECTOR<float>& template_fluxes);


void MakePhotometryFittingGrid(const char* catalog_filename,
                               const char* filterlist_filename,
                               const char* templatelist_filename,
                               const char *output_filename);

/// this is not implemented yet
void MakeRedshiftFittingGrid(const char* gridspec_filename,
                             const char* filterlist_filename,
                             const char* templatelist_filename,
                             const char* output_grid_filename);

/*!
 * @brief FindBestFits Finds the combination of reddening, reddening law, and template that best
 * fits the integrated filter fluxes, given in a catalog.
 * @param catalog_filename is the name of a text file containing flux information for a collection
 * with one observed source on each line of the file.
 * @param model_grid_filename is the name of the model grid that was pre-computed using
 * MakePhotometryFittingGrid()
 * @param output_filename is the name of the output file, containing the best-fitting parameters for
 * each object in the catalog.
 */
void FindBestFits(const char* catalog_filename,
                  const char* model_grid_filename,
                  const char* output_filename);

/*!
 * @brief UniqueRedshiftsInCatalog reads a catalog file and returns a list of unique redshifts.
 * @param catalog_file is the name of the catalog of interest
 * @return A list of unique redshifts, sorted from low to high. Redshifts are considered unique if
 * they differ by more than 0.0005.
 */
VECTOR<float> UniqueRedshiftsInCatalog(const char* catalog_file);

/*!
 * @brief RegularizeFiltersAndTemplates Interpolates and resamples template SEDs and filters so that
 * the filters can be integrated against the templates
 * @param filterlist
 * @param templatelist
 * @param output
 */
void ResampleFiltersAndTemplates(const char* filterlist,
                                 const char* templatelist,
                                 ResampledData& output,
                                 const uint Rmax = 150,
                                 const float dzmin = 0.0023);

/*!
 * @brief WriteTemplates saves the wavelenths and fluxes of a normalized template to a file, rather
 * than integrating the template against filters.
 * @param output_file is a file containing the wavelengths and fluxes for each object listed in the
 * input catalog
 * @param resampled_templates is the output of ResampleFiltersAndTemplates
 * @param n_entries is the number of entries in the input catalog
 * @param catalog_data is the data in the input catalog, already converted to real_t.
 */
void WriteTemplates(const std::string output_file,
                    const ResampledData& resampled_templates,
                    const int n_entries,
                    real_t** catalog_data);

/*!
 * @brief ScaleEmissionLines Modifies regions of the template SED by multiplying emission lines by
 * a constant factor.
 * @param n_wavelengths the number of samples in the template (wavelength-flux pairs)
 * @param fluxes the fluxes of each sample, in Jy.
 * @param wavelengths is an array containint the wavelengths of each sample (in Angstroms).
 * @param line_correction is the factor by which the fluxes are multiplied
 * @return An array containing the modified template
 *
 * Multiplies the following regions by line_correction.
 *      [OII]3727:  [3715, 3742]
 *      H-beta:     [4850, 4873]
 *      [OIII]4959: [4950, 4970]
 *      [OIII]5007: [4995, 5020]
 *      H-alpha/[NII]: [6540, 6595]
 */
double* ScaleEmissionLines(const uint n_wavelengths,
                           const double* fluxes,
                           const double* wavelengths,
                           const float line_correction);

#endif // GRID_TOOLS_H
