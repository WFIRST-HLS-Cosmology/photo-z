/*! \file
 *
 * \brief Contains declarations of the SED manipulation functions.
*/

#ifndef SED_H
#define SED_H

typedef double real_t;

/*!
 * @brief allen_k implements the extinction curve from allen (MW 1976)
 * @param wave observed wavelengths (A)
 * @param k relevant attenuation
 * @param Np  number of bands
 */
void allen_k(real_t* wave, real_t* k, int Np);

/*!
 * @brief prevot_k implements the extinction curve from prevot et al (SMC 1984)
 * @param wave observed wavelengths (A)
 * @param k relevant attenuation
 * @param Np number of bands
 */
void prevot_k(real_t* wave, real_t* k, int Np);

/*!
 * @brief calz_k implements the extinction curve from calzetti et al (SB 2000)
 * @param wave observed wavelengths (A)
 * @param k relevant attenuation
 * @param Np number of bands
 */
void calz_k(real_t* wave, real_t* k, int Np);

/*!
 * @brief seaton_k implements the extinction curve from seaton (MW 1979)
 * @param wave observed wavelengths (A)
 * @param k relevant attenuation
 * @param Np number of bands
 */
void seaton_k(real_t* wave, real_t* k, int Np);

/*!
 * @brief fitzpatrick_k implements the extinction curve from fitzpatrick (LMC 1986)
 * @param wave observed wavelengths (A)
 * @param k relevant attenuation
 * @param Np number of bands
 */
void fitzpatrick_k(real_t* wave, real_t* k, int Np);

/*!
 * @brief make_dustvec GENERATES A DUST VECTOR FOR A GIVEN EXTINCTION LAW
 * @param wave wavelengths (A)
 * @param flux fluxes (Jy)
 * @param law reddening law
 * @param Np number of points (bands)
 */
void make_dustvec(real_t* wave, real_t* flux, int law, int Np);

/*!
 * @brief madau_teff1 COMPUTES OPTICAL DEPTH MODIFICATION DUE TO H-CLOUDS ON LAMBDA=912-1216
 * (LYMAN-ALPHA AND LYMAN LIMIT) PHOTONS (Madau 1995).
 * (i.e., corrections from foreground H-clouds, which occurs for the Lyman series)
 * @param wavelength observed wavelength (A)
 * @param z redshift
 * @return effective optical depth
 */
const real_t inline madau_teff1(real_t wavelength, real_t z);

/*!
 * @brief madau_teff2 COMPUTES OPTICAL DEPTH MODIFICATION DUE TO H-CLOUDS ON LAMBDA<912
 * (BLUEWARD OF LYMAN LIMIT) PHOTONS (Madau 1995)
 * @param x observed wavelength (A)
 * @param z redshift
 * @return effective optical depth (converted from optical depth)
 */
const real_t inline madau_teff2(real_t x, real_t z);

/*!
 * @brief madau_teff APPLIES THE EFFECTIVE TRANSMISSION FROM INTERVENING H-CLOUDS USING MADAU'S
 * CORRECTIONS FOR A GIVEN WAVEVECTOR
 * @param x observed wavelength (A)
 * @param z redshift
 * @return effective transmission
 */
const real_t inline madau_teff(real_t x, real_t z);

/*!
 * @brief apply_madau APPLIES THE EFFECTIVE TRANSMISSION FROM INTERVENING H-CLOUDS USING MADAU'S
 * CORRECTIONS FOR A GIVEN WAVEVECTOR
 * @param wave restframe wavelengths (A)
 * @param flux fluxes (Jy)
 * @param z observed redshift
 * @param Np number of points
 */
void inline apply_madau(const real_t* wave,
                        real_t* flux,
                        const real_t z,
                        const int Np) __attribute__ ((hot));

/*!
 * @brief waveidx COMPUTES WAVELENGTH INDICES---THE NUMBER OF WAVELENGTH BINS LOCATED BETWEEN TWO LIMITS,
 GIVEN A SET LOG SPACING
 * @param lmin minimum wavelength
 * @param dlogwave log spacing
 * @param l wavelength
 * @return THE NUMBER OF WAVELENGTH BINS
 */
int waveidx(real_t lmin, real_t dlogwave, real_t l);


/// SETS THE OPTIMAL SPACING TO BE USED FOR SAMPLING, IN PREPARATION FOR RESAMPLE_SEDS
/// note that this just initializes everything - the resample_seds actually does all the work
// inputs are SED wavelengths (double array; A),
// number of points (per SED),
// number of SEDs, size of redshift step, max resolution of filter (A),
// number of wavelengths that should be sampled, minimum wavelength (A), log wavelength spacing
void set_sampling(real_t** Swave,
                  int* Ns,
                  int Nsed,
                  real_t dz,
                  real_t Rmax,
                  int* outNwave,
                  real_t* outlmin,
                  real_t* outdlogwave);


/// TAKES INPUT SEDS (AT ARBITRARY WAVELENGTHS) AND TURNS THEM INTO SEDS SAMPLED AT THE DESIRED
// WAVELENGTHS
// inputs are SED wavelengths (double array; A), SED fluxes (double array; Jy), number of points
// (per SED), number of SEDs, output SED wavelengths, output frequencies (Hz),
// output frequency spacing (Hz), output fluxes (Jy), number of total wavelengths sampled,
// minimum wavelength, log wavelength spacing
void resample_seds(real_t** Swave,
                   real_t** Sflux,
                   int* Ns,
                   int Nsed,
                   real_t* outSwave,
                   real_t* outnu,
                   real_t* outdnu,
                   real_t** outSflux,
                   int Nwave,
                   real_t lmin,
                   real_t dlogwave);


/// PUTS FILTERS ON LOG GRID (PREPARATION FOR RESAMPLE_FILTERS)
// inputs are output SED wavelengths (A), minimum wavelength, log wavelength spacing,
// filter wavelengths (double array; A), filter transmission (double array),
// number of filters, number of wavelengths per filter, number of wavelengths down to start
// sampling the filter from, number of elements/post-sampled wavelengths in filter
void setup_resample_filters(real_t* outSwave,
                            real_t lmin,
                            real_t dlogwave,
                            real_t** Fwave,
                            real_t** Ftrans,
                            int NF,
                            int* NFwave,
                            int* Fstart,
                            int* F_nelem);


/// MAKES REDDENING VECTORS THAT CAN BE SUBSEQUENTLY APPLIED TO SEDS
// inputs are output SED wavelengths (A), output reddening vectors (double array; k), number
// of rlaws, and number of desired wavelengths to be sampled
void make_reddening_vec(real_t* outSwave,
                        real_t** outRvec,
                        int Nrlaw,
                        int Nwave);


/// TAKES INPUT FILTERS (AT ARBITRARY WAVELENGTHS) AND TURNS THEM INTO FILTERS SAMPLED AT THE
// DESIRED WAVELENGTHS (IDENTICAL TO SEDS)
// inputs are output SED wavelengths (A), output filter wavelengths (double array; A),
// output filter transmission (double array), number of filters, number of wavelengths per filter,
// number of wavelengths down to start sampling each filter from, output number of elements
// per filter, output filter transmissions (double array)
void resample_filters(real_t* outSwave,
                      real_t** Fwave,
                      real_t** Ftrans,
                      int NF,
                      int* NFwave,
                      int* Fstart,
                      int* outFnelem,
                      real_t** outFtrans);


/// QUICKLY INTEGRATES THE SED (ALONG WITH THE ASSOCIATED REDDENING VECTOR AND EXTINCTION DUE TO
// H-CLOUDS) OVER THE ASSOCIATED FILTERS TO GIVE OUTPUT PHOTOMETRY
// inputs are log wavelength spacing, SED wavelengths (A), frequencies (Hz),
// frequency spacings (Hz), SED fluxes (Jy), number of points, number of elements per filter,
// starting index per filter, filter transmissions, number of filters,
// the associated reddening vector, redshift, E(B-V) value, output fluxes (Jy),
// output errors in flux (Jy)
void integrate_sed(real_t dlogwave, // log wavelength spacing
                   const real_t* Swave, // SED wavelengths (A)
                   const float* dnu_nu, // frequency spacings / frequencies
                   real_t* Sflux, // SED fluxes (Jy)
                   int Ns,  // number of points (in what? the SED?)
                   const int* Fnelem,   // number of elements per filter
                   const int* Fstart,   // starting index per filter
                   real_t** Ftrans, // filter transmissions
                   int NF,        // number of filters
                   const real_t* Rvec,  // reddening vector
                   real_t z,
                   real_t Ebv,
                   real_t* cache,
                   real_t* FluxOut); // output fluxes (Jy)


/// Creates a redshifted, reddened SED, starting with a de-reddened, rest-frame SED template.
void ComputeSed(const real_t* Swave, // SED wavelengths
                real_t* Sflux,   // inout SED fluxes
                const int Ns,    // number of flux-wavelength entries in the template
                const real_t* Rvec, // reddening vector
                real_t z,        // redshifts
                real_t Ebv,      // E(B-V)
                real_t* fluxes); // template fluxes

#endif // SED_H

