/*! \file

  Contains definitions of the SED manipulation functions.
 */

#include "sed.hpp"
#include "fitting_tools.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <omp.h>

using namespace std;


#define MAXLINE 1024 // maximum length of line

// INITIALIZING FUNCTIONS THAT MODIFY DATA BASED ON DIFFERENT EXTINCTION LAWS
// FOR GIVEN INPUT PARAMETERS


void allen_k(real_t* wave, real_t* k, int Np)
{
    real_t Rv = 3.1; // slope of the dust extinction law, = Av/E(B-V) or total/selective extinction

    real_t  wa[20] = {1000.0,
                      1110.0,
                      1250.0,
                      1430.0,
                      1670.0,
                      2000.0,
                      2220.0,
                      2500.0,
                      2850.0,
                      3330.0,
                      3650.0,
                      4000.0,
                      4400.0,
                      5000.0,
                      5530.0,
                      6700.0,
                      9000.0,
                      10000.0,
                      20000.0,
                      100000.0
                     }; // discrete wavelength ranges where the data has been observed

    real_t ka[20] = {4.20,
                     3.70,
                     3.30,
                     3.00,
                     2.70,
                     2.80,
                     2.90,
                     2.30,
                     1.97,
                     1.69,
                     1.58,
                     1.45,
                     1.32,
                     1.13,
                     1.00,
                     0.74,
                     0.46,
                     0.38,
                     0.11,
                     0.00
                    }; // corresponding attenuations at the given wavelength

    // the gsl_spline group of commands allows GSL to do a cubic spline interpolation of the data
    gsl_spline* interpolation = gsl_spline_alloc (gsl_interp_cspline, 20);
    gsl_spline_init(interpolation, wa, ka, 20);
    gsl_interp_accel* accelerator =  gsl_interp_accel_alloc();

    for (int i = 0; i < Np; ++i)
    {
        if ((wave[i] >= wa[0]) && (wave[i] <= wa[19]))
        {
            // assign values via interpolation and shift based on Rv
            k[i] = gsl_spline_eval(interpolation, wave[i], accelerator) + Rv;
        }
        else
        {
            if (wave[i] < wa[0])
            {
                // bluer wavelengths get attenuated by constant * delta_lambda (dex) + offset
                k[i] = 4.791107 * ( log(wa[0]) - log(wave[i]) ) + ka[0] + Rv;
            }
            else
            {
                k[i] = 0; // redder wavelengths set to 0
            }
        }
    }

    gsl_spline_free(interpolation);
    gsl_interp_accel_free(accelerator);
}


void prevot_k(real_t* wave, real_t* k, int Np)
{
    real_t Rv = 2.72;

    double wa[39] = {1275,
                     1330,
                     1385,
                     1435,
                     1490,
                     1545,
                     1595,
                     1647,
                     1700,
                     1755,
                     1810,
                     1860,
                     1910,
                     2000,
                     2115,
                     2220,
                     2335,
                     2445,
                     2550,
                     2665,
                     2778,
                     2890,
                     2995,
                     3105,
                     3704,
                     4255,
                     5291,
                     12500,
                     16500,
                     22000,
                     24000,
                     26000,
                     28000,
                     30000,
                     32000,
                     34000,
                     36000,
                     38000,
                     40000
                    };

    double ee[39] = {13.54,
                     12.52,
                     11.51,
                     10.8,
                     9.84,
                     9.28,
                     9.06,
                     8.49,
                     8.01,
                     7.71,
                     7.17,
                     6.9,
                     6.76,
                     6.38,
                     5.85,
                     5.3,
                     4.53,
                     4.24,
                     3.91,
                     3.49,
                     3.15,
                     3.0,
                     2.65,
                     2.29,
                     1.81,
                     1.0,
                     0.0,
                     -2.02,
                     -2.36,
                     -2.47,
                     -2.51,
                     -2.55,
                     -2.59,
                     -2.63,
                     -2.67,
                     -2.71,
                     -2.75,
                     -2.79,
                     -2.83
                    };

    gsl_spline* interpolation = gsl_spline_alloc (gsl_interp_cspline, 39);
    gsl_spline_init(interpolation, wa, ee, 39);
    gsl_interp_accel* accelerator =  gsl_interp_accel_alloc();

    for (int i = 0; i < Np; ++i)
    {
        if ((wave[i] <= 34500) && (wave[i] > 1275))
        {
            // assign values and moving by offsets
            k[i] = gsl_spline_eval(interpolation, wave[i], accelerator) + Rv;
        }
        else
        {
            if (wave[i] > 34500)
            {
                k[i] = 0; // no redward extinction
            }
            else
            {
                // blueward values attenuated by constant * delta_lambda (dex) + offset + Rv shift
                k[i] = 24.151865 * (log(wa[0]) - log(wave[i])) + ee[0] + Rv;
            }
        }
    }

    gsl_spline_free(interpolation);
    gsl_interp_accel_free(accelerator);
}


void calz_k(real_t* wave, real_t* k, int Np)
{
    real_t Rv = 4.05, w, iwave;

    for (int i = 0; i < Np; ++i)
    {
        w = wave[i] * 1.0e-4; // convert to microns
        iwave = 1.0 / w; // invert

        if (w < 0.63)
        {
            k[i] = (2.659 * (-2.156 + (1.509 * iwave) - (0.198 * pow(iwave, 2)) + (0.011 * pow(iwave, 3)))) + Rv; // redward portion, cubic fit w/ offset
        }
        else
        {
            k[i] = (2.659 * (-1.857 + (1.040 * iwave))) + Rv; // blueward portion, linear fit w/ offset
        }

        if (k[i] < 0)
        {
            k[i] = 0; // setting negative extinction to 0
        }
    }
}


void seaton_k(real_t* wave, real_t* k, int Np)
{
    real_t Rv = 3.1, w, iwave;

    // declaring useful constants
    real_t lo = 4.595; // lorentzian location parameter
    real_t gamma = 1.051; // lorentzian scale parameter
    real_t c1 = -0.38;
    real_t c2 = 0.74;
    real_t c3 = 3.96;
    real_t c4 = 0.26;

    real_t* ak = (real_t*) malloc(sizeof(real_t) * Np); // allocates memory
    allen_k(wave, ak, Np); // applies allen_k rlaw for given input parameters

    for (int i = 0; i < Np; ++i)
    {
        w = wave[i] * 1.0e-4; // convert to microns
        iwave = 1.0 / w; // invert

        if (iwave >= 5.9)
        {
            k[i] = c1 + c2 * iwave + c3 * pow(iwave, 2) / (pow((pow(iwave, 2) - pow(lo, 2)), 2)
                                                           + pow(gamma, 2) * pow(iwave, 2))
                    + c4 * (0.539 * pow((iwave - 5.9), 2) + 0.0564 * pow((iwave - 5.9), 3)); // fit to bluer data
            k[i] += Rv; // applying Rv offset
        }
        else
        {
            if ((iwave < 5.9) && (iwave >= 2.74))
            {
                k[i] = c1 + c2 * iwave + c3 * pow(iwave, 2) / (pow((pow(iwave, 2) - pow(lo, 2)), 2) + pow(gamma, 2) * pow(iwave, 2)); // fit to redder data
                k[i] += Rv; // applying Rv offset
            }
            else
            {
                k[i] = ak[i]; // for redder data, keeps allen_k values
            }
        }
    }

    free(ak);
}


void fitzpatrick_k(real_t* wave, real_t* k, int Np)
{
    real_t Rv = 3.1, w, iwave;

    real_t lo = 4.608;
    real_t gamma = 0.994;
    real_t c1 = -0.69;
    real_t c2 = 0.89;
    real_t c3 = 2.55;
    real_t c4 = 0.50;

    real_t* ak = (real_t*) malloc(sizeof(real_t) * Np);
    allen_k(wave, ak, Np); // applying allen_k rlaw for given input parameters

    for (int i = 0; i < Np; ++i)
    {
        w = wave[i] * 1.0e-4; // convert to microns
        iwave = 1.0 / w; // invert

        if (iwave >= 5.9)
        {

            k[i] = c1 + c2 * iwave + c3 / (pow((iwave - pow(lo, 2) * w), 2) + pow(gamma, 2)) + c4 * (0.539 * (pow((iwave - 5.9), 2)) + 0.0564 * (pow((iwave - 5.9), 3))); // fit to bluer data
            k[i] += Rv; // applying Rv offset
        }
        else
        {
            if ((iwave < 5.9) && (iwave >= 3.0))
            {
                k[i] = c1 + c2 * iwave + c3 / (pow((iwave - pow(lo, 2) * w), 2) + pow(gamma, 2)); // fit to redder data
                k[i] += Rv; // applying Rv offset
            }
            else
            {
                k[i] = ak[i]; // for redder data, keep allen_k values
            }
        }
    }

    free(ak);
}


void make_dustvec(real_t* wave, real_t* flux, int law, int Np)
{
    real_t* k = (real_t*) malloc(sizeof(real_t) * Np);

    switch (law)
    {
    case 1:
        {
            prevot_k(wave, k, Np);
            break;
        }

    case 2:
        {
            calz_k(wave, k, Np);
            break;
        }

    case 3:
        {
            seaton_k(wave, k, Np);
            break;
        }

    case 4:
        {
            allen_k(wave, k, Np);
            break;
        }

    case 5:
        {
            fitzpatrick_k(wave, k, Np);
            break;
        }

    default:
        {
            cerr << "an unknown reddening law has been specified\n";
            exit(1);
        }
    }

    for (int i = 0; i < Np; ++i)
    {
        flux[i] = k[i];
    }

    free(k);
}


const real_t inline madau_teff1(real_t wavelength, real_t z)
{
    real_t teff1 = 0.0;
    const real_t z_plus_one = 1.0 + z;
    const real_t x_power = pow(wavelength,3.46);

    // lyman alpha line (2->1)
    constexpr const float lambda1 = 1216.0; // declaring relevant restframe wavelength, starting from Lya line
    constexpr const float power_invlambda1 = pow(lambda1, -3.46);

    if (wavelength < lambda1 * z_plus_one)
    {
        teff1 += 0.0037 * x_power * power_invlambda1;
    }
    else return teff1; // because if x >= lambda1 * z_plus_one, it's going to also going to fail all
                       // of the following tests because these are in order of decreasing lambda

    // lyman beta (3->1)
    constexpr const float lambda2 = 1026.0;
    constexpr const float power_invlambda2 = pow(lambda2, -3.46);

    if (wavelength < lambda2 * z_plus_one)
    {
        teff1 += 0.00177 * x_power * power_invlambda2;
    }
    else return teff1;

    // lyman gamma (4->1)
    constexpr const float lambda3 = 973.0;
    constexpr const float power_invlambda3 = pow(lambda3, -3.46);

    if (wavelength < lambda3 * z_plus_one)
    {
        teff1 += 0.00106 * x_power * power_invlambda3;
    }
    else return teff1;

    // lyman (5->1)
    constexpr const float lambda4 = 950.0;
    constexpr const float power_invlambda4 = pow(lambda4, -3.46);

    if (wavelength < lambda4 * z_plus_one)
    {
        teff1 += 0.000584 * x_power * power_invlambda4;
    }
    else return teff1;

    // lyman (6->1)
    constexpr const float lambda5 = 938.1;
    constexpr const float power_invlambda5 = pow(lambda5, -3.46);

    if (wavelength < lambda5 * z_plus_one)
    {
        teff1 += 0.00044 * x_power * power_invlambda5;
    }
    else return teff1;

    // lyman (7->1)
    constexpr const float lambda6 = 931.0;
    constexpr const float power_invlambda6 = pow(lambda6, -3.46);

    if (wavelength < lambda6 * z_plus_one)
    {
        teff1 += 0.00040 * x_power * power_invlambda6;
    }
    else return teff1;

    // lyman (8->1)
    constexpr const float lambda7 = 926.5;
    constexpr const float power_invlambda7 = pow(lambda7, -3.46);

    if (wavelength < lambda7 * z_plus_one)
    {
        teff1 += 0.00037 * x_power * power_invlambda7;
    }
    else return teff1;

    // lyman (9->1)
    constexpr const float lambda8 = 923.4;
    constexpr const float power_invlambda8 = pow(lambda8, -3.46);

    if (wavelength < lambda8 * z_plus_one)
    {
        teff1 += 0.00035 * x_power * power_invlambda8;
    }
    else return teff1;

    // lyman (10->1)
    constexpr const float lambda9 = 921.2;
    constexpr const float power_invlambda9 = pow(lambda9, -3.46);

    if (wavelength < lambda9 * z_plus_one)
    {
        teff1 += 0.00033 * x_power * power_invlambda9;
    }
    else return teff1;

    // lyman (11->1)
    constexpr const float lambda10 = 919.6;
    constexpr const float power_invlambda10 = pow(lambda10, -3.46);

    if (wavelength < lambda10 * z_plus_one)
    {
        teff1 += 0.00032 * x_power * power_invlambda10;
    }
    else return teff1;

    // lyman (12->1)
    constexpr const float lambda11 = 918.4;
    constexpr const float power_invlambda11 = pow(lambda11, -3.46);

    if (wavelength < lambda11 * z_plus_one)
    {
        teff1 += 0.00031 * x_power * power_invlambda11;
    }

    // assumed that higher order contributions are essentially negligible at this point
    return teff1; // effective optical depth (NOT TRANSMISSION)
}


const real_t inline madau_teff2(real_t x, real_t z)
{
    real_t teff2 = 0.0;

    const float zlambda = 912.0 * (1.0 + z); // declaring the lyman limit
    const constexpr float inv_912 = 1.0 / 912.0;

    if (x < zlambda)
    {
        const float xc = x * inv_912; // the (1+z) factor assuming this was a 912 line
        const float xc_cubed = xc * xc * xc;
        const float xem = 1.0 + z; // the actual observed (1+z) factor
        teff2 = 0.25 * xc_cubed * (pow(xem, 0.46) - pow(xc, 0.46));
        teff2 += teff2 + 9.4 * pow(xc, 1.5) * (pow(xem, 0.18) - pow(xc, 0.18));
        teff2 -= 0.7 * xc_cubed * (pow(xc, -1.32) - pow(xem, -1.32));
        teff2 -= 0.023 * (pow(xem, 1.68) - pow(xc, 1.68));
    }

    // avoid negative numbers that might result from approximation:

    return (teff2 > 0.0) ? teff2 : 0.0; // effective optical depth
}


const real_t inline madau_teff(real_t x, real_t z)
{
    // add 912 - 1216 (teff1) and < 912 (teff2) portions for optical depth
    const real_t tf = madau_teff1(x, z) + madau_teff2(x, z);

    // return effective TRANSMISSION (converted from optical depth)
    return exp(-tf);
}


void apply_madau(const real_t* wave, real_t* flux, const real_t z, const int Np)
{
    // multiply flux by effective transmission
    for (int i = 0; i < Np; ++i)
    {
        flux[i] *= madau_teff(wave[i] * (1.0 + z), z);
    }
}


int waveidx(real_t lmin, real_t dlogwave, real_t l)
{
    int i;
    real_t idx, nu, nuMax;
    const real_t c = 299792458.0; // speed of light (m/s)

    // as all the integrals and whatnot are done in frequency space, all values are converted to nu whenever possible
    nuMax = c / (lmin * 1e-10); // max nu, based on input minimum lambda (in Hz)
    nu = c / (l * 1e-10); // nu at given lambda (Hz)

    // figuring out how many discrete wavelength channels, separated by dlogwave, are spanned by lmin-l (nuMax-nu)
    idx = (log10(nuMax) - log10(nu)) / dlogwave;

    i = idx;

    return (i); // returns this spacing/number
}


void set_sampling(real_t** Swave,
                  int* Ns,
                  int Nsed,
                  real_t dz,
                  real_t Rmax,
                  int* outNwave,
                  real_t* outlmin,
                  real_t* outdlogwave)
{
    int Nwave = 0;
    int Npoints = 10; // number of points used to determine spacing (i.e. coarseness)
    real_t lmin, lmax, dlogwave;
    // const real_t c = 299792458.0;

    // determine the spacing of wavelength grid, using either redshift spacing (dz)
    // or the highest desired filter resolution (Rmax)
    if ((1.0 / dz) > Rmax)
    {
        dlogwave = 1.0 / ((1.0 / dz) * Npoints);
    }
    else
    {
        dlogwave = 1.0 / (Rmax * Npoints);
    }

    // figuring out min/max wavelength

    // declaring large starting values
    lmin = 9e99;
    lmax = -9e99;

    // goes through all SEDs, figures out min/max wavelengths
    for (int i = 0; i < Nsed; ++i)
    {
        if (Swave[i][0] < lmin)
        {
            lmin = Swave[i][0];
        }

        if (Swave[i][Ns[i]] > lmax)
        {
            lmax = Swave[i][Ns[i] - 1];
        }

        /* fprintf(stderr,"lmax = %lf %d\n",lmax,Ns[i]); */
    }

    // shift by one percent to avoid edge effects
    lmax = 0.99 * lmax;
    lmin = 1.01 * lmin;

    // number of wavelengths between max and min, spaced by dlogmax

    Nwave = (log10(lmax) - log10(lmin)) / dlogwave;

    // set the number of wavelengths, minimum wavelength, and log spacing
    *outNwave = Nwave;
    *outlmin = lmin;
    *outdlogwave = dlogwave;
}


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
                   real_t dlogwave)
{
    int i = 0;
    int j = 0;
    real_t lmax, l;
    real_t nuMax, nu, lnu;
    const real_t c = 299792458.0; // in (m/s)

    if (lmin == 0.0)
    {
        cerr << "in sed.cpp, lmin cannot be zero\n";
        exit(1);
    }

    real_t dlambda = pow(10, dlogwave) - 1;

    // shifting over to frequency space, declaring edge
    nuMax = c / (lmin * 1e-10);

    // set spacing in equal dlog(nu)
    for (int i = 0; i < Nwave; ++i)
    {
        lnu = log10(nuMax) - i * dlogwave; // stepping down in log(nu), up in log(lambda)
        nu = pow(10, lnu); // shifting back to non-log space
        l = c * 1e10 / nu; // converting back to wavelength (A)
        outSwave[i] = l;
        outnu[i] = nu;
        outdnu[i] = dlambda * nu;
    }

    for (int i = 0; i < Nsed; ++i)
    {
        // setting up the spline
        gsl_spline* Sinterp = gsl_spline_alloc (gsl_interp_cspline, Ns[i]);
        gsl_spline_init(Sinterp, Swave[i], Sflux[i], Ns[i]);
        gsl_interp_accel* Sacc =  gsl_interp_accel_alloc();

        // running through all desired wavelengths
        for (int j = 0; j < Nwave; ++j)
        {
            l = outSwave[j]; // grab the relevant wavelength

            // check if the wavelength runs past the blueward edge of the SED's given wavelengths
            if (l < Swave[i][0])
            {

                outSflux[i][j] = Sflux[i][0]; // if so, just assign the output flux at given wavelength to be equal to flux at blueward edge of SED
            }

            // otherwise, check to see if wavelength runs past the redward edge of the SED's given wavelengths
            else
            {
                if (l > Swave[i][Ns[i] - 1])
                {
                    outSflux[i][j] = Sflux[i][Ns[i] - 1]; // if so, assign output flux at the given wavelength equal to flux at redward edge of SED
                }

                // if flux is within bounds, assign flux based on spline interpolation at given wavelength
                else
                {
                    outSflux[i][j] = gsl_spline_eval(Sinterp, l , Sacc);
                }
            }

            // set negative fluxes to 0
            if (outSflux[i][j] < 0)
            {
                outSflux[i][j] = 0;
            }
        }

        // free up spline interpolation for next SED
        gsl_spline_free(Sinterp);
        gsl_interp_accel_free(Sacc);
    }
}


void setup_resample_filters(real_t* outSwave,
                            real_t lmin,
                            real_t dlogwave,
                            real_t** Fwave,
                            real_t** Ftrans,
                            int NF,
                            int* NFwave,
                            int* Fstart,
                            int* F_nelem)
{
    for (int i = 0; i < NF; ++i)
    {
        const int istart = waveidx(lmin, dlogwave, Fwave[i][0]); // get the number of wavelengths between lmin and the start of a given filter for dlogwave
        const int iend = waveidx(lmin, dlogwave, Fwave[i][NFwave[i] - 1]); // get the number of wavelengths between lmin and the end of a given filter for dlogwave
        F_nelem[i] = iend - istart; // takes the difference to get the number of elements in the filter
        Fstart[i] = istart; // tells the program how far along Swave to start sampling for a given filter
    }
}


void make_reddening_vec(real_t* outSwave, real_t** outRvec, int Nrlaw, int Nwave)
{
    for (int i = 1; i <= Nrlaw; ++i)
    {
        make_dustvec(outSwave, outRvec[i - 1], i, Nwave);
    }
}


void resample_filters(real_t* outSwave,
                      real_t** Fwave,
                      real_t** Ftrans,
                      int NF,
                      int* NFwave,
                      int* Fstart,
                      int* outFnelem,
                      real_t** outFtrans)
{
    // lets loop over the filters and do some integrals!
    for (int i = 0; i < NF; ++i)
    {
        // set up the filter spline
        gsl_spline* Finterp = gsl_spline_alloc (gsl_interp_cspline, NFwave[i]);
        gsl_spline_init(Finterp, Fwave[i], Ftrans[i], NFwave[i]);
        gsl_interp_accel* Facc =  gsl_interp_accel_alloc();

        for (int j = 0; j < outFnelem[i]; ++j)
        {
            const uint k = j + Fstart[i]; // run through the filter, with an amount offset by the appropriate Fstart

            // check if the SED wavelength at the appropriate element is less than the actual wavelength at the blue edge of the filter
            if (outSwave[k] < Fwave[i][0])
            {
                outFtrans[i][j] = 0; // if so, set output transmission to 0
            }

            // otherwise, check if the SED wavelength at the appropriate element is greater than the actual wavelength at the red edge of the filter
            else
            {
                if (outSwave[k] > Fwave[i][NFwave[i] - 1])
                {
                    outFtrans[i][j] = 0; // if so, set output transmission to 0
                }

                // otherwise, set the transmission via the best-fit spline to the filter at the appropriate wavelength
                else
                {
                    outFtrans[i][j] = gsl_spline_eval(Finterp, outSwave[k], Facc);
                }
            }

            // getting rid of any negative values
            if (outFtrans[i][j] < 0)
            {
                outFtrans[i][j] = 0;
            }

            /* fprintf(stderr,"%d %d %lf %lf\n",j,k,outSwave[k],outFtrans[i][j]); */
        }

        // freeing up spline
        gsl_spline_free(Finterp);
        gsl_interp_accel_free(Facc);
    }
}


void integrate_sed(real_t dlogwave,
                   const real_t* Swave,
                   const float* dnu_nu,
                   real_t* Sflux,
                   int Ns,
                   const int* Fnelem,
                   const int* Fstart,
                   real_t** Ftrans,
                   int NF,
                   const real_t* Rvec,
                   real_t z,
                   real_t Ebv,
                   real_t* Sflux_corrected,
                   real_t* FluxOut)
{
    constexpr const real_t invc = 1.0 / 1.299792458; // (m/s)

    // for each point, apply reddening vector
    for (int i = 0; i < Ns; ++i)
    {
        const real_t s_wavei = Swave[i] * 1e-10; // angstroms to meters

        const real_t sflux = Sflux[i];

        Sflux_corrected[i] = (sflux > 0.0) ? Sflux[i] * pow(10, -0.4 * Rvec[i] * Ebv) * s_wavei * s_wavei * invc
                                           : 0.0;  // any negative values are set to 0.0
    }

    apply_madau(Swave, Sflux_corrected, z, Ns); // apply extra extinction for H-clouds based on madau (1995, 1996)

    const int ioff = log10(1 + z) / dlogwave; // applying an offset for redshift (as this is in log space, redshift is additive!)

    // for each filter
    for (int i = 0; i < NF; ++i)
    {
        // integrate each filter

        const int filter_elements_i = Fnelem[i];

        const int starting_offset = Fstart[i] - ioff;

        const real_t* f_transfer_i = Ftrans[i];

        real_t sum = 0;
        real_t weight = 0;

        if (filter_elements_i + starting_offset >= Ns)
        {
            std::cerr << i << "\n";
            std::cerr << filter_elements_i << "\n";
            std::cerr << starting_offset << "\n";
            exit(101);
        }

        for (int j = 0; j < filter_elements_i; ++j)
        {
            // shift starting point for integration based on filter and redshift (k = j + Fstart[i] - ioff)
            const int k = j + starting_offset;
            const real_t norm = f_transfer_i[j] * dnu_nu[k];
            sum += Sflux_corrected[k] * norm; // iteratively performs integral
            weight += norm; // overall normalization
        }

        FluxOut[i] = sum / weight; // normalizing integral to get appropriate flux
    }
}

void ComputeSed(const real_t* Swave,
                real_t* Sflux,
                const int Ns,
                const real_t* Rvec,
                real_t z,
                real_t Ebv,
                real_t* fluxes)
{
    constexpr const real_t invc = 1.0 / 1.299792458; // (m/s)

    // for each point, apply reddening vector
    for (int i = 0; i < Ns; ++i)
    {
        const real_t s_wavei = Swave[i] * 1e-10; // angstroms to meters

        const real_t sflux = Sflux[i];

        fluxes[i] = (sflux > 0.0) ? Sflux[i] * pow(10, -0.4 * Rvec[i] * Ebv) * s_wavei * s_wavei * invc
                                           : 0.0;  // any negative values are set to 0.0
    }

    apply_madau(Swave, fluxes, z, Ns); // apply extra extinction for H-clouds based on madau (1995, 1996)

}
