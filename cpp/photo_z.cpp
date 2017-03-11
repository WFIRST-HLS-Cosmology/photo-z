/*! \file
 *
 *  Contains the main function for the redshift estimation code (photo-z).
 *
 * The binary expects (1) a model_grid file, created using the output of make-grid (compiled from
 * make_grid.cpp), (2) a file containing a list of SEDs that will be examined, (3) the name of the
 * file in which to write the output data (estimated redshifts for each SED in argument 2), and,
 * optionally, the name of a file that will hold probability distribution functions for each input
 * SED.
 *
 * */

#include "fitting_tools.hpp"

#include <omp.h>

#include <iostream>
#include <cstdlib>

using namespace std;


int main(int argc, char* argv[])
{
    if (argc < 5 || argc > 6)
    {
        cerr << "Usage: " << argv[0] << " model_grid SED_catalog template_probability_file output_redshift_file [pdfs_file]\n";
        return -1;
    }

    const bool user_wants_to_ouput_the_pdfs = (argc == 6);

    const char* model_grid_file = argv[1];
    const char* object_sed_file = argv[2];
    const char* template_probability_file = argv[3];
    const char* output_file = argv[4];
    const char* pdfs_file = (user_wants_to_ouput_the_pdfs) ? argv[5] : "";

    cerr << "Starting the redshift estimation program, " << argv[0] << "\n";

    // read the input files and compute probability distributions

    RedshiftPDF redshift_pdf(model_grid_file, object_sed_file, template_probability_file);

    // compute and output moments of the probability distributions and save the results. Output -9
    // if the probability distribution was not useful (for instance, if it consisted only of zeros).

    redshift_pdf.SaveRedshifts(output_file);

    // if the pdfs_file argument is provided on the command line, save the pdfs.

    if (user_wants_to_ouput_the_pdfs) { redshift_pdf.SavePDFs(pdfs_file); }

    cerr << argv[0] << " is finished!\n";

    return 0;
}

/*! \mainpage SPHEREx Photo-z Documentation

\section welcome Welcome

This document contains details about the source code for the photometric redshift estimation pipeline
for SPHEREx.

If you wish to understand the code as quickly as possible, I recommend using the Qt Creator
integrated development environment in addition to this documentation. See the following section for
more information.

\subsection QtCreator Using Qt Creator

If you use the Qt Creator integrated development environment to work on the code (recommended), you
can integrate this source code documentation, which makes it even easier to understand the source.

\par Requirements

First, you will need to install Doxygen, Graphviz, and a \f$\mbox{\LaTeX}\f$ system, such as
TeXLive on your local system, in addition to Qt Creator and qt4-dev-tools.

\par Step 1:

In the \c src/ directory, issue the command

\verbatim
doxygen
\endverbatim

This command will create the latest version of this documentation in the \c
Docs/html/ directory.

\par Step 2:

Navigate to \c Docs/html/ and run the command

\verbatim
qhelpgenerator index.qhp
\endverbatim

This generates the file \c index.qch, which can be used by Qt Creator to connect this documentation
with the source.

\par Step 3:

Open Qt Creator and navigate to Tools > Options > Help. Then click the Documentation tab. Click
Add... and locate the file \c Docs/html/index.qch.

\par Step 4:

Restart Qt Creator and load the photo-z project file, \c src/photoz.creator Now you can quickly
access the documentation for any documented class, function, variable, typedef, etc. by hovering over the name
of the entity in the source code and pressing the F1 key.

\par Note:

Using Qt Creator 2.8 or later and Doxygen 1.8.3.1 or later is recommended, since earlier versions
did not support enough C++11 syntax to be useful for all parts of the source. If you wish to
use Eclipse or XCode as your IDE, you can modify the Doxygen file, \c src/Doxyfile.

*/
