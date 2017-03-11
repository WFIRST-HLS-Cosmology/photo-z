/*! \file

Contains the main function for the catalog photometry fitting code (fitcat).
*/

#include "grid_tools.hpp"

#include <cstdlib>
#include <string>
#include <iostream>
#include <fmt/format.h>

using namespace fmt::literals;

int main(int argc, char* argv[])
{    
    std::string mkgrid_usage = "{} mkgrid catalog filter_list template_list output_file"_format(argv[0]);

    std::string fit_usage = "{} fit catalog grid output_file"_format(argv[0]);

    if (argc == 2)
    {
        std::string task_str(argv[1]);

        if (task_str == "fit") { fmt::print("\nUsage: {}\n", fit_usage); exit(1); }

        if (task_str == "mkgrid") { fmt::print("\nUsage: {}\n", mkgrid_usage); exit(1); }

    }
    else if (argc < 5 || argc > 6)
    {
        fmt::print("\nUsage:\n\n {}\n\nOR\n\n {}\n", mkgrid_usage, fit_usage);
        exit(1);
    }

    std::string task_str(argv[1]);

    if (task_str == "fit")
    {
        const char* catalog_filename = argv[2];
        const char* model_grid_filename = argv[3];
        const char* output_filename = argv[4];

        FindBestFits(catalog_filename, model_grid_filename, output_filename);
    }
    else if (task_str == "mkgrid")
    {
        const char* catalog_filename = argv[2];
        const char* filterlist_filename = argv[3];
        const char* templatelist_filename = argv[4];
        const char* output_filename = argv[5];

        MakePhotometryFittingGrid(catalog_filename,
                                  filterlist_filename,
                                  templatelist_filename,
                                  output_filename);
    }
    else
    {
        fmt::print_colored(fmt::RED, "Error: The task must be 'fit' or 'mkgrid.' You entered {}.\n", task_str);
        exit(1);
    }

    return 0;
}




