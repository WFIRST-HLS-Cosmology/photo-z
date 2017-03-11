/*! \file
 *
 * Definitions of functions for loading data from text files.
*/

#include "basic_utils.hpp"

#include <boost/algorithm/string.hpp>

#include <fstream>
#include <iostream>

StringList ReadLines(const char* filename)
{
    std::ifstream input_fs(filename);

    if ( !input_fs.is_open() )
    {
        std::cerr << "could not open " << filename << std::endl;
        exit(1);
    }

    std::string line;

    StringList list_of_lines;

    std::cerr << "Reading " << filename << std::endl;

    while(std::getline(input_fs, line))
    {
        boost::trim(line);

        if (line != "") { list_of_lines.push_back(line); }
    }

    return list_of_lines;
}


void ReadTwoColumnFiles(const char* file_list,
                    VECTOR<ColumnData>& column1,
                    VECTOR<ColumnData>& column2)
{
    std::string file_list_str(file_list);

    // check whether the absolute or relative path was entered

    std::string filepath;

    if (file_list_str.find('/') != std::string::npos) // if the file_list_string contains "/"
    {
        uint end_of_path = file_list_str.find_last_of('/');

        filepath = file_list_str.substr(0, end_of_path) + "/";
    }
    else // relative path is the local / current directory
    {
        std::cerr << "Error: When specifying a file containing a list of files, please provide "
                     "the path. Examples: /path/to/file.list or ../relative/path/file.list.\n";
        exit(1);
    }

    // open the list of files and store contents

    StringList file_list_contents = ReadLines(file_list);

    for (auto& file : file_list_contents)
    {
        // open file, store lines, parse the lines to fill templates

        // we are only interested in the last entry on the line; all other text is ignored

        StringList filename_strings;

        boost::split(filename_strings, file,
                     boost::algorithm::is_any_of("\t "),
                     boost::algorithm::token_compress_on);

        std::string filename = filepath + filename_strings.back(); // the final column

        StringList file_contents = ReadLines(filename.c_str());

        ColumnData col1;
        ColumnData col2;

        for (auto& line : file_contents)
        {
            StringList strings;

            boost::split(strings, line,
                         boost::algorithm::is_any_of("\t "),
                         boost::algorithm::token_compress_on);

            col1.push_back(std::stod(strings[0])); // column 1
            col2.push_back(std::stod(strings[1])); // column 2
        }

        column1.push_back(col1);
        column2.push_back(col2);
    }
}
