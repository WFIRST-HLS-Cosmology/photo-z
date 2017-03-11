/*! \file
 *
 * Provides basic tools for reading data stored in text files.
*/

#ifndef BASIC_UTILS_H
#define BASIC_UTILS_H

#ifdef HAS_FOLLY
    #include <folly/FBVector.h>
    #define VECTOR folly::fbvector
#else
    #include <vector>
    #define VECTOR std::vector
#endif

#include <string>

typedef VECTOR<std::string> StringList;
typedef VECTOR<double> ColumnData;

/*! Reads a text file and stores the file's contents as a vector of the non-empty strings.
 *  \param filename is the name of the file that will be read */

/*!
 * Reads a text file and stores the file's contents as a vector of the non-empty strings.
 * \param filename is the name of a text that will be read.
 * \return Returns a StringList containing one entry for each non-blank line of the file.
 */
StringList ReadLines(const char* filename);


/*! Reads all of the files listed in the file_list and saves the content of each file as vectors of
 * double precision floating point numbers. column1 contains the numbers in the first column of
 * all of the files, such that column1[n][m] contains the mth entry in the firsth column of the nth
 * file listed in the file_list. column2 contains the corresponding value from the second column of
 * each file. \todo this should be replaced with a more general function that allows for more data
 * locality once the other functions have been modified. */
void ReadTwoColumnFiles(const char* file_list,
                        VECTOR<ColumnData>& column1,
                        VECTOR<ColumnData>& column2);


#endif // BASIC_UTILS_H
