#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines);

// clean and create a directory
void clean_create_directory(const std::string& dirname);

#endif // UTILITIES_H