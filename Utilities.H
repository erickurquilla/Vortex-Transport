#ifndef UTILITIES_H
#define UTILITIES_H

#include <vector>
#include <string>

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines);

// clean and create a directory
void clean_create_directory(const std::string& dirname);

class Element;

// write output data of all elements for the given step (parallelized with OpenMP)
void write_output(Element* elements, const int& total_num_elements, const int& step_num);

#endif // UTILITIES_H