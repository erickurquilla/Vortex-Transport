#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Utilities.H"
#include "Element.H"
#include "Profiling.H"

PROFILE_DECLARE("writeToFile");
PROFILE_DECLARE("clean_create_directory");
PROFILE_DECLARE("write_output");

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines) {
    PROFILE_SCOPE("writeToFile");
    // Open the file for writing
    std::ofstream outfile(filename);

    // Check if the file is open
    if (!outfile.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return;
    }

    // Write each line to the file
    for (const std::string& line : lines) {
        outfile << line << std::endl;
    }

    // Close the file
    outfile.close();
}

void clean_create_directory(const std::string& dirname){

    PROFILE_SCOPE("clean_create_directory");

    // Clean directory
    std::string dirPath = dirname;
    std::string command = "rm -rf " + dirPath; // Remove directory and its contents
    int status = system(command.c_str());

    if (status == 0) {
        // std::cout << "Directory cleaned successfully: " << dirPath << std::endl;
    } else {
        std::cerr << "Failed to clean directory: " << dirPath << std::endl;
        exit(EXIT_FAILURE);
    }

    // Now directory to store simulation for current step
    command = "mkdir -p " + dirPath;
    status = system(command.c_str());

    if (status == 0) {
        // std::cout << "Directory created successfully: " << dirPath << std::endl;
    } else {
        std::cerr << "Failed to create directory: " << dirPath << std::endl;
        exit(EXIT_FAILURE);
    }

}

// write output data of all elements for the given step (parallelized with OpenMP)
void write_output(Element* elements, const int& total_num_elements, const int& step_num){

    PROFILE_SCOPE("write_output");

    #ifdef VORTEX_USE_OPENMP
    #pragma omp parallel for default(shared) schedule(runtime)
    #endif
    for (int i = 0; i < total_num_elements; ++i) {
        elements[i].write_data(step_num); // write data
    }

}

