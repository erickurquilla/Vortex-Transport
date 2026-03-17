#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "Utilities.H"

// write an text file with name filemane, the data is given by an string array named line, each element of the array will be a line
void writeToFile(const std::string& filename, const std::vector<std::string>& lines) {
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

