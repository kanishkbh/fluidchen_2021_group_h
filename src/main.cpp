#include <iostream>
#include <string>

#include"Communication.hpp"
#include "Case.hpp"

int main(int argn, char **args) {

    if (argn > 1) {
        std::string file_name{args[1]};
        int rank = 1; // Dummy , to get rank from MPI 
        Case problem(file_name, atoi(args[1]),rank);
        problem.simulate();

    } else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
