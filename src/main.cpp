#include <iostream>
#include <string>
#include <mpi.h>
#include "Case.hpp"

int main(int argn, char **args) {

    if (argn > 2) {
        
        // Initialize MPI
        MPI_Init(&argn, &args+1); // this needs to change

        // Get number of processors
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // Get this processor rank in this communication world
        int my_rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

        std::string file_name{args[2]};
        Case problem(file_name, size, my_rank);
        problem.simulate();

        // Finalize MPI
        MPI_Finalize();

    } 
    else {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }
}
