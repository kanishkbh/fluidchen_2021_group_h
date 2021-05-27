#include <iostream>
#include <string>
#include <mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {

    // Initialize MPI
    MPI_Init(nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argn > 1) {
        std::string file_name{args[1]};
        Case problem(file_name, argn, args);
        problem.simulate();

    } else if (rank == 0) {
        std::cout << "Error: No input file is provided to fluidchen." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    }

    MPI_Finalize();
}
