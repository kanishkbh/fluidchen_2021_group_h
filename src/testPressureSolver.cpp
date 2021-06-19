#include <iostream>
#include <string>
#include <chrono>
#include <mpi.h>

#include "PressureTestCase.hpp"

int main(int argn, char **args) {

    // Initialize MPI
    MPI_Init(nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    if (argn > 1 && rank == 0) {
        std::string file_name{args[1]};
        PressureTestCase problem(file_name, argn, args);
        
        
        
        auto t_start = std::chrono::steady_clock::now();
        auto&& res = problem.pressure_solve(500);
        auto t_end  = std::chrono::steady_clock::now();
        if (rank == 0)
        {
            std::cout << "Time elapsed : " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << " ms." << std::endl;
            for (double r : res)
                std::cout << r << std::endl;
        }

    } else if (rank == 0) {
        std::cout << "Error: No input file is provided to test_pressure_solver. Creates a Lid Driven Cavity of size N." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen exemple_cases/Lid.dat" << std::endl;
    }

    MPI_Finalize();
}
