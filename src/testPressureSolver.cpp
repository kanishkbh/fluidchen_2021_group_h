#include <iostream>
#include <string>
#include <chrono>
#include <mpi.h>
#include <map>
#include <memory>

#include "PressureTestCase.hpp"

int main(int argn, char **args) {

    const double tol = 1e-5;

    // Initialize MPI
    MPI_Init(nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);

    if (argn > 1) {
        std::string file_name{args[1]};
        PressureTestCase problem(file_name, argn, args);

        std::map<std::string, std::shared_ptr<PressureSolver>> solvers;
        solvers["SOR"] = std::make_shared<SOR>(1.7);
        solvers["Jacobi"] = std::make_shared<Jacobi>(1.0);
        solvers["SD"] = std::make_shared<SD>();
        solvers["CG"] = std::make_shared<CG>();
        if (num_processes == 1)
        {
            solvers["CG_GS"] = std::make_shared<CG_GS>();
        }

        for (auto& pair : solvers) {
            if (rank == 0)
                {
                    std::cout << "Launching solver " << pair.first << "... ";
                }
            auto t_start = std::chrono::steady_clock::now();
            auto&& res = problem.pressure_solve(tol, *pair.second);
            auto t_end  = std::chrono::steady_clock::now();
            if (rank == 0)
                {
                    std::cout << " Needed " << res.size() << " iterations and " 
                    << std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count() << " ms." << std::endl;
                }
            

        }
          



    } else if (rank == 0) {
        std::cout << "Error: No input file is provided to test_pressure_solver. Creates a Lid Driven Cavity of size N." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen exemple_cases/Lid.dat" << std::endl;
    }

    MPI_Finalize();
}
