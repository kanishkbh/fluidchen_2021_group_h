#include <iostream>
#include <string>
#include <chrono>
#include <mpi.h>

#include "Case.hpp"

int main(int argn, char **args) {

    // Initialize MPI
    MPI_Init(nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    if (argn > 1 && rank == 0) {
        int N = std::atoi(args[1])
        
        /* Initialize the problem */
        auto grid = Grid(_geom_name, domain, _left_neighbor_rank != -1, _right_neighbor_rank != -1,
                 _top_neighbor_rank != -1, _bottom_neighbor_rank != -1);

        auto field = Fields(nu, dt, tau, grid.domain().size_x, grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);
        auto discretization = Discretization(domain.dx, domain.dy, gamma);
        auto pressure_solver = std::make_unique<SOR>(omg);

        auto t_start = std::chrono::steady_clock::now();
        problem.simulate();
        auto t_end  = std::chrono::steady_clock::now();
        if (rank == 0)
            std::cout << "Time elapsed : " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count() << " ms." << std::endl;

    } else if (rank == 0) {
        std::cout << "Error: No input file is provided to test_pressure_solver. Creates a Lid Driven Cavity of size N." << std::endl;
        std::cout << "Example usage: /path/to/fluidchen 50" << std::endl;
    }

    MPI_Finalize();
}
