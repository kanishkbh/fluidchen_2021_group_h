#include "PressureTestCase.hpp"
#include "Enums.hpp"
#include "PressureSolver.hpp"
#include "Communication.hpp"

#include <mpi.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cassert>

PressureTestCase::PressureTestCase(std::string file_name, int argn, char **args) : Case(file_name, argn, args) {

}

std::vector<double> PressureTestCase::pressure_solve(double tol, PressureSolver& solver) {
    double t = 0.0;
    double dt = _field.calculate_dt(_grid);
    int timestep = 0;
    int output_id = 0;
    double output_counter = 0.0;
    double res;

    /* Get the total number of fluid cells to normalize the residuals */
    int LOCAL_NUMBER_OF_CELLS = _grid.fluid_cells().size();
    int TOTAL_NUMBER_OF_CELLS_REDUCED;
    Communication::communicate_sum_int(&LOCAL_NUMBER_OF_CELLS, &TOTAL_NUMBER_OF_CELLS_REDUCED);
    // Initialization
    std::vector<double> out;



    // Apply the Boundary conditions
    for (auto &boundary_ptr : _boundaries) {
        boundary_ptr->apply(_field);
    }

    // Update temperatures
    // (Temperature is still used in calculate_fluxes, but since T is initialized to a constant, it doesn't matter)
    if (_use_energy) {
        _field.calculate_T(_grid);
        Communication::communicate_all(_field.t_matrix(), MessageTag::T);
    }

    // Fluxes (with *new* temperatures)
    _field.calculate_fluxes(_grid);

    // Communicate F and G
    Communication::communicate_all(_field.f_matrix(), MessageTag::F);
    Communication::communicate_all(_field.g_matrix(), MessageTag::G);

    // Poisson Pressure Equation
    _field.calculate_rs(_grid);
    /* First step of the solver */
    solver.init(_field, _grid, _boundaries);

    /*Reset pressure (otherwise the benchmark is broken !) */
    _field.p_matrix() = Matrix<double>(_field.p_matrix().imax(), _field.p_matrix().jmax());

    do {
        res = solver.solve(_field, _grid, _boundaries);
            /* Compute TOTAL residual */
            double res_reduction;
            Communication::communicate_sum_double(&res, &res_reduction);
            res = res_reduction / TOTAL_NUMBER_OF_CELLS_REDUCED;
            res = std::sqrt(res);

            out.push_back(res);
    } while (out.size() < _max_iter && res > tol);

    return out;
}