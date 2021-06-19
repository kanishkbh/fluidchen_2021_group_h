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

PressureTestCase::PressureTestCase(std::string file_name, int argn, char **args) {

    
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    // _geom_name = file_name.substr(0, file_name.length()-4) + _geom_name;
    
    std::ifstream file(file_name);
    double nu;           /* viscosity   */
    double UI;           /* velocity x-direction */
    double VI;           /* velocity y-direction */
    double TI;           /* Initial temperature */
    double beta;         /* Thermal expansion coefficient */
    double alpha{0};     /* Heat diffusivity */
    double PI;           /* pressure */
    double GX;           /* gravitation x-direction */
    double GY;           /* gravitation y-direction */
    double xlength;      /* length of the domain x-dir.*/
    double ylength;      /* length of the domain y-dir.*/
    double dt;           /* time step */
    int imax;            /* number of cells x-direction*/
    int jmax;            /* number of cells y-direction*/
    double gamma;        /* uppwind differencing factor*/
    double omg;          /* relaxation factor */
    double tau;          /* safety factor for time step*/
    int itermax;         /* max. number of iterations for pressure per time step */
    double eps;          /* accuracy bound for pressure*/
    double UIN;          /* Inlet horizontal velocity */
    double VIN;          /* Inlet vertical velocity */
    double in_temp;      /* Inlet (Dirichlet) Temperature */
    int use_pressure{0}; /* If non-zero, use pressure BC instead of inflow velocity*/
    std::string energy_eq; /* If "on", enable heat transfer */

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "geo_file") { 
                    file >> _geom_name;
                    _geom_name = file_name.substr(0, file_name.find_last_of('/')+1) + _geom_name;
                }
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "in_temp") file >> in_temp;
                if (var == "TI") file >> TI;
                if (var == "alpha") file >> alpha;
                if (var == "beta") file >> beta;
                if (var == "T_IN") file >> _T_IN;
                if (var == "moving_wall_temp") file >> _moving_wall_temp;
                if (var == "wall_temp_3") file >> _fixed_wall_temp[cell_type::FIXED_WALL_3];
                if (var == "wall_temp_4") file >> _fixed_wall_temp[cell_type::FIXED_WALL_4];
                if (var == "wall_temp_5") file >> _fixed_wall_temp[cell_type::FIXED_WALL_5];
                if (var == "wall_temp_6") file >> _fixed_wall_temp[cell_type::FIXED_WALL_6];
                if (var == "wall_temp_7") file >> _fixed_wall_temp[cell_type::FIXED_WALL_7];
                if (var == "use_pressure_input") file >> use_pressure;
                if (var == "PIN") file >> _P_IN;
                if (var == "energy_eq") file >> energy_eq;
                if (var == "iproc") file >> _iproc;
                if (var == "jproc") file >> _jproc;
            }
        }
    }
    else {
        std::cerr << "Couldn't open file " << file_name << ". Aborting." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Check number of processes matches. If only one process, no partition.
    int num_processes;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    assert (_iproc * _jproc == num_processes || num_processes == 1);
    if (num_processes == 1) {
        _iproc = 1;
        _jproc = 1;
    }

    
    file.close();

    // Update the boolean for pressure input and heat equation
    if (_rank == 0) {
        if (energy_eq == "on") {
            std::cout << "Enabling heat transfer." << std::endl;
        } else {
            std::cout << "Heat transfer disabled." << std::endl;
        }
    }
    _use_energy = energy_eq == "on";
    _use_pressure_input = (use_pressure != 0);

    //-----------------------------------------------------------------------------------------------------------
    std::map<int, double> wall_vel;
    //Should be deprecated soon
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }
//-----------------------------------------------------------------------------------------------------------
    // Set file names for output directory
    set_file_names(file_name);
//-----------------------------------------------------------------------------------------------------------
    // Build up the domain
    Domain domain;
    domain.dx = xlength / (double)imax;
    domain.dy = ylength / (double)jmax;
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;
    build_domain(domain, imax, jmax);
    //-----------------------------------------------------------------------------------------------------------
    // Load the geometry file
    _grid = Grid(_geom_name, domain, _left_neighbor_rank != -1, _right_neighbor_rank != -1,
                 _top_neighbor_rank != -1, _bottom_neighbor_rank != -1);

    //-----------------------------------------------------------------------------------------------------------    
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);
//-----------------------------------------------------------------------------------------------------------
    _discretization = Discretization(domain.dx, domain.dy, gamma);
//-----------------------------------------------------------------------------------------------------------

/* TODO : PICK SOLVER  */
    _pressure_solver = std::make_unique<CG>();
//-----------------------------------------------------------------------------------------------------------
    _max_iter = itermax;
//-----------------------------------------------------------------------------------------------------------
    _tolerance = eps;
//-----------------------------------------------------------------------------------------------------------
    // Construct boundaries 

    _wall_velocity = 0; //TODO
    _u_in = UIN;
    _v_in = VIN;
    _p_i = PI;
    setupBoundaryConditions();

}

std::vector<double> PressureTestCase::pressure_solve(unsigned N) {
    double t = 0.0;
    double dt = _field.calculate_dt(_grid);
    int timestep = 0;
    int output_id = 0;
    double output_counter = 0.0;
    double res;

    /* Get the total number of fluid cells to normalize the residuals */
    int LOCAL_NUMBER_OF_CELLS = _grid.fluid_cells().size();
    int TOTAL_NUMBER_OF_CELLS_REDUCED;
    MPI_Allreduce(&LOCAL_NUMBER_OF_CELLS, &TOTAL_NUMBER_OF_CELLS_REDUCED, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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
        communicate_all(_field.t_matrix(), MessageTag::T);
    }

    // Fluxes (with *new* temperatures)
    _field.calculate_fluxes(_grid);

    // Communicate F and G
    communicate_all(_field.f_matrix(), MessageTag::F);
    communicate_all(_field.g_matrix(), MessageTag::G);

    // Poisson Pressure Equation
    _field.calculate_rs(_grid);
    /* First step of the solver */
    _pressure_solver->init(_field, _grid, _boundaries);

    for (int i = 0; i < N; ++i) {
        res = _pressure_solver->solve(_field, _grid, _boundaries);
            /* Compute TOTAL residual */
            double res_reduction;
            MPI_Allreduce(&res, &res_reduction, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            res = res_reduction / TOTAL_NUMBER_OF_CELLS_REDUCED;
            res = std::sqrt(res);

            // Apply the Boundary conditions (only on pressure)
            for (auto& boundary_ptr : _boundaries) {
                boundary_ptr->apply(_field, true);
            }

            communicate_all(_field.p_matrix(), MessageTag::P);
            out.push_back(res);
    }


    return out;
}