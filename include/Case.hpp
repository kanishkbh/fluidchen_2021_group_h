#pragma once

#include <memory>
#include <string>
#include <vector>

#include "Boundary.hpp"
#include "Discretization.hpp"
#include "Domain.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"
#include "Datastructures.hpp"

/**
 * @brief Class to hold and orchestrate the simulation flow.
 *
 */
class Case {
  public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    Case(std::string file_name, int no_of_processors, int my_rank);

    /**
     * @brief Main function to simulate the flow until the end time.
     *
     * Calculates the fluxes
     * Calculates the right hand side
     * Solves pressure
     * Calculates velocities
     * Outputs the solution files
     */
    void simulate();

  private:
    /// Plain case name without paths
    std::string _case_name;
    /// Output directiory name
    std::string _dict_name;
    /// Geometry file name
    std::string _geom_name{"NONE"};
    /// Relative input file path
    std::string _prefix;

    /// Simulation time
    double _t_end;
    /// Solution file outputting frequency
    double _output_freq;

    Processor _this_processor; 
    Fields _field;
    Grid _grid;
    Discretization _discretization;
    std::unique_ptr<PressureSolver> _pressure_solver;
    std::vector<std::unique_ptr<Boundary>> _boundaries;

    // Boundary conditions parameters
    double _wall_velocity, _u_in, _v_in, _p_i;
    // Temperature BC
    double _T_IN{0}, _moving_wall_temp{0};
    std::map<cell_type, double> _fixed_wall_temp;
    // Inflow parameters
    bool _use_pressure_input = false;
    double _P_IN{0};
    // Use energy ?
    bool _use_energy = false;

    /// Solver convergence tolerance
    double _tolerance;
    

    /// Maximum number of iterations for the solver
    int _max_iter;

    /**
     * @brief Creating file names from given input data file
     *
     * Extracts path of the case file and creates code-readable file names
     * for outputting directory and geometry file.
     *
     * @param[in] input data file
     */
    void set_file_names(std::string file_name);

    /**
     * @brief Solution file outputter
     *
     * Outputs the solution files in .vtk format. Ghost cells are excluded.
     * Pressure is cell variable while velocity is point variable while being
     * interpolated to the cell faces
     *
     * @param[in] Timestep of the solution
     */
    void output_vtk(int t, int my_rank = 0);

    /** 
     *  @brief Solution log outputter
     * 
     *  Outputs simulation data in a text file, containing the history of pressure iterations and time stepping.
     */
    void output_simulation_logs(const std::vector<int>& pressure_iter, const std::vector<double>& dts);

    void build_domain(Domain &domain, int imax_domain, int jmax_domain,const int& local_imin,const int& local_jmin,
                      const int& local_imax,const int& local_jmax);

    // Utility function to fill _boundaries
    void setupBoundaryConditions();
};
