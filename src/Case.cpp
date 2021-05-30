#include "Case.hpp"
#include "Enums.hpp"
#include "PressureSolver.hpp"

#include <mpi.h>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>


//  std::string lid_path = "../example_cases/LidDrivenCavity/LidDrivenCavity.pgm";
//-----------------------------------------------------------------------------------------------------------

Case::Case(std::string file_name, int argn, char **args) {
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
            _use_energy = true;
        } else {
            std::cout << "Heat transfer disabled." << std::endl;
        }
    }
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
    _grid = Grid(_geom_name, domain);

//-----------------------------------------------------------------------------------------------------------    
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);
//-----------------------------------------------------------------------------------------------------------
    _discretization = Discretization(domain.dx, domain.dy, gamma);
//-----------------------------------------------------------------------------------------------------------
    _pressure_solver = std::make_unique<SOR>(omg);
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
//-----------------------------------------------------------------------------------------------------------
void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    /* _geom_name handled in constructor */


    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}
//-----------------------------------------------------------------------------------------------------------
/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculat the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
//-----------------------------------------------------------------------------------------------------------

void Case::simulate() {
    double t = 0.0;
    double dt = _field.calculate_dt(_grid);
    int timestep = 0;
    int output_id = 0;
    double output_counter = 0.0;
    std::ofstream logger;

    /* Initialize the logger */
    if (_rank == 0) {
        std::string outputname = _dict_name + '/' + _case_name + "_log.txt";

        logger = std::ofstream(outputname);
        if (!logger.is_open())
            std::cerr << "Couldn't open the file " << outputname << ". Simulation will run without logs" << std::endl;
        else
            std::cerr << "Starting to log in file " << outputname << std::endl;

        logger << "# iter_number; time ; dt; pressure_iterations; pressure_residual" << std::endl;
    }

    /* Main loop */
    while (t < _t_end) {
          
        // Apply the Boundary conditions
        for (auto& boundary_ptr : _boundaries) {
            boundary_ptr->apply(_field);
        }
        
        //Update temperatures
        // (Temperature is still used in calculate_fluxes, but since T is initialized to a constant, it doesn't matter)
        if (_use_energy)
            _field.calculate_T(_grid);
        

        // Fluxes (with *new* temperatures)
        _field.calculate_fluxes(_grid);

        // Poisson Pressure Equation
        _field.calculate_rs(_grid); 

        double res;
        unsigned iter = 0;

        do {
            res = _pressure_solver->solve(_field, _grid, _boundaries);
            // Apply the Boundary conditions (only on pressure)
            for (auto& boundary_ptr : _boundaries) {
                boundary_ptr->apply(_field, true);
            }
            ++iter;
        } while (res > _tolerance && iter < _max_iter);
                
        
        // Update velocity
        _field.calculate_velocities(_grid);

        // Update logging data and, if enough time has elapsed since the last VTK write ("dt_value" on the .dat file),
        // output the current state.
        if (_rank == 0)
            logger << timestep << "; " << t << "; " << dt << "; " << iter << "; " << res << std::endl;
            
        if (output_counter >= _output_freq)
            {
                output_vtk(output_id++, _rank);
                output_counter -= _output_freq;
            }


        // Update time elapsed since last VTK write, total time and dt 
        output_counter += dt;
        t += dt;
        timestep += 1;
        dt = _field.calculate_dt(_grid);



    }

}



void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Temperature Array (Worksheet 2): implemented analogous to pressure array
    vtkDoubleArray *Temperature = vtkDoubleArray::New();
    if (_use_energy) {
        Temperature->SetName("temperature");
        Temperature->SetNumberOfComponents(1);
    }

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Geometry Array (Worksheet 2)
    // TODO: see if there's a IntArray instead of DoubleArray
    vtkDoubleArray *Geometry = vtkDoubleArray::New();
    Geometry->SetName("geometry");
    Geometry->SetNumberOfComponents(1);

    // Print pressure and geometry from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);

            float geometry = (int)_grid.cell(i, j).type(); // worksheet 2
            Geometry->InsertNextTuple(&geometry);          // worksheet 2
            // TODO: For Efficiency, Geometry need not be printed afresh in every time step.
            // Take it out of here and put in Case constructor.
        }
    }

    // Print temperature
    if (_use_energy) {
        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {

                double temperature = _field.t(i, j);        // worksheet 2
                Temperature->InsertNextTuple(&temperature); // worksheet 2
            }
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Temperature to Structured Grid (worksheet 2)
    if (_use_energy)
        structuredGrid->GetCellData()->AddArray(Temperature);

    // Add Geometry to Structured Grid (worksheet 2)
    structuredGrid->GetCellData()->AddArray(Geometry);

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {

    /* Partition the domain in the main process */
    if (_rank == 0) {
        int Lx = imax_domain / _iproc;
        int Ly = jmax_domain / _jproc;
        int rx = imax_domain % _iproc;
        int ry = jmax_domain % _jproc;

        for (int ip = 0; ip < _iproc; ++ip) {
            for (int jp = 0; jp < _jproc; ++jp)
                {
                    // Assign to process ip + iproc * jp the bounds of the reigion (ip, jp)
                    // The will be used to build the geometry and become each case's imin/max and jmin/max on the matrices
                    int target_rank = ip + _iproc * jp;

                   
                        

                    /**
                     * Handling rounding :
                     * Let N be the number of internal cells on a row (i.e. not counting the outer boundaries)
                     * Let P = iproc. Then N = L*P + r with L = floor(N/P) and r = N % P
                     * The r first subdomains will have L+1 true cells (+ ghosts) and the (P-r) lasts will have L cells.
                     * Similar reasoning for vertical boundaries
                     */

                    //0-5 : imin, imax, jmin, jmax, size_x, size_y of the current rank, in terms of the global grid
                    int boundary_data[10]; 
                    
                    boundary_data[0] = ip * Lx + std::min(rx, ip);
                    boundary_data[1] = boundary_data[0] + Lx + 2 + (ip < rx ? 1 : 0); //Add domain_size + 1 to the left boundary
                    boundary_data[2] = jp * Ly + std::min(ry, jp);
                    boundary_data[3] = boundary_data[2] + Ly + 2 + (jp < ry ? 1 : 0);
                    boundary_data[4] = Lx + (ip < rx ? 1 : 0);
                    boundary_data[5] = Ly + (jp < ry ? 1 : 0);

                     // Assign neighbor ranks. Default value of -1 no neighbor (i.e. true border)
                     // 6-9 : left, right, top, bottom
                    if (ip == 0)
                        boundary_data[6] = -1;
                    else
                        boundary_data[6] = target_rank - 1;

                    if (ip == _iproc - 1)
                        boundary_data[7] = -1;
                    else
                        boundary_data[7] = target_rank + 1;

                    if (jp == 0)
                        boundary_data[8] = -1;
                    else
                        boundary_data[8] = target_rank + _iproc;

                    if (jp == _jproc - 1)
                        boundary_data[9] = -1;
                    else
                        boundary_data[9] = target_rank - _iproc;

                    MPI_Send(boundary_data, 10, MPI_INT, target_rank, MessageTag::DOMAIN, MPI_COMM_WORLD);
                }
        }

    }


    // Read local values from the master rank
    int boundary_data[10];
    MPI_Recv(boundary_data, 10, MPI_INT, 0, MessageTag::DOMAIN, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    domain.imin = boundary_data[0];
    domain.imax = boundary_data[1];
    domain.jmin = boundary_data[2];  
    domain.jmax = boundary_data[3];
    domain.size_x = boundary_data[4];
    domain.size_y = boundary_data[5];

    _left_neighbor_rank = boundary_data[6];
    _right_neighbor_rank = boundary_data[7];
    _top_neighbor_rank = boundary_data[8];
    _bottom_neighbor_rank = boundary_data[9];

    // Print the partition
    std::cout << "(" << _rank << ") works on x-domain " << domain.imin << '-' << domain.imax 
    << " and on y-domain " << domain.jmin << '-' << domain.jmax 
    << ". Domain size is " << domain.size_x << 'x' << domain.size_y << std::endl;
}

void Case::setupBoundaryConditions() {
    if (_geom_name.compare("NONE") == 0) {
        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(
                std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }
    } else {

        // General case of BCs

        /**
         *  Temperature BCs :
         *  - If wall_temp_K is -1, we interpret it as an adiabatic wall, otherwise it's a Dirichlet fixed temperature BC.
         *  - Inflow is always Dirichlet
         *  - Outflow is adiabatic (?)
         * */
        if (not _grid.inflow_cells().empty()) {

            if (_use_pressure_input) {
                _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.inflow_cells(), _P_IN));
            }
            else {
                _boundaries.push_back(std::make_unique<InflowBoundary>(_grid.inflow_cells(), _u_in, _v_in));
            }
            _boundaries.push_back(std::make_unique<TemperatureDirichlet>(_grid.inflow_cells(), _T_IN));
        }

        if (not _grid.outflow_cells().empty()) {
            _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.outflow_cells(), _p_i));
            _boundaries.push_back(std::make_unique<TemperatureAdiabatic>(_grid.outflow_cells()));
        }

        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), _wall_velocity));
            if (std::fabs(_moving_wall_temp + 1) < 0.001) //Float comparison
                _boundaries.push_back(std::make_unique<TemperatureAdiabatic>(_grid.moving_wall_cells()));
            else
                _boundaries.push_back(std::make_unique<TemperatureDirichlet>(_grid.moving_wall_cells(), _moving_wall_temp));

        }
        // We don't want to update inner parts of obstacles : check if the cell has 1 or 2 fluid neighbors
        // Wel also split the temperature BCs depending on the wall ID
        std::vector<Cell *> fixed_outer_walls;
        for (auto cell_ptr : _grid.fixed_wall_cells()) {
            if (cell_ptr->borders().size() > 2)
                throw std::runtime_error("A boundary cell has too may fluid cells as neighbors !\n Check input geometry file for Forbidden cells.(Three or more neighbours)");
            else if (cell_ptr->borders().size() > 0)
                fixed_outer_walls.push_back(cell_ptr);
            
        }
        if (not fixed_outer_walls.empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(fixed_outer_walls));
            std::map<cell_type, std::vector<Cell *>> outer_walls_by_ID;
            for (auto cell_ptr : fixed_outer_walls) {
                // An empty vector is default constructed
                outer_walls_by_ID[cell_ptr->type()].push_back(cell_ptr);
            }

            //Create a BC for each (non empty) wall ID
            for (auto& pair : outer_walls_by_ID) {
                cell_type ID = pair.first;
                if (std::fabs(_fixed_wall_temp[ID] + 1) < 0.001)
                    {
                        _boundaries.push_back(std::make_unique<TemperatureAdiabatic>(pair.second));
                    }
                else
                    {
                        _boundaries.push_back(std::make_unique<TemperatureDirichlet>(pair.second, _fixed_wall_temp[ID]));
                    }
            }
        }


    }
}


void Case::communicate_right(Matrix<double>& x, int tag) {
    if (_left_neighbor_rank == -1 && _right_neighbor_rank == -1)
        return; //No horizontal communication

    std::vector<double> out = x.get_col(x.imax()-1);
    std::vector<double> in(x.jmax());

    if (_left_neighbor_rank == -1) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_right_neighbor_rank == -1)
        MPI_Recv(in.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag,
                     in.data(), x.jmax(), MPI_DOUBLE,
                     _left_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the left)
    for (int j = 0; j < x.jmax(); ++j)
        x(0, j) = in[j];
}

void Case::communicate_left(Matrix<double>& x, int tag) {
    if (_left_neighbor_rank == -1 && _right_neighbor_rank == -1)
        return; //No horizontal communication

    std::vector<double> out = x.get_col(0);
    std::vector<double> in(x.jmax());

    if (_right_neighbor_rank == -1) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_left_neighbor_rank == -1)
        MPI_Recv(in.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag,
                     in.data(), x.jmax(), MPI_DOUBLE,
                     _right_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the right)
    auto i = x.imax()-1;
    for (int j = 0; j < x.jmax(); ++j)
        x(i, j) = in[j];
}

void Case::communicate_top(Matrix<double>& x, int tag) {
    if (_top_neighbor_rank == -1 && _bottom_neighbor_rank == -1)
        return; //No vertical communication

    std::vector<double> out = x.get_row(x.jmax()-1);
    std::vector<double> in(x.imax());

    if (_bottom_neighbor_rank == -1) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_top_neighbor_rank == -1)
        MPI_Recv(in.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag,
                     in.data(), x.imax(), MPI_DOUBLE,
                     _bottom_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the bottom)
    for (int i = 0; i < x.imax(); ++i)
        x(i, 0) = in[i];
}

void Case::communicate_bottom(Matrix<double>& x, int tag) {
    if (_top_neighbor_rank == -1 && _bottom_neighbor_rank == -1)
        return; //No vertical communication

    std::vector<double> out = x.get_row(0);
    std::vector<double> in(x.imax());

    if (_top_neighbor_rank == -1) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_top_neighbor_rank == -1)
        MPI_Recv(in.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag,
                     in.data(), x.imax(), MPI_DOUBLE,
                     _top_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the top)
    auto j = x.jmax() - 1;
    for (int i = 0; i < x.imax(); ++i)
        x(i, j) = in[i];
}