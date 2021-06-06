#include "Case.hpp"
#include "Enums.hpp"
#include "PressureSolver.hpp"
#include"Communication.hpp" 

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <cmath> // for sqrt

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

Case::Case(std::string file_name, int no_of_processors, int my_rank) {
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
            }
        }
    }

    else {
        std::cerr << "Couldn't open file " << file_name << ". Aborting." << std::endl;
        exit(EXIT_FAILURE);
    }

    file.close();

    // Update the boolean for pressure input and heat equation
    if (energy_eq == "on") {
        std::cout << "Enabling heat transfer." << std::endl;
        _use_energy = true;
    } else {
        std::cout << "Heat transfer disabled." << std::endl;
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
    /// Establish Communication : setup the current processor 
    // find closest integer factors of no_of_processors
    
    //-----------------------------------------------------------------------------------------------------------
    int jproc, iproc;
    int sr = std::sqrt(no_of_processors);
    while( no_of_processors % sr != 0 && sr > 0) { sr--; }
    if (!sr) 
        std::runtime_error("no of processors couldn't be divided into sub-domains.");  
    else {
        jproc = sr;
        iproc = sr/jproc;
        std::cerr << "Domain is divided into " << iproc << " by " << jproc << " sub-domains." << std::endl;
    }
    //-----------------------------------------------------------------------------------------------------------
    // Creat an object to represent this processor
    _this_processor   = Processor(iproc, jproc, my_rank);

    // Status update
    std::cerr << "Processor (" << _this_processor.ip() << "," << _this_processor .jp() << ") ready." << std::endl;
    //-----------------------------------------------------------------------------------------------------------
    /// ip,jp are processor (or sub-domain) indices in x and y direction respectively
    int ip = _this_processor .ip(); 
    int jp = _this_processor .jp();

    /// Calculate size of each sub-domain (no of cells)
    /// Calculating the indices of cells in geom_data (global) that correspond to this processor (sub-domain)
    
    int local_size_x,local_size_y;
    int local_igeom_min; 
    int local_jgeom_min; 
    int local_igeom_max; 
    int local_jgeom_max; 
    /// Including cases for single processor run/ 
    // TODO : Optimize , see above , needless computation. 
    // If you want a serial version. 
    if(no_of_processors == 1){
        local_size_x = imax; 
        local_size_y = jmax;
        local_igeom_min = 0;
        local_jgeom_min = 0;
        local_igeom_max = imax+2;  
        local_jgeom_max = jmax+2; 
    }
    // Parallel zone 
    else{
        local_size_x = (imax + 2)/iproc;      // note: imax,jmax is only the interior cells. +2 is required to include the boundary cells.
        local_size_y = (jmax + 2)/jproc
        local_igeom_min = local_size_x * ip;   
        local_jgeom_min = local_size_y * jp;
        local_igeom_max = local_size_x * (ip+1); 
        local_jgeom_max = local_size_y * (jp+1);
    
    /// Give remainder cells to the last sub-domains 
    // (when total no of is cells not divisible by the no processors in each direction)
   
        if(ip == iproc-1)
            local_size_x += (imax+2)%iproc;  
        if(jp == jproc-1)
            local_size_y += (jmax + 2)%jproc;
        }
    }
    
    
    //----------------------------------------------------------------------------------------------------------
    // Build up the domain and include variables for local domain info as well
    Domain domain;
    domain.dx = xlength / (double)imax;
    domain.dy = ylength / (double)jmax;
    
    // I don't see the point of these 
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;  
    

    build_domain(domain,local_size_x, local_size_y,local_igeom_min,local_jgeom_min,local_igeom_max,local_jgeom_max);
    // TODO: also put local_size_x and local_size_y into domain
    //-----------------------------------------------------------------------------------------------------------
    // Load the geometry file
    _grid = Grid(_geom_name, domain,_this_processor ); 

    // Are halo cells only fluid cells. - no, they can be anything. Simply depends on what the corresponding "real" cells in the adjacent sub-domain are.
    //-----------------------------------------------------------------------------------------------------------    
    
    _field = Fields(nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, TI, alpha, beta, GX, GY);
    
    //-----------------------------------------------------------------------------------------------------------
    
    _discretization = Discretization(domain.dx, domain.dy, gamma);
    
    //-----------------------------------------------------------------------------------------------------------
    
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
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

    /* Initialize the logger */
    std::string outputname =
        _dict_name + '/' + _case_name + "_log.txt";

    std::ofstream logger(outputname);
    if (!logger.is_open())
        std::cerr << "Couldn't open the file " << outputname << ". Simulation will run without logs" << std::endl;
    else
        std::cerr << "Starting to log in file " << outputname << std::endl;

    logger << "# iter_number; time ; dt; pressure_iterations; pressure_residual" << std::endl;


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

        _this_processor.communicate(_grid,_field.f_matrix());
        _this_processor.communicate(_grid,_field.g_matrix());
        
        // Poisson Pressure Equation
        _field.calculate_rs(_grid); 

        double rloc,res;
        unsigned iter = 0;
        int comm_fluid_cells = 0; 
        bool flag_next_iteration = true; 
        do {
            rloc = _pressure_solver->solve(_field, _grid, _boundaries);
            // Communicate the pressure data 
            _this_processor.communicate(_grid,_field.p_matrix());
            // Find global res and decide on the next iteration 
            int num_cells = _grid.fluid_cells().size(); 
            MPI_Allreduce(&rloc,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(&num_cells,&comm_fluid_cells,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD); 
            if(my_rank==0){ 
                res = sqrt(res/comm_fluid_cells); 
                if(res<_tolerance)
                    flag_next_iteration = false;
                    MPI_Bcast(&flag_next_iteration,1,MPI_C_BOOL,MPI_COMM_WORLD);
            }
            // Apply the Boundary conditions (only on pressure)
            for (auto& boundary_ptr : _boundaries) {
                boundary_ptr->apply(_field, true);
            }
            ++iter;
        } while (flag_next_iteration && iter < _max_iter);
                
        
        // Update velocity
        _field.calculate_velocities(_grid);

        // Communicate the velocity for the next iteration 
        _this_processor.communicate(_grid,_field.u_matrix()); 
        _this_processor.communicate(_grid,_field.v_matrix());
        // Update logging data and, if enough time has elapsed since the last VTK write ("dt_value" on the .dat file),
        // output the current state.
        logger << timestep << "; " << t << "; " << dt << "; " << iter << "; " << res << std::endl;
        if (output_counter >= _output_freq)
            {
                output_vtk(output_id++);
                output_counter -= _output_freq;
            }


        // Update time elapsed since last VTK write, total time and dt 
        output_counter += dt;
        t += dt;
        timestep += 1;
        
        // Reduce time step : min of all the dt from different processes.  
        dt = _field.calculate_dt(_grid);
        dt = _this_processor.reduce_min(dt);

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

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain, const int& local_imin,
                        const int &local_jmin, const int& local_imax, const int& local_jmax ) {
    

    // How to read this ? 
    // Incase of a parallel implementation : imax_domain (number of ALLsubdomian cells - including boundaries)) 
    // Incase of serial implementation :     imax_domain (number of inner cells - exculding boundaries)                                         

    // Meant to be local 
    domain.imin = 0;
    domain.jmin = 0;
    domain.imax = imax_domain+2;
    domain.jmax = jmax_domain+2;
    
    /// Global indices used to parse the gom file.  
    domain.local_igeom_min = local_imin;
    domain.local_igeom_max = local_imax; 
    domain.local_jgeom_min = local_jmin;
    domain.local_jgeom_max = local_jmax;  
    
    // Local to a specific processor
    domain.size_x = imax_domain;
    domain.size_y = jmax_domain;
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
