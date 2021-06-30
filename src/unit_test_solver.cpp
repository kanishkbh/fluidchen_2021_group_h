#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <memory>

#include <eigen3/Eigen/Dense>
#include "Discretization.hpp"
#include "Fields.hpp"
#include "Domain.hpp"
#include "Grid.hpp"
#include "PressureSolver.hpp"
#include "Boundary.hpp"

// Function declarations (definitions after main)
void set_A_matrix(Matrix<double>*, int, int, int, double, double);
void set_A_matrix(Eigen::MatrixXd*, int, int, int, double, double);
void set_Boundary_Conditions(std::vector<std::unique_ptr<Boundary>>& _boundaries, Grid& _grid);

/*
MAIN
- set up Fields and PressureSolver objects
- set up A-matrix
- set up RHS-vector according to F and G
- add known pressure values (boundaries of field.p matrix) to RHS vector
- solve using a standard solver
- solve using our own iterative solver
- compare and print pass/fail
*/
int main(int argn, char **args) {

    // if (argn > 1) {
        // std::string file_name{args[1]};
        // Case problem(file_name, argn, args);
        // problem.simulate();
        
        // dummy values
        double nu = 0.01;           /* viscosity   */
        double UI = 0;           /* velocity x-direction */
        double VI = 0;           /* velocity y-direction */
        double TI = 0;           /* Initial temperature */
        double beta = 0;         /* Thermal expansion coefficient */
        double alpha{0};            /* Heat diffusivity */
        double PI = 0;              /* pressure */
        double GX = 0;              /* gravitation x-direction */
        double GY = 0;              /* gravitation y-direction */
        double xlength;             /* length of the domain x-dir.*/
        double ylength;             /* length of the domain y-dir.*/
        double dt = 0.05;           /* time step */
        int imax = 5;               /* number of cells x-direction*/
        int jmax = 7;               /* number of cells y-direction*/
        double gamma = 0.5;        /* uppwind differencing factor*/
        double omg = 1.7;          /* relaxation factor */
        double tau = 0.5;          /* safety factor for time step*/
        int itermax = 500;         /* max. number of iterations for pressure per time step */
        double eps = 0.001;        /* accuracy bound for pressure*/
        double UIN;                 /* Inlet horizontal velocity */
        double VIN;                 /* Inlet vertical velocity */
        double in_temp;             /* Inlet (Dirichlet) Temperature */
        double dx = 0.5, dy = 1;    /// TODO: set realistic values for dx,dy
        std::string geom_name = "NONE";


        // opening files
        std::string pressure_file_path = "../test_matrices/pressure_matrix.txt";
        std::ifstream pressure_file;
        pressure_file.open(pressure_file_path);
        std::string x_velocity_file_path = "../test_matrices/x_velocity_matrix.txt";
        std::ifstream x_velocity_file(x_velocity_file_path);
        std::string y_velocity_file_path = "../test_matrices/y_velocity_matrix.txt";
        std::ifstream y_velocity_file(y_velocity_file_path);

        // asserts
        assert(pressure_file.is_open());    std::cout << "opened "<< pressure_file_path << std::endl;
        assert(x_velocity_file.is_open());  std::cout << "opened "<< x_velocity_file_path << std::endl;
        assert(y_velocity_file.is_open());  std::cout << "opened "<< y_velocity_file_path << std::endl;

        // READ FROM FILES
        std::stringstream ss, ss2, ss3;
        ss << pressure_file.rdbuf();
        ss2 << x_velocity_file.rdbuf();
        ss3 << y_velocity_file.rdbuf();
        // First line : size
        ss >> imax >> jmax;
        ss2 >> imax >> jmax; // reading to advance the buffer to next line
        ss3 >> imax >> jmax; // reading to advance the buffer to next line

        //-------------------------------------------------------------------------------------------

        // INITIALIZATIONS
        Discretization discretization(dx, dy, gamma);
        Fields field(nu, dt, tau, imax, jmax, UI, VI, PI, TI, alpha, beta, GX, GY);
        std::cout << "initialized discretization and fields.\n";
        std::vector<std::vector<double>> x_velocity_matrix(imax-2, std::vector<double>(jmax-2, 0)); // not used
        std::vector<std::vector<double>> y_velocity_matrix(imax-2, std::vector<double>(jmax-2, 0)); // not used
        std::vector<std::vector<double>> pressure_matrix  (imax-2, std::vector<double>(jmax-2, 0)); // not used
        
        Domain domain;
        domain.imax = imax;
        domain.jmax = jmax;
        domain.imin = 0;
        domain.jmin = 0;
        domain.size_x = imax-2;
        domain.size_y = jmax-2;
        domain.domain_size_x = imax-2;
        domain.domain_size_y = jmax-2;
        domain.dx = dx; domain.dy = dy;
        std::cout << "initialized domain.\n";

        Grid grid(geom_name, domain, false, false, false, false);
        std::cout << "initialized grid.\n";
        std::cout << "grid.imax,jmax = " << grid.imax() << "," << grid.jmax() << std::endl;

        //-------------------------------------------------------------------------------------------

        // Read p,u,v matrices from respective files
        std::cout << "\nreading the files into matrices and fields...\n";
        for (int col = 0; col < imax; ++col) {
            for (int row = 0; row < jmax; ++row) {

                double p; ss >> p;
                double u; ss2 >> u;
                double v; ss3 >> v;
                
                field.p(row,col) = p;
                field.u(row,col) = u;
                field.v(row,col) = v;
                
                if(col > 0 && row > 0 && col < imax-1 && row < jmax-1) { // collecting only interior points
                    pressure_matrix[col-1][row-1] = p;    // not used anywhere
                    x_velocity_matrix[col-1][row-1] = u;  // not used anywhere
                    y_velocity_matrix[col-1][row-1] = p;  // not used anywhere
                }
                else { // Set boundaries of fields to zero (for now)
                    field.p(row,col) = 0.0; 
                    field.u(row,col) = 0.0;
                    field.v(row,col) = 0.0;
                }
            }
        }
        std::cout<< "matrices and fields initialization complete.\n";

        // PRINTS FOR DEBUGGING
        // print pressure matrix
        // std::cout << "pressure_matrix" << std::endl;
        // std::cout << "-------------------------------" << std::endl;
        // std::cout << std::setiosflags(std::ios::fixed);
        // std::cout << std::setprecision(3);
        // for(int r=0;r<jmax-2;r++) {
        //     for(int c=0;c<imax-2;c++) {
        //         std::cout << pressure_matrix[r][c] << ' ';
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout<<"----------------------"<<std::endl;

        // print field.p
        std::cout << "field.p: " <<imax << "x" << jmax << std::endl;
        std::cout << "-------------------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        for(int j=0;j<jmax;j++) {
            for(int i=0;i<imax;i++) {
                std::cout << field.p(i,j) << ' ';
            }
            std::cout << std::endl;
        }
        std::cout<<"----------------------"<<std::endl;

        // // print x-velocity matrix
        // std::cout << "x_velocity_matrix" << std::endl;
        // std::cout << "-------------------------------" << std::endl;
        // std::cout << std::setiosflags(std::ios::fixed);
        // std::cout << std::setprecision(3);
        // for(int r=0;r<jmax-2;r++) {
        //     for(int c=0;c<imax-2;c++) {
        //         std::cout << x_velocity_matrix[r][c] << ' ';
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout<<"----------------------"<<std::endl;

        // print field.u
        std::cout << "field.u: " << imax << "x" << jmax << std::endl;
        std::cout << "-------------------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        for(int j=0;j<jmax;j++) {
            for(int i=0;i<imax;i++) {
                std::cout << field.u(i,j) << ' ';
            }
            std::cout << std::endl;
        }
        std::cout<<"----------------------"<<std::endl;

        // // print y-velocity matrix
        // std::cout << "y_velocity_matrix" << std::endl;
        // std::cout << "-------------------------------" << std::endl;
        // std::cout << std::setiosflags(std::ios::fixed);
        // std::cout << std::setprecision(3);
        // for(int r=0;r<jmax-2;r++) {
        //     for(int c=0;c<imax-2;c++) {
        //         std::cout << y_velocity_matrix[r][c] << ' ';
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout<<"----------------------"<<std::endl;

        // print field.v
        std::cout << "field.v: " << imax << "x" << jmax << std::endl;
        std::cout << "-------------------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        for(int j=0;j<jmax;j++) {
            for(int i=0;i<imax;i++) {
                std::cout << field.v(i,j) << ' ';
            }
            std::cout << std::endl;
        }
        std::cout<<"----------------------"<<std::endl;
        
        //-------------------------------------------------------------------------------------------

       // SETUP LHS "A-MATRIX" FOR THE STANDARD MATRIX SOLVER (using 1/1/-4/1/1 stencil)
        const int rmax = (imax-2)*(jmax-2);
        const int cmax = (imax-2)*(jmax-2);
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(rmax, cmax);   
        set_A_matrix(&A, rmax, cmax, imax-2, dx, dy);

        // Matrix<double> A(rmax, cmax, 0);                 // left here for debugging
        // set_A_matrix(&A, rmax, cmax, imax-2, dx, dy);    // left here for debugging
        // Matrix<double> A(25, 25, 0);                     // left here for debugging
        // set_A_matrix(&A, 25, 25, 5, dx, dy);             // left here for debugging
        // Eigen::MatrixXd A = Eigen::MatrixXd::Zero(25, 25);// left here for debugging
        // set_A_matrix(&A, 25, 25, 5, dx, dy);              // left here for debugging

        //-------------------------------------------------------------------------------------------

        std::cout << "calling field.calculate_fluxes...\n";
        field.calculate_fluxes(grid);
        std::cout << "calling field.calculate_rs...\n";
        field.calculate_rs(grid);

        // // print field.f (for debugging)
        // std::cout << "\nfield.f" << std::endl;
        // std::cout << "-------------------------------" << std::endl;
        // std::cout << std::setiosflags(std::ios::fixed);
        // std::cout << std::setprecision(3);
        // for(int j=0;j<jmax;j++) {
        //     for(int i=0;i<imax;i++) {
        //         std::cout << field.f(i,j) << ' ';
        //     }
        //     std::cout << std::endl;
        // }
        // std::cout<<"----------------------"<<std::endl;

        //-------------------------------------------------------------------------------------------

        // SETUP RHS VECTOR FOR THE STANDARD MATRIX SOLVER (just unwrap the interior values of fields.rs)
        // std::vector<double> rhs_vec((imax-2)*(jmax-2), 0);
        std::cout << "Setting up rhs_vec (Eigen::VectorXd)...\n";
        Eigen::VectorXd rhs_vec = Eigen::VectorXd::Zero((imax-2)*(jmax-2)) ;
        int k=0; // iterator for rhs_vec
        for (int j=1; j<jmax-1; j++) {
            for(int i=1; i<imax-1; i++) {
                rhs_vec[k] = field.rs(i, j); // TODO: declare and define field object
                k++;
            }
        }
        std::cout << "rhs_vec set-up complete.\n";

        
        // FOR DEBUGGING
        // print field.rs
        std::cout << "\nfield.rs: " << imax << "x" << jmax << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        for(int j=0;j<jmax;j++) {
            for(int i=0;i<imax;i++) {
                std::cout << field.rs(i,j) << ' ';
            }
            std::cout << std::endl;
        }
        std::cout << "-------------------------------" << std::endl;

        // // print rhs_vec
        std::cout << "\nrhs_vector (for standard solver): ("<< rhs_vec.size() << ")" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        // for(int k=0; k < (imax-2)*(jmax-2); k++) {
        //     std::cout << rhs_vec[k] << ' ';
        // }
        std::cout << rhs_vec.transpose() << std::endl;
        std::cout << "-------------------------------" << std::endl;

        //-------------------------------------------------------------------------------------------
        
        /// TODO: add known (boundary) p values in rhs_vec (using zero boundary values for now)

        //-------------------------------------------------------------------------------------------

        // Call a direct solver from Eigen library
        Eigen::VectorXd solution_vec = A.colPivHouseholderQr().solve(rhs_vec);
        std::cout << "\nSolution using Eigen::colPivHouseholderQr: ("<< solution_vec.size() << ")" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        std::cout << solution_vec.transpose() << std::endl;
        std::cout << "----------------------" << std::endl;

        //-------------------------------------------------------------------------------------------

        /// Initialize boundaries (PressureSolver::solve() needs it but doesn't use it)
        std::vector<std::unique_ptr<Boundary>> boundaries;
        assert(geom_name.compare("NONE") == 0);
        set_Boundary_Conditions(boundaries, grid);
        std::cout << "boundaries set-up complete.\n";

        // Initialize and Call a pressure solver
        std::unique_ptr<PressureSolver> pressure_solver = std::make_unique<SOR>(omg);
        // pressure_solver->init(field, grid, boundaries);  // SOR doesn't have init()
        double res;
        unsigned iter = 0;
        auto _tolerance = eps;
        auto _max_iter = itermax;
        do {
            res = pressure_solver->solve(field, grid, boundaries);
            // Apply the Boundary conditions (only on pressure)
            // for (auto& boundary_ptr : _boundaries) {
            //     boundary_ptr->apply(_field, true);
            // }
            ++iter;
        } while (res > _tolerance && iter < _max_iter);
        std::cout << "Pressure solver complete." << std::endl;

        // SETUP SOLUTION VECTOR FROM field.p (just unwrap the interior values of fields.p)
        Eigen::VectorXd p_vec = Eigen::VectorXd::Zero((imax-2)*(jmax-2)) ;
        k=0; // iterator for p_vec
        for (int j=1; j<jmax-1; j++) {
            for(int i=1; i<imax-1; i++) {
                p_vec[k] = field.p(i, j); // TODO: declare and define field object
                k++;
            }
        }
        std::cout << "\np_vector using SOR: (" << p_vec.size() << ")" << std::endl;
        std::cout << "----------------------" << std::endl;
        std::cout << std::setiosflags(std::ios::fixed);
        std::cout << std::setprecision(3);
        std::cout << p_vec.transpose() << std::endl;
        std::cout << "-------------------------------" << std::endl;

        //-------------------------------------------------------------------------------------------


    // } else {
    //     std::cout << "Error: No input file is provided." << std::endl;
    //     // std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    // }
}

//-------------------------------------------------------------------------------------------

// NOT USED. Function definition (overload 1: use Datastructures/Matrix<double>)
void set_A_matrix(Matrix<double> *A, int rmax, int cmax, int imax, double dx, double dy) {

    std::cout << "\nSetting up A-Matrix, size: " << rmax << "x" << cmax << "...\n";
    int r,c; // iterators

    // diagonal elements
    for(c=0,r=0; r<rmax; r++,c++) {
        (*A)(c,r) = -4;
    }

    // off diagonal upper
    for(r=0, c=imax; c < cmax; c++,r++) {
        (*A)(c,r) = 1 / (dy*dy);
        // std::cout << "c = " << c << ", r = " << r << std::endl;
    }

    // off diagonal lower
    for(c=0, r=imax; r<rmax; c++,r++) {
        (*A)(c,r) =1/(dy*dy);
    }

    // superdiagonal
    for(r=0, c=1; c<cmax; c++,r++) {
        if (c % imax != 0)
            (*A)(c,r) = 1/(dx*dx);
    }

    // sub diagonal
    for(r=1,c=0; r<rmax; r++,c++) {
        if (r % imax != 0)
            (*A)(c,r) = 1/(dx*dx);
    }
    std::cout << "A-matrix set-up complete.\n"<<std::endl;

    
    // print matrix (for debugging)
    std::cout << "\nA-Matrix:\n";
    auto f = std::cout.flags();
    std::cout << "-------------------------------" << std::endl;
    std::cout << std::setiosflags(std::ios::fixed);
    std::cout << std::setprecision(1);
    for(r=0;r<rmax;r++) {
        for(c=0;c<cmax;c++) {
            std::cout << (*A)(c,r) << ' ';
        }
        std::cout << std::endl;
    }
    std::cout<<"----------------------"<<std::endl;
    std::cout.flags(f);

}

//-------------------------------------------------------------------------------------------

// Function definition (overload 2: uses Eigen::matrix)
void set_A_matrix(Eigen::MatrixXd *A, int rmax, int cmax, int imax, double dx, double dy) {

    std::cout << "\nSetting up A-Matrix (Eigen::MatrixXd), size: " << rmax << "x" << cmax << "...\n";
    int r,c; // iterators

    // diagonal elements
    for(c=0,r=0; r<rmax; r++,c++) {
        (*A)(r,c) = -4;
    }

    // off diagonal upper
    for(r=0, c=imax; c < cmax; c++,r++) {
        (*A)(r,c) = 1 / (dy*dy);
        // std::cout << "c = " << c << ", r = " << r << std::endl;
    }

    // off diagonal lower
    for(c=0, r=imax; r<rmax; c++,r++) {
        (*A)(r,c) =1/(dy*dy);
    }

    // superdiagonal
    for(r=0, c=1; c<cmax; c++,r++) {
        if (c % imax != 0)
            (*A)(r,c) = 1/(dx*dx);
    }

    // sub diagonal
    for(r=1,c=0; r<rmax; r++,c++) {
        if (r % imax != 0)
            (*A)(r,c) = 1/(dx*dx);
    }
    std::cout << "A-matrix set-up complete.\n"<<std::endl;

    /*
    // print matrix (for debugging)
    std::cout << "\nEigen::A-Matrix:\n";
    std::cout << "-----------------------" << std::endl;
    std::cout << std::setiosflags(std::ios::fixed);
    std::cout << std::setprecision(1);
    std::cout << (*A) << std::endl;
    std::cout << "-------------------------" << std::endl;
    */
}

//-------------------------------------------------------------------------------------------

void set_Boundary_Conditions(std::vector<std::unique_ptr<Boundary>>& _boundaries, Grid& _grid) {
     
        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
                std::cout << "added moving walls to boundaries\n";
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
            std::cout << "added fixed walls to boundaries\n";
        }
    
}
//-------------------------------------------------------------------------------------------
