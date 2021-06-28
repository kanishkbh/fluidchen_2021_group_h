#include <iostream>
#include <string>
#include "include/Discretization.hpp"
#include "include/Fields.hpp"

// Function declaration (definition after main)
void set_A_matrix(Matrix<double>*, int, int, int, double, double);


// MAIN
// - set up Fields and PressureSolver objects
// - set up A-matrix
// - set up RHS-vector according to F and G
// - add known pressure values (boundaries of field.p matrix) to RHS vector
// - solve using a standard solver
// - solve using our own iterative solver
// - compare and print pass/fail
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
        double alpha{0};     /* Heat diffusivity */
        double PI = 0;           /* pressure */
        double GX = 0;           /* gravitation x-direction */
        double GY = 0;           /* gravitation y-direction */
        double xlength;      /* length of the domain x-dir.*/
        double ylength;      /* length of the domain y-dir.*/
        double dt = 0.05;           /* time step */
        int imax = 5;            /* number of cells x-direction*/
        int jmax = 7;            /* number of cells y-direction*/
        double gamma;        /* uppwind differencing factor*/
        double omg;          /* relaxation factor */
        double tau = 0.5;          /* safety factor for time step*/
        int itermax;         /* max. number of iterations for pressure per time step */
        double eps;          /* accuracy bound for pressure*/
        double UIN;          /* Inlet horizontal velocity */
        double VIN;          /* Inlet vertical velocity */
        double in_temp;      /* Inlet (Dirichlet) Temperature */
        
        Fields field(nu, dt, tau, imax, jmax, UI, VI, PI, TI, alpha, beta, GX, GY);
        // TODO:    complete and check fields initialization
        //          esp. check that field.p size is imax+1, jmax+1

        double dx = 0.5, dy = 1; // TODO: set realistic values for dx,dy
        int rmax = imax*jmax;
        int cmax = imax*jmax;
        Matrix<double> A(rmax, cmax, 0);
        set_A_matrix(&A, rmax, cmax, imax, dx, dy);

        // RHS vector
        std::vector<double> rhs_vec;
        int k=0; // iterator for rhs_vec
        for (int j=1; j<jmax-1; j++) {
            for(int i=1; i<imax-1; i++) {
                rhs_vec[k] = field.rs(i, j); // TODO: declare and define field object
                k++;
            }
        }
        // add known (boundary) p values in rhs_vec
        // [TODO]

        // TODO: call a direct solver from some library

        // TODO: declare and call a pressure solver

        // TODO: compare solutions of library solver and our solver, print test pass/fail

    // } else {
    //     std::cout << "Error: No input file is provided." << std::endl;
    //     // std::cout << "Example usage: /path/to/fluidchen /path/to/input_data.dat" << std::endl;
    // }
}

// FUnction definition
void set_A_matrix(Matrix<double> *A, int rmax, int cmax, int imax, double dx, double dy) {

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

    // print matrix
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