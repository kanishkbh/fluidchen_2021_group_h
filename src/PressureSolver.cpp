#include "PressureSolver.hpp"

#include <cmath>
#include <iostream>

SOR::SOR(double omega) : _omega(omega) {}

double SOR::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {

    double imax = grid.imax();
    double jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();

    double coeff = _omega / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy))); // = _omega * h^2 / 4.0, if dx == dy == h

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) = (1.0 - _omega) * field.p(i, j) +
                        coeff * (Discretization::sor_helper(field.p_matrix(), i, j) - field.rs(i, j));
    }

    double res = 0.0;
    double rloc = 0.0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
        rloc += (val * val);
    }
    {
        // Compute the sum of errors, and normalize in the main loop.
        // This is because each process needs to know the total number of cells to divide the total residual

        // res = rloc / (grid.fluid_cells().size());
        // res = std::sqrt(res);
    }

    

    
    return rloc;
}

double CG::init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    // 1) Initialization : Compute residual matrix. Direction is residual at first
    double imax = grid.imax();
    double jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();


    residual = Matrix<double>(field.p_matrix().imax(), field.p_matrix().jmax());
    // Reset pressure (Neumann is trivially true then) then apply Dirichlet BC if any
    for (int i =0; i < residual.imax(); ++i) {
        for (int j = 0; j < residual.jmax(); ++j) {
            field.p(i, j) = 0;
        }
    }

    for (auto& boundary_ptr : boundaries) {
            boundary_ptr->apply(field, true);
        }

    /* Compute residual = "b - Ax" but x is 0 in fluid cells => b 
       Search direction is initialized to residual too
    */
    residual = field.rs_matrix();
    direction = residual;
    a_direction = direction;

    /* Square residual of the previous step is reused. Let's compute it for frist iteration */
    /* Also needeed is the product A*d */
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = field.rs(i, j);
        square_residual += (val * val);
        a_direction(i, j) = Discretization::laplacian(direction, i, j);
    }


}

double CG::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    
    /*Compute alpha =  square_residual / (d*A*d) */
    // TODO : ADD COMMUNICATION

    double alpha = square_residual;
    double dAd = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = direction(i, j) * a_direction(i, j);
        dAd += val;
    }

    alpha /= dAd;
    std::cout << "Alpha : " << alpha << std::endl;

    /* Update x and r. Also compute new square residual */
    double new_square_res = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) += alpha * direction(i, j);
        residual(i, j) -= alpha * a_direction(i, j);

        new_square_res += (residual(i, j) * residual(i, j));
        
    }

    double beta = new_square_res / square_residual;
    square_residual = new_square_res;

    /* Compute new direction, then A*direction */

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = residual(i, j) + beta*direction(i, j);
        direction(i, j) = val;
    }

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        a_direction(i, j) = Discretization::laplacian(direction, i, j);
    }

    return square_residual;
}