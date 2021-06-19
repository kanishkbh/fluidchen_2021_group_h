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

    /* Square residual of the previous step is reused. Let's compute it for frist iteration */
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = field.rs(i, j);
        square_residual += (val * val);
    }


}

double CG::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    // 1) Initialization : Compute residual matrix. Direction is residual at first
    // Loop :
    /**
     * Update with alpha = r²/(dAd) (x = x + alpha*d)
     * r = r - alpha*A*d
     * beta = new_r²/old_r²
     * new direction : new_r + beta*d_i
     * */

}