#include "PressureSolver.hpp"
#include "Communication.hpp"

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

void CG::init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    // 1) Initialization : Compute residual matrix. Direction is residual at first
    double imax = grid.imax();
    double jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();


    residual = Matrix<double>(field.p_matrix().imax(), field.p_matrix().jmax());

    /* Compute residual = "b - Ax" but x is 0 in fluid cells => b 
       Search direction is initialized to residual too
    */
    residual = field.rs_matrix();

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        // Residual = b-Ax
        residual(i, j) -= Discretization::laplacian(field.p_matrix(), i, j);
    }

    direction = residual;
    a_direction = direction;

    /* Square residual of the previous step is reused. Let's compute it for frist iteration */
    /* Also needeed is the product A*d */
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = residual(i, j);
        square_residual += (val * val);
        a_direction(i, j) = Discretization::laplacian(direction, i, j);

    }

    double total_residual;
    //Dot product over all of the grid
    
    Communication::communicate_sum_double(&square_residual, &total_residual);
    square_residual = total_residual;


}

double CG::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    
    /*Compute alpha =  square_residual / (d*A*d) */

    double alpha = square_residual;
    double dAd = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        double val = direction(i, j) * a_direction(i, j);
        dAd += val;
    }

    /* Dot product over all the grid */
    double sum_dAd = dAd;
    Communication::communicate_sum_double(&dAd, &sum_dAd);
 
    alpha /= sum_dAd;

    /* Update x and r. Also compute new square residual */
    double new_square_res = 0;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) += alpha * direction(i, j);
        
    }

    /* Critical : communicate pressure ! */
    Communication::communicate_all(field.p_matrix(), MessageTag::P);
    Communication::communicate_all(direction, MessageTag::P);

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        residual(i, j) -= alpha * a_direction(i, j);
        new_square_res += (residual(i, j) * residual(i, j));
        
    }

    double total_residual;
    double local_res_to_return = new_square_res;
    //Dot product over all of the grid
    Communication::communicate_sum_double(&new_square_res, &total_residual);
    new_square_res = total_residual;

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

    return local_res_to_return;
}

void SD::init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    a_residual = Matrix<double>(field.p_matrix().imax(), field.p_matrix().jmax());
}

double SD::solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
    
    //Compute residual, A*residuall. Descent : x = x + alpha*res
    //where alpha = r.r/rAr

    residual = field.rs_matrix();
    double square_res{0}, rAr{0};

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        residual(i, j) -= Discretization::laplacian(field.p_matrix(), i, j);
        square_res += (residual(i, j) * residual(i, j));
        
    }

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        a_residual(i, j) = Discretization::laplacian(residual, i, j);
        rAr += (residual(i, j) * a_residual(i, j));
    }

    //TODO : sum
    double total_rr, total_rar;
    Communication::communicate_sum_double(&square_res, &total_rr);
    Communication::communicate_sum_double(&rAr, &total_rar);
    double alpha = total_rr/total_rar;

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        field.p(i, j) += alpha * residual(i,j);
    }


    return square_res;
}