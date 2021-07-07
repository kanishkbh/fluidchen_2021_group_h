#include "PressureSolver.hpp"
#include "Communication.hpp"
// #define DEBUG 
#include <cmath>
#include <iostream>
static int iter = 1; 
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
        a_direction(i, j) = Discretization::boundary_aware_laplacian(direction, i, j, grid);

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

    /* Update border pressure */
    for (auto& boundary_ptr : boundaries) {
                boundary_ptr->apply(field, true);
            }

    /* Critical : communicate pressure ! */
    Communication::communicate_all(field.p_matrix(), MessageTag::P);


    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        residual(i, j) -= alpha * a_direction(i, j);      
        new_square_res += (residual(i, j) * residual(i, j));
        
    }

    Communication::communicate_all(residual, MessageTag::P);

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


    Communication::communicate_all(direction, MessageTag::P);

    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        a_direction(i, j) = Discretization::boundary_aware_laplacian(direction, i, j, grid);
    }

     Communication::communicate_all(a_direction, MessageTag::P);

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
    square_res = 0;
    residual = field.rs_matrix();
    for (auto currentCell : grid.fluid_cells()) {
        int i = currentCell->i();
        int j = currentCell->j();

        residual(i, j) -= Discretization::laplacian(field.p_matrix(), i, j);
        square_res += (residual(i, j) * residual(i, j));
        
    }


    return square_res;
}


void CG_Richardson::init(Fields& field, Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries)
{
    // Initialize residual and direction 
    int imax = grid.imax();
    int jmax = grid.jmax(); 
    double dx = grid.dx(); 
    double dy = grid.dy(); 

    residual = Matrix<double>(field.p_matrix().imax(),field.p_matrix().jmax()); //
    // Ax = b 
    // b is now in a matrix form 
    // r = b - Ax 
    // x = 0 matrix 
    // r0 = b;
    residual = field.rs_matrix(); 
    // In PPE, x aka p is no longer zero 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i(); 
        int j = cell->j(); 
        residual(i,j) = residual(i,j) -  Discretization::laplacian(field.p_matrix(),i,j); 
    }
    direction = residual; 
    adirection = direction; 

    // Parameters to compute alpha 
    // residual(it-1)
    // direction(it)
    // adirection(it)
    // \alpha = r(i,j)*r(i,j) / (direction(i,j) * adirection(i,j)); 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        double val = residual(i,j); 
        square_residual += (val*val);
        adirection(i,j) = Discretization::laplacian(direction,i,j);  
    }
    double total_residual;
    Communication::communicate_sum_double(&square_residual,&total_residual); 
    square_residual = total_residual; 
}

double CG_Richardson::solve(Fields& field,Grid& grid, const std::vector<std::unique_ptr<Boundary>>& boundaries)
{
#ifdef DEBUG
    iter ++;
#endif
    // Compute alpha 
        /// Compute d*Ad
    double alpha = square_residual;  
#ifdef DEBUG 
    if (iter%100 == 0) 
        std::cerr << "denominator " << square_residual<< std::endl; 
#endif 
    double dAd = 0; 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        double val = direction(i,j) * adirection(i,j);  
        dAd += val; 
    }
    double total_denom = dAd; 
#ifdef DEBUG 
    if (iter%100 == 0) 
        std::cerr << "denominator " << total_denom << std::endl; 
#endif 
    Communication::communicate_sum_double(&dAd,&total_denom); 
    alpha = alpha/total_denom; 
#ifdef DEBUG
    if(iter % 100 == 0)
        std::cerr << "alpha " << alpha << std::endl;
#endif 
    // Update x and residual
    double new_residual_squared = 0; 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        field.p(i,j) = field.p(i,j) +  alpha*direction(i,j); 
    }
    // Communicate to all processes 
    Communication::communicate_all(field.p_matrix(),MessageTag::P); 
    Communication::communicate_all(direction,MessageTag::P);

    // Update residuals 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= alpha * adirection(i,j); 
        new_residual_squared += (residual(i,j)*residual(i,j)); 
    }

    double total_residual=0;
    double return_val = new_residual_squared;
    Communication::communicate_sum_double(&new_residual_squared,&total_residual); 
    new_residual_squared = total_residual; 
    // Compute beta 
    double beta = new_residual_squared / square_residual; 
    square_residual = new_residual_squared;
    // Compute direction
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        double val = residual(i,j) + beta*direction(i,j); 
        direction(i,j) = val;  
    }
    // 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        adirection(i,j) = Discretization::laplacian(direction,i,j); 
    }
    return return_val;
}

void CG_Jacobi::init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries)
{
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy(); 

    residual = Matrix<double>(imax+2,jmax+2); 
    cond_residual = residual; 
    residual = field.rs_matrix();
    // As we know the pressure is non homogenous 
    
    // r0 = b - A p 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= Discretization::laplacian(field.p_matrix(),i,j);  
    }
    // q1 = P^-1 r 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        cond_residual(i,j) = Discretization::jacobi(residual,i,j); 
    }
    // \sigma_0 = <r0,q1> 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        sigma += residual(i,j) * cond_residual(i,j);  
    }
    direction = cond_residual; 
    adirection = direction; 
    double total_sigma = 0; 
    Communication::communicate_sum_double(&sigma,&total_sigma); 
    sigma = total_sigma; 
}

double CG_Jacobi::solve(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) 
{


    // 1)  Compute Alpha 
    double alpha; 
    double qAq;
    
    // Aq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        adirection(i,j) = Discretization::laplacian(direction,i,j);    
    }
    // aTAq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        qAq += direction(i,j) * adirection(i,j); 
    }

    // Communicate to all processes before computing alpha 
    double total_denom = qAq; // TODO : Change to qAq if this explodes; 
    Communication::communicate_sum_double(&qAq,&total_denom); 
    alpha = sigma / (total_denom ) ;

    // 2) Update p and r 
        // p = p + alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        field.p(i,j) += alpha*direction(i,j); 
    }
    // Communicate pressure to all iterations 
    Communication::communicate_all(field.p_matrix(),MessageTag::P); 
    Communication::communicate_all(direction,MessageTag::P); 

    // Update residual 
    // r =r - alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= alpha * adirection(i,j); 
    }

    // r_bar = P^-1 * r 
    // Condition  the residual 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        cond_residual(i,j) = Discretization::jacobi(residual,i,j);
    }

    // Compute beta 
    double sigma_new = 0; 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        sigma_new += residual(i,j) * cond_residual(i,j); 
    }
    double return_val = sigma_new; 
    double total_sigma = 0;
    Communication::communicate_sum_double(&sigma_new,&total_sigma); 
    double beta = total_sigma / sigma; 
    sigma = sigma_new; 

    // Compute  new directions 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        direction(i,j) = cond_residual(i,j) + beta*(direction(i,j)); 
    }

#ifdef DEBUG 
    if(iter % 100)
    {
        std::cerr << "alpha = " << alpha << std::endl;
        std::cerr << "sigma = " << sigma << std::endl; 
        std::cerr << "beta = " << beta << std::endl;
        std::cerr << "sigma new " << sigma_new << std::endl;
        iter++;
    }
#endif 
    return fabs(sigma_new); 
}
