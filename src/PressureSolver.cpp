#include "PressureSolver.hpp"
#include "Communication.hpp"
// #define DEBUG 
#include <cmath>
#include <iostream>
#include <algorithm>
#include <array>

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

    Communication::communicate_all(residual, MessageTag::P);

    // apply boundary conditions on residual

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
    Communication::communicate_all(a_direction, MessageTag::P);

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
        // std::cout << "cell (" << i << ',' << j << ") residual(ij): " << residual(i,j) << std::endl;

    }
    // q1 = P^-1 r 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        cond_residual(i,j) = Discretization::jacobi(residual,i,j); 
        // std::cout << "cell (" << i << ',' << j << ") residual(ij), Disc::jacobi: " << (double)residual(i,j) << ',' << (double)Discretization::jacobi(residual,i,j) <<std::endl;

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
        #ifdef DEBUG
        // std::cout << "cell ( "<< i << ", " << j << "): qAq term = " << direction(i,j) * adirection(i,j) << std::endl;
        #endif
    }
    #ifdef DEBUG
        std::cout << "total qAq  = " << qAq << std::endl;
    #endif


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
    #ifdef DEBUG
        std::cout << "After Disc::jacobi.: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
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
    // if(iter % 100)
    {
        std::cerr << "alpha = " << alpha << std::endl;
        std::cerr << "sigma = " << sigma << std::endl; 
        std::cerr << "beta = " << beta << std::endl;
        std::cerr << "sigma new " << sigma_new << std::endl;
        std::cout << std::endl;
        iter++;
    }
#endif 
    return fabs(sigma_new); 
}

//-------------------------------------------------------------------------------

void CG_GS::init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries)
{
    #ifdef DEBUG
        std::cout << "\nI AM CG_GS::INIT !!!" << std::endl;
    #endif
    /* CG_GS initialization steps:
    1. r_0 = b - Ap
    2. r_cap_0 = M.inv * r_0
    3. sigma = r_0 dot (M.inv * r_0) = r_0 dot r_cap_0
    4. d_0 = r_cap_0
    */
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy(); 

    residual = Matrix<double>(imax+2, jmax+2); 
    cond_residual = Matrix<double>(imax+2, jmax+2, 0.0); 

    // As we know the pressure is non homogenous 
    
    // 1. r0 = b - A p 
    residual = field.rs_matrix();
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= Discretization::laplacian(field.p_matrix(),i,j);  
        // std::cout << "CG_GS::init(): cell (" << i << ',' << j << ") pressure(ij): " << field.p(i,j) << std::endl;
    }
    
    // 2. Solve P*rc = r 
    /*    (3-step process for GS preconditioner P = (L+D)*D.inv*(D+U))
            1. solving the system (D+L)*r' = r for r′using forward-substitution
            2. solving D.inv*r'' = r' for r'' by multiplying r' with D
            3. solving(D+U)*rc = r'' for  rc using backward-substitution.
        ref: https://www5.in.tum.de/pub/klimenko_idp.pdf [page 12]
    */
    double diag = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
    double diag_inv = 1.0 / diag;
    #ifdef DEBUG
        std::cout << " diag, diag_inv: \t" << diag << ", " << diag_inv << std::endl;   // = h^2 / 4.0, if dx == dy == h
    #endif

    // 2.1 Forward Substitution : (L+D)*r' = r
    for(int j = 0; j < jmax+2; j++) {
        for(int i = 0; i < imax+2; i++) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                diag_inv = 1 / Discretization::diagonal_term(i,j,grid);
                cond_residual(i,j) = diag_inv * ( residual(i,j) - Discretization::GS_Forward_Sub(cond_residual,i,j) );
                // std::cout << "cell (" << i << ',' << j << ") residual(ij), Disc: " << residual(i,j) << ',' << Discretization::GS_Forward_Sub(cond_residual,i,j) <<std::endl;
            }
        }    
    }
    #ifdef DEBUG
        std::cout << "After Forward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Forward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif
    
    // 2.2 Diagonal : D.inv * r'' = r' => r'' = D*r'
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        cond_residual(i,j) = Discretization::diagonal_term(i,j,grid) * cond_residual(i,j); 
    }

    // 2.3 Backward Substitution: (D+U)*r_cap = r''
    for (int i = imax+1; i >= 0; --i) {
        for (int j = jmax+1; j >= 0; --j) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                diag_inv = 1 / Discretization::diagonal_term(i,j,grid);
                cond_residual(i,j) = diag_inv * ( cond_residual(i,j) - Discretization::GS_Backward_Sub(cond_residual,i,j) );
                // if(i+1 == imax+1) std::cout << "\tcond_residual(imax+1,j) should be zero: "<< cond_residual(i+1,j) << std::endl;
                // if(j+1 == jmax+1) std::cout << "\tcond_residual(i,jmax+1) should be zero: "<< cond_residual(i,j+1) << std::endl;
                // if(i-1 == 0) std::cout << "\tcond_residual(0,j) should be zero: "<< cond_residual(i-1,j) << std::endl;
                // if(j-1 == 0) std::cout << "\tcond_residual(i,0) should be zero: "<< cond_residual(i,j-1) << std::endl;
            }
        }
    }
    #ifdef DEBUG
        std::cout << "After Backward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Backward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // sigma_0 = <r0,q1> 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        sigma += residual(i,j) * cond_residual(i,j);  
    }
    direction = cond_residual; 
    adirection = direction; 
    double total_sigma = 0; 
    Communication::communicate_sum_double(&sigma, &total_sigma); 
    sigma = total_sigma; 
}

double CG_GS::solve(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) 
{
    #ifdef DEBUG
        std::cout << "\nI AM CG_GS::SOLVE !!!" << std::endl;
    #endif
    /* Note: Terminology --> direction <=> M.inv*r <=> "q" <=> "cond_residual"
    1. Compute alpha = r.q / q.Aq = sigma / q.Aq
    2. Update p(ij) = p_old(ij) + alpha*q(ij)
    3. Update r(ij) = r_old(ij) - alpha*A.q(ij)
    4. Condition the residual ( Solve P*q = r )
    5. Compute new directions ( d(ij) = r_cap(ij) + beta * d_old(ij) )
    */
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy();
    double residual_magnitude_square = 0;
    
    // 1)  Compute Alpha 
    double alpha; 
    double qAq{0};
    
    // 1.1) Aq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        adirection(i,j) = Discretization::boundary_aware_laplacian(direction,i,j, grid);    
    }
    // 1.2) q_T * Aq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        qAq += direction(i,j) * adirection(i,j); 
        #ifdef DEBUG
            // std::cout << "cell ( "<< i << ", " << j << "): qAq term = " << direction(i,j) * adirection(i,j) << std::endl;
        #endif
    }
    #ifdef DEBUG
        std::cout << "total qAq  = " << qAq << std::endl;
    #endif

    // 1.3) Communicate to all processes before computing alpha 
    double total_denom = qAq; // TODO : Change to qAq if this explodes (?); 
    Communication::communicate_sum_double(&qAq,&total_denom); 

    alpha = sigma / total_denom;
    #ifdef DEBUG
        std::cerr << "alpha=sigma/denom = " << sigma << " / " << total_denom << std::endl;
    #endif

    // 2) Update p //   p = p + alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        field.p(i,j) += alpha*direction(i,j); 
    }
    // Communicate pressure to all iterations 
    Communication::communicate_all(field.p_matrix(),MessageTag::P); 
    Communication::communicate_all(direction,MessageTag::P); 

    // 3) Update residual //    r = r - alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= alpha * adirection(i,j); 
        residual_magnitude_square += residual(i,j) * residual(i,j);
    }

    // 4) Compute beta
    //    4.1) Condition  the residual: //         (Solve P*q1 = r )
    //         (3-step process for GS preconditioner P = (L+D)*D.inv*(D+U)) described in CG_GS::init()
    double diag = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));  /// TODO: remove this
    double diag_inv = 1.0 / diag;   /// TODO: remove this
    #ifdef DEBUG
        std::cout << " dx,dy,diag, diag_inv: \t" << dx << ", " << dy << ", " << diag << ", " << diag_inv << std::endl;   // = h^2 / 4.0, if dx == dy == h
    #endif
    
    // 4.1.1 Forward Substitution : (L+D)*r' = r
    for(int j = 0; j < jmax+2; j++) {
        for(int i = 0; i < imax+2; i++) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
            diag_inv = 1 / Discretization::diagonal_term(i,j,grid);
            cond_residual(i,j) = diag_inv * ( residual(i,j) - Discretization::GS_Forward_Sub(cond_residual,i,j) );
            }
        }    
    }
    #ifdef DEBUG
        std::cout << "After Forward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Forward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // 4.1.2 Diagonal : D.inv * r'' = r' => r'' = D*r'
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        cond_residual(i,j) = Discretization::diagonal_term(i,j,grid) * cond_residual(i,j); 
    }
    #ifdef DEBUG
        std::cout << "After Diagonal mult.: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    // 4.1.3 Backward Substitution: (D+U)r_cap = r''
    for (int i = imax+1; i >= 0; --i) {
        for (int j = jmax+1; j >= 0; --j) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                diag_inv = 1 / Discretization::diagonal_term(i,j,grid);
                cond_residual(i,j) = diag_inv * ( cond_residual(i,j) - Discretization::GS_Backward_Sub(cond_residual,i,j) );
                // if(i = imax-1) std::cout << "\tcond_residual(imax,j) should be zero: "<< cond_residual(i+1,j) << std::endl;
                // if(i = jmax-1) std::cout << "\tcond_residual(i,jmax) should be zero: "<< cond_residual(i,j+1) << std::endl;
            }
        }
    }
    #ifdef DEBUG
        std::cout << "After Backward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Backward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // 4.2) Compute denominator = r dot r_cap 
    double sigma_new = 0; 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        sigma_new += residual(i,j) * cond_residual(i,j); 
    }
    double return_val = sigma_new; 
    double total_sigma_new = 0;
    Communication::communicate_sum_double(&sigma_new,&total_sigma_new); 

    // 4.3) Compute beta
    double beta = total_sigma_new / sigma; 
    sigma = sigma_new; //  should be total_sigma_new?

    // 5) Compute new directions 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        direction(i,j) = cond_residual(i,j) + beta*(direction(i,j)); 
    }

#ifdef DEBUG 
    // if(iter % 100)
    {
        // std::cerr << "iteration = " << iter << std::endl;
        std::cerr << "alpha = " << alpha << std::endl;
        // std::cerr << "sigma = " << sigma << std::endl; 
        std::cerr << "beta = " << beta << std::endl;
        // std::cerr << "sigma new " << sigma_new << std::endl;
        std::cerr << std::endl;
        iter++;
    }
#endif 
    return residual_magnitude_square; 
}
//-------------------------------------------------------------------------------

void CG_DILU::init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries)
{
    #ifdef DEBUG
        std::cout << "\nI AM CG_DILU::INIT !!!" << std::endl;
    #endif
    /* CG_GS initialization steps:
    1. r_0 = b - Ap
    2. r_cap_0 = M.inv * r_0
    3. sigma = r_0 dot (M.inv * r_0) = r_0 dot r_cap_0
    4. d_0 = r_cap_0
    */
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy(); 

    residual = Matrix<double>(imax+2, jmax+2); 
    cond_residual = Matrix<double>(imax+2, jmax+2, 0.0); 

    // As we know the pressure is non homogenous 
    
    // 1. r0 = b - A p 
    residual = field.rs_matrix();
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= Discretization::laplacian(field.p_matrix(),i,j);  
        // std::cout << "cell (" << i << ',' << j << ") residual(ij): " << residual(i,j) << std::endl;
    }
    
    // 2. Solve P*rc = r 
    /*    (3-step process for GS preconditioner P = (L+D)*D.inv*(D+U))
            1. solving the system (D+L)*r' = r for r′using forward-substitution
            2. solving D.inv*r'' = r' for r'' by multiplying r' with D
            3. solving(D+U)*rc = r'' for  rc using backward-substitution.
        ref: https://www5.in.tum.de/pub/klimenko_idp.pdf [page 12]
    */

    // SETUP DIAGONALS VECTOR
    // ref: http://www.netlib.org/templates/templates.pdf [page 51]
    double diag = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
    double diag_inv = 1.0 / diag;    
    double x_coeff = 1 / (dx*dx);
    double y_coeff = 1 / (dy*dy);
    int k{0}; // iterator
    #ifdef DEBUG
        std::cout << " diag, diag_inv: \t" << diag << ", " << diag_inv << std::endl;   // = h^2 / 4.0, if dx == dy == h
    #endif

    // std::array<double, imax*jmax> D_array;
    D_vector.reserve(imax*jmax);
        std::cout << "Size of D_vector: " << D_vector.size();

    // std::fill(D_vector.begin(), D_vector.end(), diag);
    for(k = 0; k < (imax*jmax); k++) {
        D_vector.push_back(diag);
    }
    for(k = 0; k < (imax*jmax); k++) {
        if ((k+1)%imax) { // this needs more conditions . Sometimes x_coeff is zero
            D_vector[k+1] -= x_coeff * x_coeff / D_vector[k,k];
        }
        if ((k+1+imax) <= imax*jmax) { // needs more conditions, sometimes y_coeff is zero
            D_vector[k+imax] -= y_coeff * y_coeff / D_vector[k,k];
        }
    }
    /// TODO: verify if D_vector is correct
    /// TODO: Store inverse of D_vector. We mostly use the inverse.

    // 2.1 Forward Substitution : (L+D)*r' = r
    for(int j = 0; j < jmax; j++) {
        for(int i = 0; i < imax; i++) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                k = imax*j + i;
                cond_residual(i,j) = 1/D_vector[k] * ( residual(i,j) - Discretization::GS_Forward_Sub(cond_residual,i,j) );
                // std::cout << "cell (" << i << ',' << j << ") residual(ij), Disc: " << residual(i,j) << ',' << Discretization::GS_Forward_Sub(cond_residual,i,j) <<std::endl;

            }
        }    
    }
    #ifdef DEBUG
        std::cout << "After Forward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Forward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif
    // 2.2 Diagonal : D.inv * r'' = r' => r'' = D*r'
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        k = imax*j + i;
        cond_residual(i,j) = D_vector[k] * cond_residual(i,j); 
    }
    // 2.3 Backward Substitution: (D+U)r_cap = r''
    for (int i = imax-1; i >= 0; --i) {
        for (int j = jmax-1; j >= 0; --j) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                k = imax*j + i;
                cond_residual(i,j) = 1/D_vector[k] * ( cond_residual(i,j) - Discretization::GS_Backward_Sub(cond_residual,i,j) );
                #ifdef DEBUG
                // if(i = imax-1) std::cout << "\tcond_residual(imax,j) should be zero: "<< cond_residual(i+1,j) << std::endl;
                // if(i = jmax-1) std::cout << "\tcond_residual(i,jmax) should be zero: "<< cond_residual(i,j+1) << std::endl;
                #endif
            }
        }
    }
    #ifdef DEBUG
        std::cout << "After Backward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Backward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // sigma_0 = <r0,q1> 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        sigma += residual(i,j) * cond_residual(i,j);  
    }
    direction = cond_residual; 
    adirection = direction; 
    double total_sigma = 0; 
    Communication::communicate_sum_double(&sigma, &total_sigma); 
    sigma = total_sigma; 
}

double CG_DILU::solve(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) 
{
    #ifdef DEBUG
        std::cout << "\nI AM CG_DILU::SOLVE !!!" << std::endl;
    #endif
    /* Note: Terminology --> direction <=> M.inv*r <=> "q" <=> "cond_residual"
    1. Compute alpha = r.q / q.Aq = sigma / q.Aq
    2. Update p(ij) = p_old(ij) + alpha*q(ij)
    3. Update r(ij) = r_old(ij) - alpha*A.q(ij)
    4. Condition the residual ( Solve P*q = r )
    5. Compute new directions ( d(ij) = r_cap(ij) + beta * d_old(ij) )
    */
    int imax = grid.imax();
    int jmax = grid.jmax();
    double dx = grid.dx();
    double dy = grid.dy(); 
    
    // 1)  Compute Alpha 
    double alpha; 
    double qAq;
    
    // 1.1) Aq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        adirection(i,j) = Discretization::laplacian(direction,i,j);    
    }
    // 1.2) q_T * Aq 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j();
        qAq += direction(i,j) * adirection(i,j); 
        #ifdef DEBUG
            // std::cout << "cell ( "<< i << ", " << j << "): qAq term = " << direction(i,j) * adirection(i,j) << std::endl;
        #endif
    }
    #ifdef DEBUG
        std::cout << "total qAq  = " << qAq << std::endl;
    #endif

    // 1.3) Communicate to all processes before computing alpha 
    double total_denom = qAq; // TODO : Change to qAq if this explodes (?); 
    Communication::communicate_sum_double(&qAq,&total_denom); 
    alpha = sigma / total_denom;
    #ifdef DEBUG
        std::cerr << "alpha=sigma/denom = " << sigma << " / " << total_denom << std::endl;
    #endif

    // 2) Update p //   p = p + alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        field.p(i,j) += alpha*direction(i,j); 
    }
    // Communicate pressure to all iterations 
    Communication::communicate_all(field.p_matrix(),MessageTag::P); 
    Communication::communicate_all(direction,MessageTag::P); 

    // 3) Update residual //    r = r - alpha*q 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        residual(i,j) -= alpha * adirection(i,j); 
    }

    /* 4) Compute beta
          4.1) Condition  the residual: //         (Solve P*q1 = r )
              ( 3-step process for GS preconditioner P = (L+D)*D.inv*(D+U)) described in CG_GS::init()
    */
    // double diag = -2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy));
    // double diag_inv = 1.0 / diag;   // = h^2 / 4.0, if dx == dy == h
    #ifdef DEBUG
        // std::cout << " dx,dy,diag, diag_inv: \t" << dx << ", " << dy << ", " << diag << ", " << diag_inv << std::endl;   // = h^2 / 4.0, if dx == dy == h
    #endif
    int k{0}; // iterator
    // 4.1.1 Forward Substitution : (L+D)*r' = r
    for(int j = 0; j < jmax; j++) {
        for(int i = 0; i < imax; i++) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
            // if (i != 1) assert(cond_residual(i-1,j) != 0);
            // if (j != 1) assert(cond_residual(i,j-1) != 0);
            k = imax*j + i; 
            cond_residual(i,j) = 1/D_vector[k] * ( residual(i,j) - Discretization::GS_Forward_Sub(cond_residual,i,j) );
            }
        }    
    }
    #ifdef DEBUG
        std::cout << "After Forward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Forward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // 4.1.2 Diagonal : D.inv * r'' = r' => r'' = D*r'
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        k = imax*j + i;
        cond_residual(i,j) = D_vector[k] * cond_residual(i,j); 
    }
    #ifdef DEBUG
        std::cout << "After Diagonal mult.: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    // 4.1.3 Backward Substitution: (D+U)r_cap = r''
    for (int i = imax-1; i >= 0; --i) {
        for (int j = jmax-1; j >= 0; --j) {
            if(grid.cell(i,j).type() == cell_type::FLUID ) {
                k = imax*j + i;
                cond_residual(i,j) = 1/D_vector[k] * ( cond_residual(i,j) - Discretization::GS_Backward_Sub(cond_residual,i,j) );
                // if(i = imax-1) std::cout << "\tcond_residual(imax,j) should be zero: "<< cond_residual(i+1,j) << std::endl;
                // if(i = jmax-1) std::cout << "\tcond_residual(i,jmax) should be zero: "<< cond_residual(i,j+1) << std::endl;
            }
        }
    }
    #ifdef DEBUG
        std::cout << "After Backward_sub: cond_residual[25,25]: " << cond_residual(25,25) << std::endl;
    #endif
    #ifdef DEBUG_2
        std::cout << "cond_residual after Backward_sub: " << std::endl;
        std::cout.precision(2);
        cond_residual.pretty_print(std::cout);
    #endif

    // 4.2) Compute denominator = r dot r_cap 
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

    // 4.3) Compute beta
    double beta = total_sigma / sigma; 
    sigma = sigma_new; 

    // 5) Compute new directions 
    for(auto cell:grid.fluid_cells())
    {
        int i = cell->i();
        int j = cell->j(); 
        direction(i,j) = cond_residual(i,j) + beta*(direction(i,j)); 
    }

#ifdef DEBUG 
    // if(iter % 100)
    {
        // std::cerr << "iteration = " << iter << std::endl;
        std::cerr << "alpha = " << alpha << std::endl;
        // std::cerr << "sigma = " << sigma << std::endl; 
        std::cerr << "beta = " << beta << std::endl;
        // std::cerr << "sigma new " << sigma_new << std::endl;
        std::cerr << std::endl;
        iter++;
    }
#endif 
    return fabs(sigma_new); 
}
