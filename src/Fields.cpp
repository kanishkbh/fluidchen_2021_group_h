#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include<cmath>

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI, double TI,
               double alpha, double beta, double gx, double gy)
    : _nu(nu), _dt(dt), _tau(tau), _alpha(alpha), _beta(beta), _U(imax + 2, jmax + 2, UI), _V(imax + 2, jmax + 2, VI),
      _P(imax + 2, jmax + 2, PI), _T(imax + 2, jmax + 2, TI), _F(imax + 2, jmax + 2, 0.0), _G(imax + 2, jmax + 2, 0.0),
      _RS(imax + 2, jmax + 2, 0.0), _gx(gx), _gy(gy) {


}

void Fields::calculate_fluxes(Grid &grid) {


    for (auto cell : grid.fluid_cells()) {
        auto i = cell->i();
        auto j = cell->j();
        //If it's a boundary, leave it to the BCs to fix its value
        if (cell->neighbour(border_position::RIGHT)->type() == cell_type::FLUID)
        {
            _F(i,j) = _U(i,j) + _dt*(_nu*(Discretization::laplacian(_U,i,j))
                                 - Discretization::convection_u(_U,_V,i,j) + _gx 
                                 - _gx * _beta * 0.5*(_T(i,j) + _T(i+1, j)));
        }
        if (cell->neighbour(border_position::TOP)->type() == cell_type::FLUID)
        {
            _G(i,j) = _V(i,j) + _dt*(_nu*(Discretization::laplacian(_V,i,j))
                                 -Discretization::convection_v(_U,_V,i,j) + _gy
                                 - _gy * _beta * 0.5 *(_T(i,j) + _T(i, j+1)));
        }
        
    }

}

void Fields::calculate_T(Grid &grid) {

    // In-place update !
    auto new_T = _T;
    

    for (auto cell : grid.fluid_cells()) {
        auto i = cell->i();
        auto j = cell->j();
        new_T(i,j) = _T(i,j) + _dt*(_alpha*(Discretization::laplacian(_T,i,j))
                                 - Discretization::convection_u_T(_U,_T,i,j) - Discretization::convection_v_T(_V,_T,i,j));
    }

    //Replace _T par updated version
    _T = std::move(new_T);

}

//-----------------------------------------------------------------------------------------------------------
// Compute gradients of force fields 
void Fields::calculate_rs(Grid &grid) {
    double idt = 1/_dt;
    double idx = 1/grid.dx();
    double idy = 1/grid.dy();

    for (auto cell : grid.fluid_cells()) {
        auto i = cell->i();
        auto j = cell->j();
        _RS(i,j) = idt*( idx*(_F(i,j)-_F(i-1,j)) + idy*(_G(i,j)-_G(i,j-1)) );
    }
}
//-----------------------------------------------------------------------------------------------------------
// Compute velocity updates 
void Fields::calculate_velocities(Grid &grid) {
    double kappa1 = _dt/grid.dx();
    double kappa2 = _dt/grid.dy();
    int jmax = grid.jmax();
    int imax = grid.imax();

    for (auto cell_ptr : grid.fluid_cells()) {
        int i = cell_ptr->i();
        int j = cell_ptr->j();
        _U(i,j) = _F(i,j) - kappa1*(_P(i+1,j)-_P(i,j));
        _V(i,j) = _G(i,j) - kappa2*(_P(i,j+1)-_P(i,j));
    }
    
}
//-----------------------------------------------------------------------------------------------------------

double Fields::calculate_dt(Grid &grid) { 

    // If the safety parameter tau is negative, this is interpreted as "use a fixed dt".

    if(_tau < 0){
        return _dt;
    } 
    else {
        double dx = grid.dx();
        double dy = grid.dy();
        double k1 = (0.5/_nu) * 1/(1/(dx*dx)+(1/dy*dy));
        double k2 = dx/(_U.max() + 1e-8); //Epsilon to ensure no division by 0
        double k3 = dy/(_V.max() + 1e-8);
        double k4 = (0.5/_alpha) * 1/((1/(dx*dx))+(1/(dy*dy)));
        _dt = _tau * std::min({k1, k2, k3, k4});
        
        return _dt;
    }
}
//-----------------------------------------------------------------------------------------------------------
double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::t(int i, int j) { return _T(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
