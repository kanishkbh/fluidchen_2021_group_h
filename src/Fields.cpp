#include "Fields.hpp"

#include <algorithm>
#include <iostream>
#include<cmath> 

Fields::Fields(double nu, double dt, double tau, int imax, int jmax, double UI, double VI, double PI)
    : _nu(nu), _dt(dt), _tau(tau) {
    _U = Matrix<double>(imax + 2, jmax + 2, UI);
    _V = Matrix<double>(imax + 2, jmax + 2, VI);
    _P = Matrix<double>(imax + 2, jmax + 2, PI);

    _F = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _G = Matrix<double>(imax + 2, jmax + 2, 0.0);
    _RS = Matrix<double>(imax + 2, jmax + 2, 0.0);
}

void Fields::calculate_fluxes(Grid &grid) {

    // Template fill away :  0.25*idy*( ( ()*() - ()*() ) + gamma*(abs()*() - abs()*() ) ); 
    // donor cell scheme. 
    double gamma = 1;   
    // Create Discretization object; 
    Discretization del(grid.dx(),grid.dy(),gamma);
    for(auto j=1;j<=jmax;j++){
        for(auto i=1;i<imax;++i){
        _F(i,j) = _U(i,j) + _dt*(_nu*(del.diffusion(_U,i,j))
                                 - del.convection_u(_U,_V,i,j) +_gx);
        }
    }
        for(auto j=1;j<=jmax;++j){
            _F(0,j) = _U(0,j);
            _F(imax,j) = _U(imax,j);
        }
    
        
    for(auto j=1;j<jmax;++j){
        for(auto i=1;i<=imax;++i){        
        _G(i,j) = _V(i,j) + _dt*(_nu*(del.diffusion(_V,i,j))
                                 -del.convection_v(_U,_V,i,j) +_gy);
        }
    }
    for(auto i=1;i<=imax;++i){
        _G(i,0) = _V(i,0);
        _G(i,jmax) = _V(i,jmax);
    }
        
}

void Fields::calculate_rs(Grid &grid) {
    double idt = 1/_dt;
    double idx = 1/grid._dx;
    double idy = 1/grid._dy;
    for(auto j=1;j<=grid.jmax();++j){
        for(auto i=1;i<=grid.imax();++i){
            _RS(i,j) = idt*( idx*(_F(i,j)-_F(i-1,j)) + idy*(_G(i,j)-_G(i,j-1)) );
        }
    } 
}

void Fields::calculate_velocities(Grid &grid) {
    double kappa = _dt/grid.dx();
    int jmax = grid.jmax();
    int imax = grid.imax();
    for(auto j=1;j<=jmax;++j){
        for (auto i=1;i<imax;++i){
            _U(i,j) = _F(i,j) - kappa*(_P(i+1,j)-_P(i,j));
        }
    }
    for(auto j=1;j<jmax;++j){
        for(auto i=1;i<=imax;++i){
            _V(i,j) = _G(i,j) - kappa*(_P(i,j+1)-P(i,j));
        }
    }
}

double Fields::calculate_dt(Grid &grid) { 

    // If the safety parameter tau is negative, it makes little sense. 
    // Hence for a negative tau, we simply have a fixed time step. 
    if(_tau<0){
        return _dt;
    } 
    else {
        double dx = grid.dx();
        double dy = grid.dy();
        double k1 = (0.5/_nu) * 1/(1/(dx*dx)+(1/dy*dy));
        double k2 = dx/(_U.max());
        double k3 = dy/(_V.max());
        _dt = _tau*std::min(k1,k2,k3);
        return _dt;
    }
}

double &Fields::p(int i, int j) { return _P(i, j); }
double &Fields::u(int i, int j) { return _U(i, j); }
double &Fields::v(int i, int j) { return _V(i, j); }
double &Fields::f(int i, int j) { return _F(i, j); }
double &Fields::g(int i, int j) { return _G(i, j); }
double &Fields::rs(int i, int j) { return _RS(i, j); }

Matrix<double> &Fields::p_matrix() { return _P; }

double Fields::dt() const { return _dt; }
