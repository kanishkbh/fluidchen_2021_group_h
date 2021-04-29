#include "Discretization.hpp"

#include <cmath>

double Discretization::_dx = 0.0;
double Discretization::_dy = 0.0;
double Discretization::_gamma = 0.0;

Discretization::Discretization(double dx, double dy, double gamma) {
    _dx = dx;
    _dy = dy;
    _gamma = gamma;
}

double Discretization::convection_u(const Matrix<double> &_U, const Matrix<double> &_V, int i, int j) {

    // loop invariants 
    double idx = 1/_dx;
    double idy = 1/_dy;
    // k3 := discretization of Dx(u**2)
    double k3 = 0.25*idx*( ( (_U(i,j) + _U(i+1,j))*(_U(i,j) + _U(i+1,j))  - (_U(i-1,j) + _U(i,j))*(_U(i-1,j) + _U(i,j)) )
        + _gamma*( abs(_U(i,j) + _U(i+1,j)) * (_U(i,j) - _U(i+1,j)) - abs(_U(i-1,j) + _U(i,j)) * (_U(i-1,j) - _U(i,j)) ) );
    // k4:= discretization of Dy(u*v)
    double k4 = 0.25*idy*( ( (_V(i,j) + _V(i+1,j)) * (_U(i,j) + _U(i,j+1)) ) -  (_V(i,j-1)+_V(i+1,j-1)) * (_U(i,j-1) + _U(i,j)) )
        + _gamma*( (abs(_V(i,j) + _V(i+1,j))*(_U(i,j)-_U(i,j+1))) - (abs(_V(i,j-1) + _V(i+1,j-1)) * (_U(i,j-1)-_U(i,j)) ) );

    return k3+k4;
}

double Discretization::convection_v(const Matrix<double> &_U, const Matrix<double> &_V, int i, int j) {

    // loop invariants xxxxxxxxxxxx Not really, this is a refactored code xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    double idx = 1/_dx;
    double idy = 1/_dy;
    // k5:= discretization of Dx(uv)
    double k5 = 0.25*idx*( ( (_U(i,j) + _U(i,j+1))*(_V(i,j) + _V(i+1,j))  - (_U(i-1,j) + _U(i-1,j+1))*(_V(i-1,j) + _V(i,j)) )  
    + _gamma*( abs(_U(i,j) + _U(i,j+1)) * (_V(i,j) - _V(i+1,j)) - abs(_U(i-1,j) + _U(i-1,j+1)) * (_V(i-1,j) - _V(i,j)) ) );
    // k6:= disretization of Dy(v**2)
    double k6 =0.25*idy*( ( (_V(i,j)+_V(i,j+1))*(_V(i,j)+_V(i,j+1)) - (_V(i,j-1)+_V(i,j))*(_V(i,j-1)+_V(i,j)) ) 
    + _gamma*(abs(_V(i,j)+_V(i,j+1))*(_V(i,j)-_V(i,j+1)) - abs(_V(i,j-1)+_V(i,j))*(_V(i,j-1)-_V(i,j)) ) );

    return k5+k6;
}

double Discretization::diffusion(const Matrix<double> &_U, int i, int j) {
        
        // loop invariants xxxxxxxxxxxx Not really, this is a refactored code xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        double idx2 = 1/(_dx*_dx);
        double idy2 = 1/(_dy*_dy);
        // k1 := discretization of D2x u 
        double k1 = idx2*(_U(i+i,j) + _U(i-1,j) - 2*_U(i,j));
        // k2 := discretization of D2y u
        double k2 = idy2*(_U(i,j+1) + _U(i,j-1) - 2*_U(i,j));

        return k1+k2;
}

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}

double Discretization::interpolate(const Matrix<double> &A, int i, int j, int i_offset, int j_offset) {}