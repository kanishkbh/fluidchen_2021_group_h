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

double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double u, u1, u2, u3, u4, u5, v, v1, v2, v3, v4, v6;
 double u2x, uv_y, uxx, uyy, v2y, uv_x, vxx, vyy;
	
	u = U(i,j);
    u1 = U(i + 1,j);
    u2 = U(i,j + 1);
    u3 = U(i - 1,j);
    u4 = U(i,j - 1);
    u5 = U(i - 1,j + 1);

    v = V(i,j);
    v1 = V(i + 1,j);
    v2 = V(i,j + 1);
    v3 = V(i - 1,j);
    v4 = V(i,j - 1);
    v6 = V(i + 1,j - 1);

	  u2x =
	    1 / _dx * ( ((u + u1)*(u + u1)/4) - ((u3 + u)*(u3 + u)/4) ) +
	    _gamma / (4 * _dx) * (fabs (u + u1) * (u - u1) -
				fabs (u3 + u) * (u3 - u));
	  uv_y =
	    1 / (4 * _dy) * ((v + v1) * (u + u2) - (v4 + v6) * (u4 + u)) +
	    _gamma / (4 * _dy) * (fabs (v + v1) * (u - u2) -
				fabs (v4 + v6) * (u4 - u));
    return u2x + uv_y;
}

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double u, u1, u2, u3, u4, u5, v, v1, v2, v3, v4, v6;
 double u2x, uv_y, uxx, uyy, v2y, uv_x, vxx, vyy;
	
	u = U(i,j);
    u1 = U(i + 1,j);
    u2 = U(i,j + 1);
    u3 = U(i - 1,j);
    u4 = U(i,j - 1);
    u5 = U(i - 1,j + 1);

    v = V(i,j);
    v1 = V(i + 1,j);
    v2 = V(i,j + 1);
    v3 = V(i - 1,j);
    v4 = V(i,j - 1);
    v6 = V(i + 1,j - 1);

    v2y =
	    1 / _dy * ( (0.5 * (v + v2))*(0.5 * (v + v2)) - (0.5 * (v4 + v))*(0.5 * (v4 + v)) ) +
	    _gamma / (4 * _dy) * (fabs (v + v2) * (v - v2) -
				fabs (v4 + v) * (v4 - v));
	uv_x =
	    1 / (4 * _dx) * ((u + u2) * (v + v1) - (u3 + u5) * (v3 + v)) +
	    _gamma / (4 * _dx) * (fabs (u + u2) * (v - v1) -
				fabs (u3 + u5) * (v3 - v));
    return v2y + uv_x;
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