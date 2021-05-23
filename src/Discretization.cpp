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
//-----------------------------------------------------------------------------------------------------------
double Discretization::convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
	
	double idx = 1/_dx;
    double idy = 1/_dy;
    // temporary variables 

    double u = U(i,j);
    double u_right = U(i+1,j);
    double u_left = U(i-1,j);
    double u_top = U(i,j+1);
    double u_bottom = U(i,j-1);
    double u_west = U(i-1,j+1);
    double v = V(i,j);
    double v_left = V(i-1,j);
    double v_right = V(i+1,j);
    double v_top = V(i,j+1);
    double v_bottom = V(i,j-1);
    double v_east = V(i+1,j-1);

	double k3 = idx*0.25 *(
                        (u+u_right)*(u+u_right) - (u_left+u)*(u_left+u) +  _gamma *( fabs(u+u_right)*(u-u_right)  - fabs(u_left + u) * (u_left - u) )  
                          );

    double k4 = idy*0.25 *(
                        (v+v_right)*(u+u_top) - (v_bottom+v_east)*(u_bottom+u) + _gamma*(fabs(v+v_right)*(u-u_top) - fabs(v_bottom+v_east)*(u_bottom - u) )
                          );

    return k3+k4;
}
//-----------------------------------------------------------------------------------------------------------

double Discretization::convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j) {
    double idx = 1/_dx;
    double idy = 1/_dy;

    double u = U(i,j);
    double u_right = U(i+1,j);
    double u_left = U(i-1,j);
    double u_top = U(i,j+1);
    double u_bottom = U(i,j-1);
    double u_west = U(i-1,j+1);
    
    double v = V(i,j);
    double v_left = V(i-1,j);
    double v_right = V(i+1,j);
    double v_top = V(i,j+1);
    double v_bottom = V(i,j-1);
    double v_east = V(i+1,j-1);

    double k5 = idx*0.25 *(
                        (u+u_top)*(v+v_right) - (u_left+u_west)*(v_left+v) + _gamma*( fabs(u+u_top)*(v-v_right) - fabs(u_left+u_west)*(v_left - v))
                          );
    double k6 = idy*0.25*(
                        (v+v_top)*(v+v_top) - (v_bottom+v)*(v_bottom+v) + _gamma * (fabs(v+v_top)*(v-v_top) - fabs(v_bottom+v)*(v_bottom - v) ) 
                         );

    return k5+k6;
}
//-----------------------------------------------------------------------------------------------------------
double Discretization::convection_u_T(const Matrix<double> &U, const Matrix<double> &T, int i, int j) {
	
	double idx = 1/_dx;
    // temporary variables 

    double u = U(i,j);
    double u_left = U(i-1,j);
    
    double t = T(i,j);
    double t_left = T(i-1,j);
    double t_right = T(i+1,j);
    
    /* See lecture notes 6.4 for derivation */
	double d_uT_1 = 0.5 * idx * (u * (t + t_right) - u_left * (t_left + t));
    double d_uT_2 = 0.5 * idx * _gamma * (fabs(u) * (t - t_right) - fabs(u_left) * (t_left - t));



    return d_uT_1 + d_uT_2;
}

//-----------------------------------------------------------------------------------------------------------
double Discretization::convection_v_T(const Matrix<double> &V, const Matrix<double> &T, int i, int j) {
	
    double idy = 1/_dy;

    double v = V(i,j);
    double v_bottom = V(i,j-1);
    
    double t = T(i,j);
    double t_top = T(i,j+1);
    double t_bottom = T(i,j-1);
    

    /* See lecture notes 6.4 for derivation */
	double d_vT_1 = 0.5 * idy * (v * (t + t_top) - v_bottom * (t_bottom + t));
    double d_vT_2 = 0.5 * idy * _gamma * (fabs(v) * (t - t_top) - fabs(v_bottom) * (t_bottom - t));



    return d_vT_1 + d_vT_2;
}
//-----------------------------------------------------------------------------------------------------------

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
//-----------------------------------------------------------------------------------------------------------

double Discretization::laplacian(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) - 2.0 * P(i, j) + P(i - 1, j)) / (_dx * _dx) +
                    (P(i, j + 1) - 2.0 * P(i, j) + P(i, j - 1)) / (_dy * _dy);
    return result;
}
//-----------------------------------------------------------------------------------------------------------
double Discretization::sor_helper(const Matrix<double> &P, int i, int j) {
    double result = (P(i + 1, j) + P(i - 1, j)) / (_dx * _dx) + (P(i, j + 1) + P(i, j - 1)) / (_dy * _dy);
    return result;
}
