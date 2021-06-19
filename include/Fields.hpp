#pragma once

#include "Datastructures.hpp"
#include "Discretization.hpp"
#include "Grid.hpp"

/**
 * @brief Class of container and modifier for the physical fields
 *
 */
class Fields {
  public:
    Fields() = default;

    /**
     * @brief Constructor for the fields
     *
     * @param[in] kinematic viscosity
     * @param[in] initial timestep size (or fixed)
     * @param[in] adaptive timestep coefficient
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] initial x-velocity
     * @param[in] initial y-velocity
     * @param[in] initial pressure
     * @param[in] initial temperature
     * @param[in] Thermal diffusivity
     * @param[in] Thermal coefficient expansion (default ~= air at 20°C)
     * @param[in] Gravity in x-direction
     * @param[in] Gravity in y-direction
     */
    Fields(double _nu, double _dt, double _tau, int imax, int jmax, double UI, double VI, double PI, double TI = 0, double alpha = 0.0, double beta=0.0034, double gx=0, double gy=0);

    /**
     * @brief Calculates the convective and diffusive pseud-fluxes in x and y
     * direction based on explicit discretization of the momentum equations
     *
     * @param[in] grid in which the fluxes are calculated
     *
     */
    void calculate_fluxes(Grid &grid);

    /**
     * @brief Update temperatures through a convection-diffusion time step
     *
     * @param[in] grid in which the temperatures are calculated
     *
     */
    void calculate_T(Grid &grid);

    /**
     * @brief Right hand side calculations using the fluxes for the pressure
     * Poisson equation
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_rs(Grid &grid);

    /**
     * @brief Velocity calculation using pressure values and F&G fields
     *
     * @param[in] grid in which the calculations are done
     *
     */
    void calculate_velocities(Grid &grid);

    /**
     * @brief Adaptive step size calculation using x-velocity condition,
     * y-velocity condition and CFL condition. If fixed dt is used, return it.
     *
     * @param[in] grid in which the calculations are done
     *
     */
    double calculate_dt(Grid &grid);

    /**
     * @brief x-velocity index based access and modify
     * */
    double &u(int i, int j);

    /**
     * @brief y-velocity index based access and modify
     * */
    double &v(int i, int j);

    /**
     * @brief pressure index based access and modify
     * */
    double &p(int i, int j);

    /**
     * @brief PPE right hand side index based access and modify
     * */
    double &rs(int i, int j);

    /**
     * @brief x-momentum pseudo-flux (F field) index based access and modify
     * */
    double &f(int i, int j);

    /**
     * @brief y-momentum pseudo-flux (G field) index based access and modify
     * */
    double &g(int i, int j);

    /**
     * @brief Temperature index based access and modify
     * */
    double &t(int i, int j);

    /**
     * @brief get current timestep (without recomputing it)
     * */
    double dt() const;

    /**
     * @brief Access a matrix reference to the pressure
     * */
    Matrix<double> &p_matrix();

    /**
     * @brief Access a matrix reference to the temperature
     * */
    Matrix<double> &t_matrix() {return _T;}

     /**
     * @brief Access a matrix reference to F
     * */
    Matrix<double> &f_matrix() {return _F;}

     /**
     * @brief Access a matrix const reference to G
     * */
    Matrix<double> &g_matrix() {return _G;}

    /**
     * @brief Access a matrix const reference to U
     * */
    const Matrix<double> &u_matrix() const {return _U;}

    /**
     * @brief Access a matrix const reference to V
     * */
    const Matrix<double> &v_matrix() const {return _V;}

    /**
     * @brief Access a matrix non-const reference to U
     * */
    Matrix<double> &u_matrix() {return _U;}

    /**
     * @brief Access a matrix non-const reference to V
     * */
    Matrix<double> &v_matrix() {return _V;}

    /**
     * @brief Access a matrix const reference to T
     * */
    const Matrix<double> &t_matrix() const {return _T;}

  private:
    /// x-velocity matrix
    Matrix<double> _U;
    /// y-velocity matrix
    Matrix<double> _V;
    /// pressure matrix
    Matrix<double> _P;
    /// Temperature matrix
    Matrix<double> _T;
    /// x-momentum flux matrix
    Matrix<double> _F;
    /// y-momentum flux matrix
    Matrix<double> _G;
    /// right hand side matrix
    Matrix<double> _RS;

    /// kinematic viscosity
    double _nu;
    /// Heat diffusivity
    double _alpha;
    /// thermal expansion coefficient
    double _beta;
    /// gravitional accelearation in x direction
    double _gx{0.0};
    /// gravitional accelearation in y direction
    double _gy{0.0};
    /// timestep size
    double _dt;
    /// adaptive timestep coefficient
    double _tau;
    
};
