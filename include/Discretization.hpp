#pragma once

#include "Datastructures.hpp"

/**
 * @brief Static discretization methods to modify the fields
 *
 */
class Discretization {
  public:
    Discretization() = default;

    /**
     * @brief Constructor to set the discretization parameters
     *
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     * @param[in] upwinding coefficient
     */
    Discretization(double dx, double dy, double gamma);

    /**
     * @brief Convection in x direction using donor-cell scheme
     * (with the global value of _gamma, from the last constructor call)
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_u(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Convection in y direction using donor-cell scheme
     * (with the global value of _gamma, from the last constructor call)
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_v(const Matrix<double> &U, const Matrix<double> &V, int i, int j);

    /**
     * @brief Heat-Convection in x direction (i.e. d(uT)/dx) using donor-cell scheme
     * (with the global value of _gamma, from the last constructor call)
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_u_T(const Matrix<double> &U, const Matrix<double> &T, int i, int j);

    /**
     * @brief Heat-Convection in y direction (i.e. d(vT)/dy) using donor-cell scheme 
     * (with the global value of _gamma, from the last constructor call)
     *
     * @param[in] x-velocity field
     * @param[in] y-velocity field
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double convection_v_T(const Matrix<double> &V, const Matrix<double> &T, int i, int j);

    /**
     * @brief Laplacian term discretization using central difference
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double laplacian(const Matrix<double> &P, int i, int j);

    /**
     * @brief Terms of laplacian needed for SOR, i.e. excluding unknown value at
     * (i,j)
     *
     * @param[in] data to be discretized
     * @param[in] x index
     * @param[in] y index
     * @param[out] result
     *
     */
    static double sor_helper(const Matrix<double> &P, int i, int j);


  private:
    static double _dx;
    static double _dy;
    static double _gamma;
};
