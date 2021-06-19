#pragma once

#include "Boundary.hpp"
#include "Fields.hpp"
#include "Grid.hpp"
#include <utility>
/**
 * @brief Abstract class for pressure Poisson equation solver
 *
 */
class PressureSolver {
  public:
    PressureSolver() = default;
    virtual ~PressureSolver() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) = 0;

  /**
   * @brief For solvers that require an initialization step. Defaults to nothing.
   * 
   * */
    virtual double init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {}
};

/**
 * @brief Successive Over-Relaxation algorithm for solution of pressure Poisson
 * equation
 *
 */
class SOR : public PressureSolver {
  public:
    SOR() = default;

    /**
     * @brief Constructor of SOR solver
     *
     * @param[in] relaxation factor
     */
    SOR(double omega);

    virtual ~SOR() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

  private:
    double _omega;
};

/**
 * @brief Conjugate gradients algorithm
 *
 */
class CG : public PressureSolver {
  public:
    CG() = default;

    virtual ~CG() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries);

    /**
     * @brief Initializes the CG algorithm
     * */
    virtual double init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) override;

  private:
    Matrix<double> residual, direction;
    double square_residual{0};
};
