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
    virtual void init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {}

    /**
     * @brief Computes the residual in the current state.
     **/
    double res(Fields &field, Grid &grid) {
        double rloc = 0.0;

        for (auto currentCell : grid.fluid_cells()) {
            int i = currentCell->i();
            int j = currentCell->j();

            double val = Discretization::laplacian(field.p_matrix(), i, j) - field.rs(i, j);
            rloc += (val * val);
        }

        return rloc;
    }
};


/**
 * @brief Do N iterations of solver A, then switch to solver B
 * */
class MixedSolver : public PressureSolver {

  public:
    MixedSolver(std::unique_ptr<PressureSolver>&& f, std::unique_ptr<PressureSolver>&& s, int times_first)
        : first(std::move(f)), second(std::move(s)), max_iter(times_first) {}
    virtual ~MixedSolver() = default;

    /**
     * @brief Solve the pressure equation on given field, grid and boundary
     *
     * @param[in] field to be used
     * @param[in] grid to be used
     * @param[in] boundary to be used
     */
    virtual double solve(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
        if (iter_counter < max_iter) {
            ++iter_counter;
            auto res = first->solve(field, grid, boundaries);
            //std::cout << "Using first solver" << std::endl;
            if (iter_counter == max_iter) {
              second->init(field, grid, boundaries);
            //  std::cout << "Initiating second solver" << std::endl;
            }
            return res;
        } else {
          //std::cout << "Using second solver" << std::endl;
            return second->solve(field, grid, boundaries);
        }
    }

    /**
     * @brief For solvers that require an initialization step. Defaults to nothing.
     *
     * */
    virtual void init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) {
        iter_counter = 0;
        first->init(field, grid, boundaries);
    }

  std::unique_ptr<PressureSolver> first, second;
  int iter_counter = 0;
  int max_iter = 0;


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
    virtual void init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) override;

  private:
    Matrix<double> residual, direction, a_direction;
    double square_residual{0};
};

/**
 * @brief Steepest Descent
 *
 */
class SD : public PressureSolver {
  public:
    SD() = default;

    virtual ~SD() = default;

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
    virtual void init(Fields &field, Grid &grid, const std::vector<std::unique_ptr<Boundary>> &boundaries) override;

  private:
    Matrix<double> residual, a_residual;
};


/**
 * Richardon CG : 
 * This method is nothing exclusive, it is simply a CG method with unit preconditioner
 * Included for theoretical purposes 
 */ 
class CG_Richardson : public PressureSolver {
  public:
    CG_Richardson() = default; 
    virtual ~CG_Richardson() = default;
    virtual void init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) override; 
    virtual double solve(Fields& field,Grid& grid, const std::vector<std::unique_ptr<Boundary>>& boundaries);

  private: 
    Matrix<double> direction,adirection,residual; 
    double square_residual{0}; 
};

class CG_Jacobi : public PressureSolver {
  public:
    CG_Jacobi() = default; 
    virtual ~CG_Jacobi() = default;
    virtual void init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) override; 
    virtual double solve(Fields& field,Grid& grid, const std::vector<std::unique_ptr<Boundary>>& boundaries);

  private: 
    Matrix<double> direction,adirection,residual,cond_residual; 
    double sigma{0}; 
};

class CG_GS : public PressureSolver {
  public:
    CG_GS() = default; 
    virtual ~CG_GS() = default;
    virtual void init(Fields& field,Grid& grid,const std::vector<std::unique_ptr<Boundary>>& boundaries) override; 
    virtual double solve(Fields& field,Grid& grid, const std::vector<std::unique_ptr<Boundary>>& boundaries);

  private: 
    Matrix<double> direction,adirection,residual,cond_residual; 
    double sigma{0}; 
};