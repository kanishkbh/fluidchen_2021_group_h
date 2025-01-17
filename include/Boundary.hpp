#pragma once

#include <vector>

#include "Cell.hpp"
#include "Fields.hpp"

/**
 * @brief Abstact of boundary conditions.
 *
 * This class patches the physical values to the given field.
 */
class Boundary {
  public:
    /**
     * @brief Main method to patch the boundary conditons to given field and
     * grid
     *
     * @param[in] Field to be applied
     */
    virtual void apply(Fields &field, bool pressure_only = false) = 0;
    virtual ~Boundary() = default;
};



/**
 * @brief Moving wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities for the given velocity parallel to the fluid,
 * Neumann for pressure
 * Use an additional BC object for temperature
 */
class MovingWallBoundary : public Boundary {
  public:
    MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity);

    virtual ~MovingWallBoundary() = default;
    virtual void apply(Fields &field, bool pressure_only = false);

  private:
    std::vector<Cell *> _cells;
    double _velocity;
};

/**
 * @brief Fixed wall boundary condition for the outer boundaries of the domain.
 * Dirichlet for velocities, which is zero, Neumann for pressure
 * Use an additional BC object for temperature
 */
class FixedWallBoundary : public MovingWallBoundary {
  public:
    FixedWallBoundary(std::vector<Cell *> cells);
    virtual ~FixedWallBoundary() = default;

  private:
    std::vector<Cell *> _cells;
};

/**
 * @brief Inlet BC, with Dirichlet non-zero condition for normal velocity, zero for tangential velocity and
 * Neumann for pressure
 */
class InflowBoundary : public Boundary {
  public:
    InflowBoundary(std::vector<Cell *> cells, double u_in, double v_in);
    virtual ~InflowBoundary() = default;
    virtual void apply(Fields &field, bool pressure_only = false);

  private:
    std::vector<Cell *> _cells;
    double _u_in, _v_in;
};

/**
 * @brief Outlet BC, Neumann homogenous for normal velocity. Default Dirichlet pressure (see forum)
 */
class OutFlowBoundary : public Boundary {
  public:
    OutFlowBoundary(std::vector<Cell *> cells, double pressure);
    virtual ~OutFlowBoundary() = default;
    virtual void apply(Fields &field, bool pressure_only = false);

  private:
    std::vector<Cell *> _cells;
    double _pressure{0};
};

/**
 * @brief Fixed temperature BC
 */
class TemperatureDirichlet : public Boundary {
  public:
    TemperatureDirichlet(std::vector<Cell *> cells, double temp);
    virtual ~TemperatureDirichlet() = default;
    virtual void apply(Fields &field, bool pressure_only = false);

  private:
    std::vector<Cell *> _cells;
    double _temp{0};
};

/**
 * @brief Zero-temperature gradient BC
 */
class TemperatureAdiabatic : public Boundary {
  public:
    TemperatureAdiabatic(std::vector<Cell *> cells);
    virtual ~TemperatureAdiabatic() = default;
    virtual void apply(Fields &field, bool pressure_only = false);

  private:
    std::vector<Cell *> _cells;
};
