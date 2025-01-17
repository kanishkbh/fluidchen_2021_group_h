#pragma once

#include <array>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "Cell.hpp"
#include "Datastructures.hpp"
#include "Domain.hpp"
#include "Enums.hpp"

/**
 * @brief Data structure holds cells and related sub-containers
 *
 */
class Grid {
  public:
    Grid() = default;

    /**
     * @brief Constructor for the Grid
     *
     * @param[in] geometry file name
     * @param[in] number of cells in x direction
     * @param[in] number of cells in y direction
     * @param[in] cell size in x direction
     * @param[in] cell size in y direction
     * @param[in] true if left cells are ghost cells
     * @param[in] true if right cells are ghost cells
     * @param[in] true if top cells are ghost cells
     * @param[in] true if bottom cells are ghost cells
     */
    Grid(std::string geom_name, Domain &domain, bool has_left_neighbor, bool has_right_neighbor, bool has_top_neighbor,
         bool has_bottom_neighbor);

    /// index based cell access
    Cell cell(int i, int j) const;

    /// access number of cells in x direction
    int imax() const;
    /// access number of cells in y direction
    int jmax() const;
    /// access number of cells in x direction including ghost cells
    int imaxb() const;
    /// access number of cells in x direction including ghost cells
    int jmaxb() const;

    /// access number of cells in x direction excluding ghost cells
    const Domain &domain() const;

    /// access cell size in x-direction
    double dx() const;
    /// access cell size in y-direction
    double dy() const;

    /**
     * @brief Access fluid cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &fluid_cells() const;

    /**
     * @brief Access inflow cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &inflow_cells() const;

    /**
     * @brief Access outflow cells
     *
     * @param[out] vector of fluid cells
     */
    const std::vector<Cell *> &outflow_cells() const;

    /**
     * @brief Access moving wall cells
     *
     * @param[out] vector of moving wall cells
     */
    const std::vector<Cell *> &moving_wall_cells() const;

    /**
     * @brief Access fixed wall cells
     *
     * @param[out] vector of fixed wall cells
     */
    const std::vector<Cell *> &fixed_wall_cells() const;

  private:
    /**@brief Default lid driven cavity case generator
     *
     * This function creates default lid driven cavity
     * case without need for a pgm file
     */
    void build_lid_driven_cavity();

    /// Build cell data structures with given geometrical data
    void assign_cell_types(std::vector<std::vector<int>> &geometry_data);
    /// Extract geometry from pgm file and create geometrical data
    void parse_geometry_file(std::string filedoc, std::vector<std::vector<int>> &geometry_data);

    Matrix<Cell> _cells;
    std::vector<Cell *> _fluid_cells; // Index 0 in PGM
    std::vector<Cell *> _inflow_cells; // Index 1 in PGM
    std::vector<Cell *> _outflow_cells; // Index 2 in PGM
    std::vector<Cell *> _fixed_wall_cells; // Index 3-7 in PGM
    std::vector<Cell *> _moving_wall_cells; // Index 8 in PGM

    Domain _domain;

    double _dx;
    double _dy;

    bool _has_left_neighbor, _has_right_neighbor, _has_top_neighbor, _has_bottom_neighbor;
};
