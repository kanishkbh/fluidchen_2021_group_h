#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

//-----------------------------------------------------------------------------------------------------------
void FixedWallBoundary::apply(Fields &field) {

    int i = 0, j = 0;
    auto w = 0;

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
        for (const auto &border : this_cell->borders()) {
            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                field.v(i, j - 1) = 0;
                field.u(i, j) = 2 * w - field.u(i, j - 1); // Average of velocities is w
                field.p(i, j) = field.p(i, j - 1);
                break;

            case border_position::TOP:
                field.v(i, j) = 0;
                field.u(i, j) = 2 * w - field.u(i, j + 1); // Average of velocities is w
                field.p(i, j) = field.p(i, j + 1);
                break;

            case border_position::LEFT:
                field.u(i - 1, j) = 0;
                field.v(i, j) = 2 * w - field.v(i - 1, j); // 0.5*(v + v [left]) = w
                field.p(i, j) = field.p(i - 1, j);
                break;

            case border_position::RIGHT:
                field.u(i, j) = 0;
                field.v(i, j) = 2 * w - field.v(i + 1, j); // v = - v[right]
                field.p(i, j) = field.p(i + 1, j);
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
        }
    }
}
//-----------------------------------------------------------------------------------------------------------
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}
//-----------------------------------------------------------------------------------------------------------
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}
//-----------------------------------------------------------------------------------------------------------
void MovingWallBoundary::apply(Fields &field) {

    int i = 0, j = 0;
    auto w = _wall_velocity[LidDrivenCavity::moving_wall_id];

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
        for (const auto &border : this_cell->borders()) {
            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                field.v(i, j - 1) = 0;
                field.u(i, j) = 2 * w - field.u(i, j - 1); // Average of velocities is w
                field.p(i, j) = field.p(i, j - 1);
                break;

            case border_position::TOP:
                field.v(i, j) = 0;
                field.u(i, j) = 2 * w - field.u(i, j + 1); // Average of velocities is w
                field.p(i, j) = field.p(i, j + 1);
                break;

            case border_position::LEFT:
                field.u(i - 1, j) = 0;
                field.v(i, j) = 2 * w - field.v(i - 1, j); // 0.5*(v + v [left]) = w
                field.p(i, j) = field.p(i - 1, j);
                break;

            case border_position::RIGHT:
                field.u(i, j) = 0;
                field.v(i, j) = 2 * w - field.v(i + 1, j); // v = - v[right]
                field.p(i, j) = field.p(i + 1, j);
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
        }
    }
}
