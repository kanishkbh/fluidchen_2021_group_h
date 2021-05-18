#include "Boundary.hpp"
#include <cmath>
#include <iostream>


//A fixed wall is just a moving wall with 0 velocity
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : MovingWallBoundary(cells, 0) {}

//-----------------------------------------------------------------------------------------------------------
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells), _velocity(wall_velocity) {}


void MovingWallBoundary::apply(Fields &field) {

    int i = 0, j = 0;
    auto w = _velocity;

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
