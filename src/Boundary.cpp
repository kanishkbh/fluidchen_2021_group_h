#include "Boundary.hpp"
#include <cmath>
#include <iostream>

// A fixed wall is just a moving wall with 0 velocity
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : MovingWallBoundary(cells, 0) {}

//-----------------------------------------------------------------------------------------------------------
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity)
    : _cells(cells), _velocity(wall_velocity) {}

void MovingWallBoundary::apply(Fields &field, bool pressure_only) {

    int i = 0, j = 0;
    auto w = _velocity;

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();

        // Different regime depending on the number of fluid neighbors. If there are 2, one is vertical and one
        // horizontal
        if (this_cell->borders().size() == 1) {
            auto border = this_cell->borders()[0];

            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                if (!pressure_only) {
                    field.v(i, j - 1) = 0;
                    field.u(i, j) = 2 * w - field.u(i, j - 1); // Average of velocities is w
                    field.g(i, j - 1) = field.v(i, j - 1);
                }

                field.p(i, j) = field.p(i, j - 1);

                break;

            case border_position::TOP:
                if (!pressure_only) {
                    field.v(i, j) = 0;
                    field.u(i, j) = 2 * w - field.u(i, j + 1); // Average of velocities is w
                    field.g(i, j) = field.v(i, j);
                }
                field.p(i, j) = field.p(i, j + 1);
                break;

            case border_position::LEFT:
                if (!pressure_only) {
                    field.u(i - 1, j) = 0;
                    field.v(i, j) = 2 * w - field.v(i - 1, j); // 0.5*(v + v [left]) = w
                    field.f(i, j) = field.u(i - 1, j);
                }

                field.p(i, j) = field.p(i - 1, j);
                break;

            case border_position::RIGHT:
                if (!pressure_only) {
                    field.u(i, j) = 0;
                    field.v(i, j) = 2 * w - field.v(i + 1, j); // v = - v[right]
                    field.f(i, j) = field.u(i, j);
                }

                field.p(i, j) = field.p(i + 1, j);

                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
        } else if (this_cell->borders().size() == 2) {
            // NB : pressure is taken as average !
            field.p(i, j) = 0;
            for (const auto &border : this_cell->borders()) {
                int i_n = this_cell->neighbour(border)->i();
                int j_n = this_cell->neighbour(border)->j();

                switch (border) {
                case border_position::BOTTOM:
                    if (!pressure_only) {
                        field.v(i, j - 1) = 0;
                        field.g(i, j - 1) = field.v(i, j - 1);
                    }

                    field.p(i, j) += 0.5 * field.p(i, j - 1);

                    break;

                case border_position::TOP:
                    if (!pressure_only) {
                        field.v(i, j) = 0;
                        field.g(i, j) = field.v(i, j);
                    }
                    field.p(i, j) += 0.5 * field.p(i, j + 1);

                    break;

                case border_position::LEFT:
                    if (!pressure_only) {
                        field.u(i - 1, j) = 0;

                        field.f(i - 1, j) = field.u(i - 1, j);
                    }
                    field.p(i, j) += 0.5 * field.p(i - 1, j);
                    break;

                case border_position::RIGHT:
                    if (!pressure_only) {
                        field.u(i, j) = 0;
                        field.f(i, j) = field.u(i, j);
                    }
                    field.p(i, j) += 0.5 * field.p(i + 1, j);

                    break;

                default:
                    throw std::runtime_error("Unknown border type !");
                    break;
                }
            }
        } else {
            throw std::runtime_error("Invalid BC : obstacle with invalid borders number");
        }
    }
}

//-----------------------------------------------------------------------------------------------------------
InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double uin, double vin)
    : _cells(cells), _u_in(uin), _v_in(vin) {}

void InflowBoundary::apply(Fields &field, bool pressure_only) {

    int i = 0, j = 0;

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
        for (const auto &border : this_cell->borders()) {
            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                if (!pressure_only) {
                    field.v(i, j - 1) = _v_in;
                    field.u(i, j) = 2 * _u_in - field.u(i, j - 1);
                    field.g(i, j - 1) = field.v(i, j - 1);
                }
                field.p(i, j) = field.p(i, j - 1);

                break;

            case border_position::TOP:
                if (!pressure_only) {
                    field.v(i, j) = _v_in;
                    field.u(i, j) = 2 * _u_in - field.u(i, j + 1);
                    field.g(i, j) = field.v(i, j);
                }
                field.p(i, j) = field.p(i, j + 1);

                break;

            case border_position::LEFT:
                if (!pressure_only) {
                    field.u(i - 1, j) = _u_in;
                    field.v(i, j) = 2 * _v_in - field.v(i - 1, j);
                    field.f(i - 1, j) = field.u(i - 1, j);
                }
                field.p(i, j) = field.p(i - 1, j);

                break;

            case border_position::RIGHT:
                if (!pressure_only) {
                    field.v(i, j) = 2 * _v_in - field.v(i + 1, j);
                    field.u(i, j) = _u_in;
                    field.f(i, j) = field.u(i, j);
                }
                field.p(i, j) = field.p(i + 1, j);
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
        }
    }
}

OutFlowBoundary::OutFlowBoundary(std::vector<Cell *> cells, double pressure) : _cells(cells), _pressure(pressure) {}

void OutFlowBoundary::apply(Fields &field, bool pressure_only) {

    int i = 0, j = 0;

    /// cycle through all cells

    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
        field.p(i, j) = _pressure;

        if (!pressure_only) {
            for (const auto &border : this_cell->borders()) {
                int i_n = this_cell->neighbour(border)->i();
                int j_n = this_cell->neighbour(border)->j();

                switch (border) {
                case border_position::BOTTOM:
                    field.v(i, j - 1) = field.v(i, j - 2);
                    break;

                case border_position::TOP:
                    field.v(i, j) = field.v(i, j + 1);

                    break;

                case border_position::LEFT:
                    field.u(i - 1, j) = field.u(i - 2, j);
                    break;

                case border_position::RIGHT:
                    field.u(i, j) = field.u(i + 1, j);
                    break;

                default:
                    throw std::runtime_error("Unknown border type !");
                    break;
                }
            }
        }
    }
}

TemperatureDirichlet::TemperatureDirichlet(std::vector<Cell *> cells, double temp) : _cells(cells), _temp(temp) {}

void TemperatureDirichlet::apply(Fields &field, bool pressure_only) {

    if (pressure_only)
        return;

    int i, j;
    /// cycle through all cells

    for (auto current_cell : _cells) {
        i = current_cell->i();
        j = current_cell->j();

        if (current_cell->borders().size() == 1) {
            const auto border = current_cell->borders()[0];

            switch (border) {

            case border_position::TOP:
                field.t(i, j) = 2 * _temp - field.t(i, j + 1);
                break;
            case border_position::BOTTOM:
                field.t(i, j) = 2 * _temp - field.t(i, j - 1);
                break;
            case border_position::LEFT:
                field.t(i, j) = 2 * _temp - field.t(i - 1, j);
                break;
            case border_position::RIGHT:
                field.t(i, j) = 2 * _temp - field.t(i + 1, j);
                break;
            default:
                std::runtime_error("Invalid border position (1) : @ Dirichlet BC");
                break;
            }
        } else if (current_cell->borders().size() == 2) {
            // NB : Temperature is taken as average => we set it to 0 then add 0.5 each time
            for (auto border : current_cell->borders()) {
                field.t(i, j) = 0;
                switch (border) {

                case border_position::TOP:
                    field.t(i, j) += _temp - 0.5 * field.t(i, j + 1);
                    break;

                case border_position::BOTTOM:
                    field.t(i, j) += _temp - 0.5 * field.t(i, j - 1);
                    break;
                case border_position::LEFT:
                    field.t(i, j) += _temp - 0.5 * field.t(i - 1, j);
                    break;
                case border_position::RIGHT:
                    field.t(i, j) += _temp - 0.5 * field.t(i + 1, j);
                    break;
                default:
                    std::runtime_error("Invalid border position (2) : @ Dirichlet BC");
                }
            }
        }

        else if (current_cell->borders().size() == 0) {
            continue;
        }

        else {
            std::runtime_error("Invalid Boundary Conditions @ TemperatureDirichlet");
        }
    }
}

TemperatureAdiabatic::TemperatureAdiabatic(std::vector<Cell *> cells) : _cells(cells) {}

void TemperatureAdiabatic::apply(Fields &field, bool pressure_only) {

    if (pressure_only)
            return;
    int i, j;

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();

        // Different regime depending on the number of fluid neighbors. If there are 2, one is vertical and one
        // horizontal
        if (this_cell->borders().size() == 1) {
            auto border = this_cell->borders()[0];

            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                field.t(i, j) = field.t(i, j - 1);
                break;

            case border_position::TOP:
                field.t(i, j) = field.t(i, j + 1);
                break;

            case border_position::LEFT:
                field.t(i, j) = field.t(i - 1, j);
                break;

            case border_position::RIGHT:
                field.t(i, j) = field.t(i + 1, j);
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
        } else if (this_cell->borders().size() == 2) {
            // NB : Temperature is taken as average => we set it to 0 then add 0.5 each time
            field.t(i, j) = 0;
            for (const auto &border : this_cell->borders()) {
                int i_n = this_cell->neighbour(border)->i();
                int j_n = this_cell->neighbour(border)->j();

                switch (border) {
                case border_position::BOTTOM:
                    field.t(i, j) += 0.5 * field.t(i, j - 1);
                    break;

                case border_position::TOP:
                    field.t(i, j) += 0.5 * field.t(i, j + 1);
                    break;

                case border_position::LEFT:
                    field.t(i, j) += 0.5 * field.t(i - 1, j);
                    break;

                case border_position::RIGHT:
                    field.t(i, j) += 0.5 * field.t(i + 1, j);
                    break;

                default:
                    throw std::runtime_error("Unknown border type !");
                    break;
                }
            }
        } else {
            throw std::runtime_error("Invalid BC : obstacle with invalid borders number");
        }
    }
}