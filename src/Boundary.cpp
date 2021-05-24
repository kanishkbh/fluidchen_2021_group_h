#include "Boundary.hpp"
#include <cmath>
#include <iostream>

// A fixed wall is just a moving wall with 0 velocity
FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : MovingWallBoundary(cells, 0) {}

//-----------------------------------------------------------------------------------------------------------
MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity)
    : _cells(cells), _velocity(wall_velocity) {}

void MovingWallBoundary::apply(Fields &field, bool pressure_only) {

    int i, j;
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
            // NB : pressure is taken as average ! => We set it to 0 then add two halves
            field.p(i, j) = 0;
            
            /// neighbor set classification:
            // what = 0 -> TR
            // what = 1 -> TL
            // what = 2 -> BR
            // what = 3 -> BL
            // what = 4 -> TB
            // what = 5 -> LR
            int what;
            if ( this_cell->is_border(border_position::TOP) && 
                  this_cell->is_border(border_position::RIGHT))
                what = 0;
            else if ( this_cell->is_border(border_position::TOP) && 
                       this_cell->is_border(border_position::LEFT))
                what = 1;
            else if ( this_cell->is_border(border_position::BOTTOM) && 
                       this_cell->is_border(border_position::RIGHT))
                what = 2;
            else if ( this_cell->is_border(border_position::BOTTOM) && 
                       this_cell->is_border(border_position::LEFT))
                what = 3;
            else if ( this_cell->is_border(border_position::TOP) && 
                       this_cell->is_border(border_position::BOTTOM))
                what = 4;
            else if ( this_cell->is_border(border_position::LEFT) && 
                       this_cell->is_border(border_position::RIGHT))
                what = 5;
            else
                throw std::runtime_error("Two-Borders classification error.");

            switch (what) {
            case 0: // TOP and RIGHT fluid neighbours
                if (!pressure_only) {
                    field.u(i, j) = 0;
                    field.v(i, j) = 0;
                    field.u(i-1, j) = -field.u(i-1,j+1);
                    field.v(i,j-1)  = -field.v(i+1,j-1);
                    field.g(i,j) = field.v(i,j);
                    field.f(i,j) = field.u(i,j);
                }
                field.p(i, j) = 0.5 * ( field.p(i,j+1) + field.p(i+1,j) );

                break;

            case 1: // TOP and LEFT fluid neighbours
                if (!pressure_only) {
                    field.u(i-1, j) = 0;
                    field.v(i, j) = 0;
                    field.u(i, j) = -field.u(i,j+1);
                    field.v(i,j-1)  = -field.v(i-1,j-1);
                    field.g(i,j) = field.v(i,j);
                    field.f(i,j) = field.u(i,j);
                }
                field.p(i,j) = 0.5 * ( field.p(i,j+1) + field.p(i+1,j) );

                break;

            case 2: // BOTTOM and RIGHT fluid
                if (!pressure_only) {
                    field.u(i-1, j) = 0;
                    field.v(i,j) = 0;
                    field.u(i,j) = -field.u(i,j+1);
                    field.v(i,j-1)  = -field.v(i-1,j-1);
                    field.g(i,j) = field.v(i,j);
                    field.f(i,j) = field.u(i,j);
                }
                field.p(i,j) = 0.5 * ( field.p(i,j+1) + field.p(i+1,j) );

                break;

            case 3: // BOTTOM and LEFT fluid
                if (!pressure_only) {
                    field.u(i,j) = 0;
                    field.v(i,j-1) = 0;
                    field.u(i-1,j) = -field.u(i-1,j-1);
                    field.v(i,j)  = -field.v(i+1,j);
                    field.g(i,j) = field.v(i,j);
                    field.f(i,j) = field.u(i,j);
                }
                field.p(i,j) = 0.5 * ( field.p(i,j-1) + field.p(i+1,j) );
                break;

            case 4: // TOP and BOTTOM fluid
                throw std::runtime_error("one cell thick wall not allowed !");
                break;

            case 5: // LEFT and RIGHT fluid
                throw std::runtime_error("one cell thick wall not allowed !"); 
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
            }
            
        }
        else if(this_cell->borders().size()==0){
            continue; 
        }
        else {
            throw std::runtime_error("Invalid BC : obstacle with invalid borders number");
        }
    }
}

//-----------------------------------------------------------------------------------------------------------
InflowBoundary::InflowBoundary(std::vector<Cell *> cells, double uin, double vin)
    : _cells(cells), _u_in(uin), _v_in(vin) {}

void InflowBoundary::apply(Fields &field, bool pressure_only) {

    int i, j;

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
        if(this_cell->borders().size() != 1) //Outflow only implemented with one fluid neighbour
            throw std::runtime_error("Outflow must only have one neighboring fluid cell.");
        //field.p(i, j) = _pressure;

        for (const auto &border : this_cell->borders()) {
            int i_n = this_cell->neighbour(border)->i();
            int j_n = this_cell->neighbour(border)->j();

            switch (border) {
            case border_position::BOTTOM:
                if (!pressure_only)
                    {
                        field.v(i, j - 1) = field.v(i, j - 2);
                        field.g(i, j - 1) = field.v(i, j - 1);
                    }
                field.p(i, j) = 2 * _pressure - field.p(i, j - 1);
                break;

            case border_position::TOP:
                if (!pressure_only)
                    {
                        field.v(i, j) = field.v(i, j + 1);
                        fielg.g(i, j) = field.v(i, j);
                    }
                field.p(i, j) = 2 * _pressure - field.p(i, j + 1);
                break;

            case border_position::LEFT:
                if (!pressure_only)
                    {
                        field.u(i - 1, j) = field.u(i - 2, j);
                        field.f(i - 1, j) = field.u(i - 1, j);
                    }
                field.p(i, j) = 2 * _pressure - field.p(i - 1, j);
                break;

            case border_position::RIGHT:
                if (!pressure_only)
                    {
                        field.u(i, j) = field.u(i + 1, j);
                        field.f(i, j) = field.u(i,j);
                    }
                field.p(i, j) = 2 * _pressure - field.p(i + 1, j);
                break;

            default:
                throw std::runtime_error("Unknown border type !");
                break;
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