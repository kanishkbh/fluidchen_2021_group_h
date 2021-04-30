#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    for (auto this_cell : _cells) {
         /// check if it is top wall
         if (this_cell->is_border(border_position::BOTTOM) == true) {
            //  wall_id = 1; 
             auto bottom_cell = this_cell->neighbour(border_position::BOTTOM);
             //set v of below cell to zero
             field.v(bottom_cell->i(), bottom_cell->j()) = 0;
             // set u to - u[below cell]
             field.u(this_cell->i(), this_cell->j()) = - field.u(bottom_cell->i(), bottom_cell->j());
             // set Pressure = P[below]
             field.p(this_cell->i(), this_cell->j()) = field.p(bottom_cell->i(), bottom_cell->j());
             // break;
         }
         else
         /// check if it is bottom wall
         if (this_cell->is_border(border_position::TOP) == true) {
            //  wall_id = 2;
            auto top_cell = this_cell->neighbour(border_position::TOP);
            // set v to zero
             field.v(this_cell->i(), this_cell->j()) = 0;
             // set u to - u[top cell]
             field.u(this_cell->i(), this_cell->j()) = - field.u(top_cell->i(), top_cell->j());
             // set Pressure to P[top]
             field.p(this_cell->i(), this_cell->j()) = field.p(top_cell->i(), top_cell->j());
              
            //  break;
         }
         else
         /// check if it is right wall
         if (this_cell->is_border(border_position::LEFT) == true) {
            //  wall_id = 3; 
            auto left_cell = this_cell->neighbour(border_position::LEFT);
            // set u[left cell] to zero
             field.u(this_cell->i()-1, this_cell->j()) = 0;
             // set v to - v[left cell]
             field.v(this_cell->i(), this_cell->j()) = - field.v(this_cell->i()-1, this_cell->j());
             // set Pressure to P[left]
             field.p(this_cell->i(), this_cell->j()) = field.p(left_cell->i(), left_cell->j());
            //  break;
         }
         else
         /// check if it is left wall
         if (this_cell->is_border(border_position::RIGHT) == true) {
            //  wall_id = 4; 
            auto right_cell = this_cell->neighbour(border_position::RIGHT);
            // set u[this cell] to zero
            field.u(this_cell->i(), this_cell->j()) = 0;
            // set v to - v[right cell]
            field.v(this_cell->i(), this_cell->j()) = - field.v(right_cell->i(), right_cell->j());
            // set Pressure to P[right]
            field.p(this_cell->i(), this_cell->j()) = field.p(right_cell->i(), right_cell->j());
            //  break;
         }
     }


}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, double wall_velocity) : _cells(cells) {
    _wall_velocity.insert(std::pair(LidDrivenCavity::moving_wall_id, wall_velocity));
}

MovingWallBoundary::MovingWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_velocity,
                                       std::map<int, double> wall_temperature)
    : _cells(cells), _wall_velocity(wall_velocity), _wall_temperature(wall_temperature) {}

void MovingWallBoundary::apply(Fields &field) {

    // wall_id = 1 for top, 2 for bottom, 3 for left, 4 for right
    // int wall_id = 0;

     for (auto this_cell : _cells) {
         /// check if it is top wall
         if (this_cell->is_border(border_position::BOTTOM) == true) {
            //  wall_id = 1; 
             auto bottom_cell = this_cell->neighbour(border_position::BOTTOM);
             //set v of below cell to zero
             field.v(bottom_cell->i(), bottom_cell->j()) = 0;
             // set u to 2*wall_velocity - u[below cell]
             field.u(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(bottom_cell->i(), bottom_cell->j());
             // set Pressure = P[below]
             field.p(this_cell->i(), this_cell->j()) = field.p(bottom_cell->i(), bottom_cell->j());
             // break;
         }
         else
         /// check if it is bottom wall
         if (this_cell->is_border(border_position::TOP) == true) {
            //  wall_id = 2;
            auto top_cell = this_cell->neighbour(border_position::TOP);
            // set v to zero
             field.v(this_cell->i(), this_cell->j()) = 0;
             // set u to 2*wall_velocity - u[top cell]
             field.u(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(top_cell->i(), top_cell->j());
             // set Pressure to P[top]
             field.p(this_cell->i(), this_cell->j()) = field.p(top_cell->i(), top_cell->j());
              
            //  break;
         }
         else
         /// check if it is right wall
         if (this_cell->is_border(border_position::LEFT) == true) {
            //  wall_id = 3; 
            auto left_cell = this_cell->neighbour(border_position::LEFT);
            // set u[left cell] to zero
             field.u(this_cell->i()-1, this_cell->j()) = 0;
             // set v to 2*wall_velocity - v[left cell]
             field.v(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(this_cell->i()-1, this_cell->j());
             // set Pressure to P[left]
             field.p(this_cell->i(), this_cell->j()) = field.p(left_cell->i(), left_cell->j());
            //  break;
         }
         else
         /// check if it is left wall
         if (this_cell->is_border(border_position::RIGHT) == true) {
            //  wall_id = 4; 
            auto right_cell = this_cell->neighbour(border_position::RIGHT);
            // set u[this cell] to zero
            field.u(this_cell->i(), this_cell->j()) = 0;
            // set v to 2*wall_velocity - v[right cell]
            field.v(this_cell->i(), this_cell->j()) = 2*_wall_velocity [LidDrivenCavity::moving_wall_id]- field.v(right_cell->i(), right_cell->j());
            // set Pressure to P[right]
            field.p(this_cell->i(), this_cell->j()) = field.p(right_cell->i(), right_cell->j());
            //  break;
         }
     }

    //  /// DO stuff
    //  switch (wall_id) {

    //      case 1: { // Top wall is moving
    //          for (int i = 1; i < cells.size(); ++i) {

    //          }
    //      }
    //  }
}
