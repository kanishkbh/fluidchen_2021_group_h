#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

void FixedWallBoundary::apply(Fields &field) {

    for (int k=0; k < _cells.size(); k++) {
         /// check if it is top wall
         if (_cells[k]->is_border(BOTTOM) == true) {
            //  wall_id = 1; 
             //set y velocity of below cell to zero
             field.v(_cells[k]->i(), _cells[k]->j()-1) = 0;
             // set x velocity to 2*wall_velocity - u[bottom cell]
             field.u(_cells[k]->i(), _cells[k]->j()) = - field.u(_cells[k]->i(), _cells[k]->j()-1);
             // set Pressure to - P[below]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i(), _cells[k]->j()-1);
             // break;
         }
         else
         /// check if it is bottom wall
         if (_cells[i]->is_border(BOTTOM) == true) {
            //  wall_id = 2;
            // set y velocity to zero
             field.v(_cells[k]->i(), _cells[k]->j()) = 0;
             // set x velocity to 2*wall_velocity - u[above cell]
             field.u(_cells[k]->i(), _cells[k]->j()) = - field.u(_cells[k]->i(), _cells[k]->j()+1);
             // set Pressure to P[above]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i(), _cells[k]->j()+1);
              
            //  break;
         }
         else
         /// check if it is right wall
         if (_cells[i]->is_border(LEFT) == true) {
            //  wall_id = 3; 
            // set u[left cell] to zero
             field.u(_cells[k]->i()-1, _cells[k]->j()) = 0;
             // set v to 2*wall_velocity - v[left cell]
             field.v(_cells[k]->i(), _cells[k]->j()) =  - field.v(_cells[k]->i()-1, _cells[k]->j());
             // set Pressure to P[left]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i()-1, _cells[k]->j());
            //  break;
         }
         else
         /// check if it is left wall
         if (elem->is_border(RIGHT) == true) {
            //  wall_id = 4; 
            // set u[this cell] to zero
            field.u(_cells[k]->i(), _cells[k]->j()) = 0;
            // set v to 2*wall_velocity - v[right cell]
            field.v(_cells[k]->i(), _cells[k]->j()) =  - field.v(_cells[k]->i()+1, _cells[k]->j());
            // set Pressure to P[right]
            field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i()+1, _cells[k]->j());
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

     for (int k=0; k < _cells.size(); k++) {
         /// check if it is top wall
         if (_cells[k]->is_border(BOTTOM) == true) {
            //  wall_id = 1; 
             //set y velocity of below cell to zero
             field.v(_cells[k]->i(), _cells[k]->j()-1) = 0;
             // set x velocity to 2*wall_velocity - u[bottom cell]
             field.u(_cells[k]->i(), _cells[k]->j()) = 2*_wall_velocity - field.u(_cells[k]->i(), _cells[k]->j()-1);
             // set Pressure to - P[below]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i(), _cells[k]->j()-1);
             // break;
         }
         else
         /// check if it is bottom wall
         if (_cells[i]->is_border(BOTTOM) == true) {
            //  wall_id = 2;
            // set y velocity to zero
             field.v(_cells[k]->i(), _cells[k]->j()) = 0;
             // set x velocity to 2*wall_velocity - u[above cell]
             field.u(_cells[k]->i(), _cells[k]->j()) = 2*_wall_velocity - field.u(_cells[k]->i(), _cells[k]->j()+1);
             // set Pressure to P[above]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i(), _cells[k]->j()+1);
              
            //  break;
         }
         else
         /// check if it is right wall
         if (_cells[i]->is_border(LEFT) == true) {
            //  wall_id = 3; 
            // set u[left cell] to zero
             field.u(_cells[k]->i()-1, _cells[k]->j()) = 0;
             // set v to 2*wall_velocity - v[left cell]
             field.v(_cells[k]->i(), _cells[k]->j()) = 2*_wall_velocity - field.v(_cells[k]->i()-1, _cells[k]->j());
             // set Pressure to P[left]
             field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i()-1, _cells[k]->j());
            //  break;
         }
         else
         /// check if it is left wall
         if (elem->is_border(RIGHT) == true) {
            //  wall_id = 4; 
            // set u[this cell] to zero
            field.u(_cells[k]->i(), _cells[k]->j()) = 0;
            // set v to 2*wall_velocity - v[right cell]
            field.v(_cells[k]->i(), _cells[k]->j()) = 2*_wall_velocity - field.v(_cells[k]->i()+1, _cells[k]->j());
            // set Pressure to P[right]
            field.p(_cells[k]->i(), _cells[k]->j()) = field.p(_cells[k]->i()+1, _cells[k]->j());
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
