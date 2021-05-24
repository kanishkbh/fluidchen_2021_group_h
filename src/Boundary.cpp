#include "Boundary.hpp"
#include <cmath>
#include <iostream>

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells) : _cells(cells) {}

FixedWallBoundary::FixedWallBoundary(std::vector<Cell *> cells, std::map<int, double> wall_temperature)
    : _cells(cells), _wall_temperature(wall_temperature) {}

//-----------------------------------------------------------------------------------------------------------
void FixedWallBoundary::apply(Fields &field) {

    // for debugging
   //  std::cout << "\nInside FixedWallBoundary::apply\n"
   //              << "matrix imax = " << field.p_matrix().imax() << std::endl
   //              << "matrix jmax = " << field.p_matrix().jmax() << std::endl;
  
    int i = 0, j = 0;
//-----------------------------------------------------------------------------------------------------------
    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
//-----------------------------------------------------------------------------------------------------------
         // TOP-WALL CELL
         if (this_cell->is_border(border_position::BOTTOM)) {
            //  std::cout << "\nTOP WALL: i = " << i << ", j = " << j << std::endl;
            field.v(i, j-1) = 0; // CHEK: is it j-1 or j ?
            field.u(i, j) = - field.u(i, j-1); // u = - u[bottom]
            field.p(i,j) = field.p(i, j-1);
         }
//-----------------------------------------------------------------------------------------------------------
         // BOTTOM-WALL CELL
         if (this_cell->is_border(border_position::TOP)) {
            //  std::cout << "\nBOTTOM WALL: i = " << i << ", j = " << j << std::endl;
            field.v(i, j) = 0;
            field.u(i, j) = - field.u(i, j+1); // u = - u[top]
            field.p(i,j) = field.p(i,j+1);
         }
//-----------------------------------------------------------------------------------------------------------
         // RIGHT-WALL CELL
         if (this_cell->is_border(border_position::LEFT)) {
            //  std::cout << "\nRIGHT WALL: i = " << i << ", j = " << j << std::endl;
            field.u(i-1, j) = 0;
            field.v(i, j) = - field.v(i-1,j); // v = - v [left]
            field.p(i,j) = field.p(i-1,j);
         }
//-----------------------------------------------------------------------------------------------------------
         else
         // LEFT-WALL CELL
         if (this_cell->is_border(border_position::RIGHT)) {
            //  std::cout << "\nLEFT WALL: i = " << i << ", j = " << j << std::endl;
            field.u(i,j) = 0;
            field.v(i,j) = - field.v(i+1,j); // v = - v[right]
            field.p(i,j) = field.p(i+1,j);
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

    // FOR DEBUGGING
   //  std::cout << "\nInside MovingWallBoundary::apply\n"
   //          << "matrix imax = " << field.p_matrix().imax() << std::endl
   //          << "matrix jmax = " << field.p_matrix().jmax() << std::endl
   //          << "wall velocity = " << _wall_velocity[LidDrivenCavity::moving_wall_id] << std::endl;

    int i = 0, j = 0;
    auto w =  _wall_velocity[LidDrivenCavity::moving_wall_id];

    /// cycle through all cells
    for (auto this_cell : _cells) {
        i = this_cell->i();
        j = this_cell->j();
//-----------------------------------------------------------------------------------------------------------
         // TOP-WALL CELL
         if (this_cell->is_border(border_position::BOTTOM)) {
            //  std::cout << "\nTOP WALL: i = " << i << ", j = " << j << std::endl
            //             << "Top wall cell-indices should be : i," << field.p_matrix().jmax()-1 << std::endl;
            field.v(i, j-1) = 0; // CHEK: is it j-1 or j ?
            field.u(i, j) = 2*w - field.u(i, j-1); // u = 2*wall_velocity - u[bottom]
            field.p(i,j) = field.p(i, j-1);
         }
//-----------------------------------------------------------------------------------------------------------
         // BOTTOM-WALL CELL
         if (this_cell->is_border(border_position::TOP)) {
            //  std::cout << "\nBOTTOM WALL: i = " << i << ", j = " << j << std::endl;
            field.v(i, j) = 0;
            field.u(i, j) = 2*w - field.u(i, j+1); // u = - u[top]
            field.p(i,j) = field.p(i,j+1);
         }
//-----------------------------------------------------------------------------------------------------------
         // RIGHT-WALL CELL
         if (this_cell->is_border(border_position::LEFT)) {
            //  std::cout << "\nRIGHT WALL: i = " << i << ", j = " << j << std::endl;
            field.u(i-1, j) = 0;
            field.v(i, j) = 2*w - field.v(i-1,j); // v = - v [left]
            field.p(i,j) = field.p(i-1,j);
         }
//-----------------------------------------------------------------------------------------------------------
         else
         // LEFT-WALL CELL
         if (this_cell->is_border(border_position::RIGHT)) {
            //  std::cout << "\nLEFT WALL: i = " << i << ", j = " << j << std::endl;
            field.u(i,j) = 0;
            field.v(i,j) = 2*w - field.v(i+1,j); // v = - v[right]
            field.p(i,j) = field.p(i+1,j);
         }
    
    }

}
//-----------------------------------------------------------------------------------------------------------
/*
// APPLY FUNCTION USING CELL::NEIGHBOUR FUNCTION - DIDN'T WORK
// LEFT HERE TO DEBUG LATER

     for (auto this_cell : _cells) {
         /// check if it is top wall
         if (this_cell->is_border(border_position::BOTTOM) == true) {
             auto bottom_cell = this_cell->neighbour(border_position::BOTTOM);
             //set v of below cell to zero
             field.v(bottom_cell->i(), bottom_cell->j()) = 0;
             // set u to 2*wall_velocity - u[below cell]
             field.u(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(bottom_cell->i(), bottom_cell->j());
             // set Pressure = P[below]
            //  field.p(this_cell->i(), this_cell->j()) = field.p(bottom_cell->i(), bottom_cell->j());
         }
         else
         /// check if it is bottom wall
         if (this_cell->is_border(border_position::TOP) == true) {
            auto top_cell = this_cell->neighbour(border_position::TOP);
            // set v to zero
             field.v(this_cell->i(), this_cell->j()) = 0;
             // set u to 2*wall_velocity - u[top cell]
             field.u(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.u(top_cell->i(), top_cell->j());
             // set Pressure to P[top]
            //  field.p(this_cell->i(), this_cell->j()) = field.p(top_cell->i(), top_cell->j());
         }
         else
         /// check if it is right wall
         if (this_cell->is_border(border_position::LEFT) == true) {
            auto left_cell = this_cell->neighbour(border_position::LEFT);
            // set u[left cell] to zero
             field.u(this_cell->i()-1, this_cell->j()) = 0;
             // set v to 2*wall_velocity - v[left cell]
             field.v(this_cell->i(), this_cell->j()) = 2*_wall_velocity[LidDrivenCavity::moving_wall_id] - field.v(this_cell->i()-1, this_cell->j());
             // set Pressure to P[left]
            //  field.p(this_cell->i(), this_cell->j()) = field.p(left_cell->i(), left_cell->j());
         }
         else
         /// check if it is left wall
         if (this_cell->is_border(border_position::RIGHT) == true) {
            auto right_cell = this_cell->neighbour(border_position::RIGHT);
            // set u[this cell] to zero
            field.u(this_cell->i(), this_cell->j()) = 0;
            // set v to 2*wall_velocity - v[right cell]
            field.v(this_cell->i(), this_cell->j()) = 2*_wall_velocity [LidDrivenCavity::moving_wall_id]- field.v(right_cell->i(), right_cell->j());
            // set Pressure to P[right]
            // field.p(this_cell->i(), this_cell->j()) = field.p(right_cell->i(), right_cell->j());
         }
     }
*/

/*
/// FOR REFERENCE: HARD CODED BOUNDARY CONDITIONS ASSUMING TOPWALL MOVING AND OTHERS FIXED.
/// (WROTE THIS WHEN APPLY FUNCN WASN'T WORKING)

    // For the left and right walls
    for(int j=1;j<=jmax;++j){
         
        // Left wall : U = 0 on ghost cell. 
        // Right wall : U(imax) = 0 on pre-ghost cell 
        _U(0,j) = 0;
        _U(imax,j) = 0;
        // Left wall :0.5*( V(0,j) + V(1,j)) = 0;  
        _V(0,j) = -_V(1,j);
        // Right wall :0.5*( V(imax+1,j) + V(imax,j)) = 0; 
        _V(imax+1,j) = -_V(imax,j);
    }
    
    // For the top and bottom walls
    for(int i=1;i<=imax;++i){
        _U(i,0) = -_U(i,1);
        _U(i,jmax+1) = 2-_U(i,jmax);
        _V(i,0) = 0; 
        _V(i,jmax) = 0;
    }
*/