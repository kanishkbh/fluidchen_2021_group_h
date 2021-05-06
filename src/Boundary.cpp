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
