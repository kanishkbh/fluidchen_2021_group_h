#pragma once

// If no geometry file is provided in the input file, lid driven cavity case
// will run by default. In the Grid.cpp, geometry will be created following
// PGM convention, which is:
// 0: fluid, 3: fixed wall, 4: moving wall
namespace LidDrivenCavity {
const int moving_wall_id = 8;
const int fixed_wall_id = 4;
const double wall_velocity = 1.0;
} // namespace LidDrivenCavity

namespace GeometryType {
    const int moving_wall = 8; 
    const int cooled_wall = 5; 
    const int heated_wall = 4;
    const int insulated_wall = 3;
    const int fluid_outlet = 2;
    const int fluid_inlet = 1;
    const int fluid_interior = 0;
}

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {
    FLUID,          // WS1: fluid interior = 0 in .pgm file
    FIXED_WALL,
    MOVING_WALL,
    DEFAULT
};
