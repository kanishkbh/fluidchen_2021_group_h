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

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT,
};

namespace MessageTag {
    enum {
        DOMAIN,
        F,
        G,
        T,
        P,
        U,
        V,
        TIMESTEP
    };

}

enum GeometryType {
    obstacle_id,
    fluid_outlet,
    fluid_inlet_v,
    fluid_inlet_u,
    inlet_temp,
    fluid_interior,
};


namespace border {
const int TOP = 0;
const int BOTTOM = 1;
const int LEFT = 2;
const int RIGHT = 3;
} // namespace border

enum class cell_type {

    FLUID,
    INFLOW,
    OUTFLOW,
    FIXED_WALL_3,
    FIXED_WALL_4,
    FIXED_WALL_5,
    FIXED_WALL_6,
    FIXED_WALL_7,
    MOVING_WALL,
    MPI_GHOST,
    DEFAULT,
};
