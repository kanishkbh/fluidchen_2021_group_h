#pragma once

#include "Datastructures.hpp"
#include <vector>
#include <mpi.h>

class Communication {

    public:
    /**
     * @brief Sends the right-part of x to the appropriate neighbor (if any) and read from the left the new value.
     * Tag should identify what field is transmitted (U, V, P, ...)
     * */

    static void communicate_right(Matrix<double>& x, int tag);

    /**
     * @brief Sends the left-part of x to the appropriate neighbor (if any) and read from the right the new value.
     * Tag should identify what field is transmitted (U, V, P, ...)
     * */
    static void communicate_left(Matrix<double>& x, int tag);

    /**
     * @brief Sends the top-part of x to the appropriate neighbor (if any) and read from the bottom the new value.
     * Tag should identify what field is transmitted (U, V, P, ...)
     * */
    static void communicate_top(Matrix<double>& x, int tag);

    /**
     * @brief Sends the bottom-part of x to the appropriate neighbor (if any) and read from the top the new value.
     * Tag should identify what field is transmitted (U, V, P, ...)
     * */
    static void communicate_bottom(Matrix<double>& x, int tag);

    // Threads of the neighboring
    static int _left_neighbor_rank, _right_neighbor_rank, _top_neighbor_rank, _bottom_neighbor_rank;

};