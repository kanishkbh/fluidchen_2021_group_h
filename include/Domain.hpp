#pragma once
#include "Enums.hpp"
#include <mpi.h>

/**
 * @brief Data structure that holds geometrical information
 * necessary for decomposition.
 *
 */
struct Domain {
    /// Minimum x index
    int imin{-1};
    /// Maximum x index
    int imax{-1};

    /// Minimum y index
    int jmin{-1};
    /// Maximum y index
    int jmax{-1};

    /// Cell length
    double dx{-1.0};
    /// Cell height
    double dy{-1.0};

    /// Number of cells in x direction
    int size_x{-1};
    /// Number of cells in y direction
    int size_y{-1};

    /// Number of cells in x direction, not-decomposed
    int domain_size_x{-1};
    /// Number of cells in y direction, not-decomposed
    int domain_size_y{-1};


    // Local data - pertains to a local processor

    int local_imin{-1};
    int local_imax{-1};

    int local_jmin{-1}; 
    int local_jmax{-1}; 
};
