#pragma once
#include <array>
#include "Enums.hpp"


/**
 * @brief An element of the matrix of processors that holds rank and neighbouring procssors.
 */

class Processor {
    
    private: 
        /// total number of processors in x and y direction
        int _iproc{1}, _jproc{1};
        /// x index of this processor
        int _ip{0};
        /// y index of this processor
        int _jp{0};
        /// Processor id (rank)
        int _id{0};

        /// Vector of bools that holds neighbor conditions. TOP-BOTTOM-LEFT-RIGHT
        std::array<bool, 4> _neighbours_bool{false, false, false, false};
        
        /// ID if neighbour processors. // TOP,BOTTOM,LEFT,RIGHT
        std::array<int, 4> _neighbours;
        
        /// Neigbour and neighbour_bool setter (called by constructor)
        void set_neighbours();

    public: 

        /**
        * @brief Constructor for Processor object
        *
        * @param[in] iproc number of processors in x direction
        * @param[in] jproc number of processors in y direction
        * @param[in] rank of this processor
        */
        Processor(int iproc, int jproc, int rank);

        /**
        * @brief Neighbour getter for given positon
        *
        * @param[in] border_position at which neighbour id is required
        * @param[out] id of the neighbouring processor
        */
        const int neighbour(border_position position) const;

        /**
        * @brief Check whether the given position (TOP/BOTTOM/LEFT/RIGHT) 
        * has a neighbour
        *
        * @param[in] border position where neigbour-processor may exist
        * @param[out] whether the given position has a neighbouring processor
        */
        bool has_neighbour(border_position position) const;

        /// Getter of x index
        int ip() const;
        /// Getter of y index
        int jp() const;
        /// Getter of processor id
        int proc_id() const;
};