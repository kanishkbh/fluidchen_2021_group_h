#pragma once
#include <array>
#include "Enums.hpp"


/**
 * @brief An element of the matrix of processors that holds rank and neighbouring procssors.
 *
 */
class Processor {
    
    private: 

        /// x index
        int _ip{0};
        /// y index
        int _jp{0};
        /// Processor id (rank)
        int _id{0};

        /// Vector of bools that holds neighbor conditions. TOP-BOTTOM-LEFT-RIGHT
        std::array<bool, 4> _neighbour_bool{false, false, false, false};
        
        /// Pointers to neighbours. // TOP -  BOTTOM - LEFT - RIGHT
        std::array<Processor *, 4> _neighbours;

    public: 

        /**
        * @brief Constructor for Processor object
        *
        * @param[in] i x index of the processor
        * @param[in] j y index of the processor
        * @param[in] id of the processor
        */
        Processor(int i, int j, int id);

        /**
        * @brief Neighbour getter for given positon
        *
        * @param[in] border position
        * @param[out] pointer to the neighboring cell
        */
        const Processor* neighbour(border_position position) const;

        /**
        * @brief Neighbour setter for given positon
        *
        * @param[in] border position
        */
        void set_neighbour(Processor *processor, border_position position);

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
        int wall_id() const;
};

// Ip and Jp getters 

std::array<unsigned,2> coordinates(int rank) {
    std::array<unsigned int,2> coord;
    // coord houses ip and jp respectively 
    coord[0] = rank /
}