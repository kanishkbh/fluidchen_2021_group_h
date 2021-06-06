#pragma once
#include <array>
#include "Enums.hpp"
#include "Fields.hpp"
#include "Datastructures.hpp"

/**
 * @brief A set of C-style arrays to serve as buffer for communicating field values
 * Example: Field_buffer p_buffer(grid, field.p); p_buffer.top_send...
 */
struct Field_buffer {
    double* top_send;
    double* top_recv;
    double* bottom_send;
    double* bottom_recv;
    double* left_send;
    double* left_recv;
    double* right_send;
    double* right_recv;
    
    /**
    * @brief Constructor to generate the C-style arrays from fields
    *
    * @param[in] grid Grid object
    * @param[in] m field variable of type Matrix<double>
    */
    Field_buffer(Grid& grid, Matrix<double>& m);

    /**
    * @brief Destructor to deallocate the pointer data members
    */
    ~Field_buffer();

    /**
    * @brief Copy data from receive buffers to respective halo cells in the given Matrix
    *
    * @param[in] grid Grid object
    * @param[in] m field variable of type Matrix<double>
    */
    void buffer_to_halo(Grid& grid, Matrix<double>& m);

};


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
        std::array<int, 4> _neighbours{MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL, MPI_PROC_NULL};
        
        /// Neigbour and neighbour_bool setter (called by constructor)
        void set_neighbours();

    public: 

        Processor() = default; 

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

        /**
        * @brief Communicate all fields (or only pressure) to adjacent processors
        *
        * @param[in] grid
        * @param[in] field
        * @param[in] pressure_only whether we want to communicate only pressure or all fields
        */
        void communicate(Grid& grid, Fields& field, bool pressure_only);

         /** 
         * @brief Reduce min 
         * 
         * @param[in] param;
         * @param[in] my_rank
         */
        double reduce_min(double param,int rank); 

          /** 
         * @brief Reduce min 
         * 
         * @param[in] param;
         * @param[in] my_rank
         */
        double reduce_sum(double param,int rank); 

        /// Getter of x index
        int ip() const;
        /// Getter of y index
        int jp() const;
        /// Getter of processor id
        int proc_id() const;
        //  Getters for proecssor neighbours 
        std::array<bool,4> get_neighbours() {return _neighbours_bool;}
};


