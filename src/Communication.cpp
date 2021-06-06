#include "Communication.hpp"

#include <iostream>

Processor::Processor(int iproc, int jproc, int rank): _iproc(iproc), _jproc(jproc), _id(rank) {
    int _jp = static_cast<int> ( rank / jproc );
    int _ip = static_cast<int> ( rank % iproc );
}

const int Processor::neighbour(border_position position) const {
    return _neighbours.at(static_cast<int>(position));
}

void Processor::set_neighbours() {
    // set TOP neighbour of this processor
    if (_jp != (_jproc-1) ) {
        int top_proc_rank = _id + _iproc;
        _neighbours.at(static_cast<int>(border_position::TOP)) = top_proc_rank;
        _neighbours_bool.at(static_cast<int>(border_position::TOP)) = true;
    }
    // set BOTTOM neighbour of this processor
    if (_jp != 0 ) {
        int bottom_proc_rank = _id - _iproc;
        _neighbours.at(static_cast<int>(border_position::BOTTOM)) = bottom_proc_rank;
        _neighbours_bool.at(static_cast<int>(border_position::BOTTOM)) = true;
    }
    // set RIGHT neighbour of this processor
    if (_ip != (_iproc-1) ) {
        int right_proc_rank = _id + 1;
        _neighbours.at(static_cast<int>(border_position::RIGHT)) = right_proc_rank;
        _neighbours_bool.at(static_cast<int>(border_position::RIGHT)) = true;
    }
    // set LEFT neighbour of this processor
    if (_ip != 0 ) {
        int left_proc_rank = _id - 1;
        _neighbours.at(static_cast<int>(border_position::LEFT)) = left_proc_rank;
        _neighbours_bool.at(static_cast<int>(border_position::LEFT)) = true;
    }
}

bool Processor::has_neighbour(border_position position) const {
    return _neighbours_bool.at(static_cast<int>(position));
}

int Processor::ip() const {
    return _ip;
}

int Processor::jp() const {
    return _jp;
}

int Processor::proc_id() const {
    return _id;
}

void Processor::communicate(Grid& grid, Fields& field, bool pressure_only = false) {

    MPI_Status status;
    
    size_t buf_size_x = grid.domain().local_size_x - 2; // TODO: inside domain there is a domain_size_x and a size_x. none of them seem to be local size.
    size_t buf_size_y = grid.domain().local_size_y - 2;

    // copying data from fields matrices to send buffers
    Field_buffer p_buffer(grid, field.p_matrix());
    Field_buffer u_buffer(grid, field.u_matrix());
    Field_buffer v_buffer(grid, field.v_matrix());
    Field_buffer t_buffer(grid, field.t_matrix());


    /* ALGORITHM (S=send, R=recv)
        1.  VERTICAL COMMUNICATION
            if (top & bottom neighbours)
                1.1 S to top, R from bottom
                1.2 S to bottom, R from top 
            else if ( Top & !bottom)
                S-R from top
            else if ( Bottom & ! Top)
                S-R from Bottom
        2. Same for Left Right - HORIZONTAL COMMUNICATION
    */

   /// Set 1: VERTICAL Communication 
    if (has_neighbour(border_position::TOP) && has_neighbour(border_position::BOTTOM)) {
        /// PRESSURE
        // 1. Send to Top, Receive from Bottom
        MPI_Sendrecv(   (void *)p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        (void *)p_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 1, 
                        MPI_COMM_WORLD, &status
                    );
        // 2. Send to Bottom , Receive from Top
        MPI_Sendrecv(   (void *)p_buffer.bottom_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 2, 
                        (void *)p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 2, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 3, 
                        (void *)u_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 3, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)u_buffer.bottom_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 4, 
                        (void *)u_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 4, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 5, 
                        (void *)v_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 5, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)v_buffer.bottom_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 6, 
                        (void *)v_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 6, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 7, 
                        (void *)t_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 7, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)t_buffer.bottom_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM),8, 
                        (void *)t_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 8, 
                        MPI_COMM_WORLD, &status
                    ); 
    }
    // Send to TOP, Receive from TOP (Bottom row of processors)
    else if (has_neighbour(border_position::TOP) && !has_neighbour(border_position::BOTTOM)) {
        /// PRESSURE // Send to Top, Receive from Top
        MPI_Sendrecv(   (void *)p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        (void *)p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 3, 
                        (void *)u_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 3, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 5, 
                        (void *)v_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 5, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 7, 
                        (void *)t_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 7, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to BOTTOM, Receive from BOTTOM (Top row of processors)
    else if (has_neighbour(border_position::TOP) && !has_neighbour(border_position::BOTTOM)) {
        /// PRESSURE // Send to Top, Receive from Top
        MPI_Sendrecv(   (void *)p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        (void *)p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 3, 
                        (void *)u_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 3, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 5, 
                        (void *)v_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 5, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 7, 
                        (void *)t_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 7, 
                        MPI_COMM_WORLD, &status
                    );
    }

   /// Set 2: HORIZONTAL Communication 
    if (has_neighbour(border_position::LEFT) && has_neighbour(border_position::RIGHT)) {
        /// PRESSURE
        // 1. Send to LEFT, Receive from RIGHT
        MPI_Sendrecv(   (void *)p_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        (void *)p_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        MPI_COMM_WORLD, &status
                    );
        // 2. Send to RIGHT , Receive from LEFT
        MPI_Sendrecv(   (void *)p_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 10, 
                        (void *)p_buffer.left_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::LEFT), 10, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 11, 
                        (void *)u_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 11, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)u_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 12, 
                        (void *)u_buffer.left_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::LEFT), 12, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 13, 
                        (void *)v_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)v_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 14, 
                        (void *)v_buffer.left_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::LEFT), 14, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 13, 
                        (void *)t_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        MPI_COMM_WORLD, &status
                    );
        MPI_Sendrecv(   (void *)t_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 14, 
                        (void *)t_buffer.left_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::LEFT), 14, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to LEFT, Receive from LEFT (for Right-most row of processors)
    else if (has_neighbour(border_position::LEFT) && !has_neighbour(border_position::RIGHT)) {
        /// PRESSURE // Send to Left, Receive from Left
        MPI_Sendrecv(   (void *)p_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        (void *)p_buffer.left_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 11, 
                        (void *)u_buffer.left_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 11, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 13, 
                        (void *)v_buffer.left_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 13, 
                        (void *)t_buffer.left_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 13, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to RIGHT, Receive from RIGHT (for Left-most row of processors)
    else if (has_neighbour(border_position::RIGHT) && !has_neighbour(border_position::LEFT)) {
        /// PRESSURE // Send to Right, Receive from Right
        MPI_Sendrecv(   (void *)p_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        (void *)p_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        MPI_COMM_WORLD, &status
                    );
        /// X-VELOCITY
        MPI_Sendrecv(   (void *)u_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 11, 
                        (void *)u_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 11, 
                        MPI_COMM_WORLD, &status
                    );
        /// Y-VELOCITY
        MPI_Sendrecv(   (void *)v_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        (void *)v_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        MPI_COMM_WORLD, &status
                    );
        /// TEMPERATURE
        MPI_Sendrecv(   (void *)t_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        (void *)t_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 13, 
                        MPI_COMM_WORLD, &status
                    );
      
    }


    // Copying data from recv buffers to fields natrices
    p_buffer.buffer_to_halo(grid, field.p_matrix());
    u_buffer.buffer_to_halo(grid, field.u_matrix());
    v_buffer.buffer_to_halo(grid, field.v_matrix());
    t_buffer.buffer_to_halo(grid, field.t_matrix());
    
    /// TODO: Check whether dynamic memory allocation-deallocation logic is correct.
    /// TODO: Clarify domain.local_size_x/y with Rahul
    
}


Field_buffer::Field_buffer(Grid& grid, Matrix<double>& m) {
/// TODO: Deal with border cases
    
    // no of elements to send/recv
    auto n_x = grid.domain().local_size_x; // TODO: inside domain there is a domain_size_x and a size_x. none of them seem to be local size.
    auto n_y = grid.domain().local_size_y;

    // Allocate C-style arrays
    top_send    = new double[n_x - 2];
    top_recv    = new double[n_x - 2];
    bottom_recv = new double[n_x - 2];
    left_send   = new double[n_x - 2];
    left_recv   = new double[n_y - 2];
    right_send  = new double[n_y - 2];
    right_recv  = new double[n_y - 2];

    int i=0, j=0;

    // copy bottom elements from m to array
    j = 1;
    for ( i = 1; i < n_x - 1; i++)
        bottom_send[i-1] = m(i,j);

    // copy top elements from m to array
    j = n_y - 1;
    for ( i = 1; i < n_x - 1; i++)
        top_send[i-1] = m(i,j);
        
    // copy left elements from m to array
    i = 1;
    for ( j = 1; i < n_y - 1; j++)
        left_send[j-1] = m(i,j);
    
    // copy right elements from m to array
    i = n_x - 1;
    for ( j = 1; j < n_y - 1; j++)
        right_send[j-1] = m(i,j);
}


Field_buffer::~Field_buffer() {
    delete[] top_send;
    delete[] top_recv;
    delete[] bottom_send;
    delete[] bottom_recv;
    delete[] left_send;
    delete[] left_recv;
    delete[] right_send;
    delete[] right_recv;
}


void Field_buffer::buffer_to_halo(Grid& grid, Matrix<double>& m) {
/// TODO: Deal with border cases

    auto n_x = grid.domain().local_size_x; /// TODO: inside domain there is a domain_size_x and a size_x. Semantically none of them seem to be local size.
    auto n_y = grid.domain().local_size_y;

    int i=0, j=0;

    // copy bottom elements from array to m
    j = 1;
    for ( i = 1; i < n_x - 1; i++)
        m(i,j) = bottom_recv[i-1];

    // copy top elements from array to m
    j = n_y - 1;
    for ( i = 1; i < n_x - 1; i++)
        m(i,j) = top_recv[i-1];
        
    // copy left elements from array to m
    i = 1;
    for ( j = 1; i < n_y - 1; j++)
        m(i,j) = left_send[j-1];
    
    // copy right elements from array to m
    i = n_x - 1;
    for ( j = 1; j < n_y - 1; j++)
        m(i,j) = right_send[j-1];
}
