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

void Processor::communicate(Grid& grid, Matrix<double>& p) {

    MPI_Status status;
    
    size_t buf_size_x = grid.domain().size_x;  
    size_t buf_size_y = grid.domain().size_y;

    // copying data from fields matrices to send buffers
    Field_buffer p_buffer(grid, p);


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
        // 1. Send to Top, Receive from Bottom
        // TODO : Experiment without void* 
        MPI_Sendrecv(   &p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        &p_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 1, 
                        MPI_COMM_WORLD, &status
                    );
        // 2. Send to Bottom , Receive from Top
        MPI_Sendrecv(   &p_buffer.bottom_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::BOTTOM), 2, 
                        &p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 2, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to TOP, Receive from TOP (Bottom row of processors)
    else if (has_neighbour(border_position::TOP) && !has_neighbour(border_position::BOTTOM)) {
        // Send to Top, Receive from Top
        MPI_Sendrecv(   &p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        &p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to BOTTOM, Receive from BOTTOM (Top row of processors)
    else if (has_neighbour(border_position::TOP) && !has_neighbour(border_position::BOTTOM)) {
        // Send to Top, Receive from Top
        MPI_Sendrecv(   &p_buffer.top_send, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        &p_buffer.top_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::TOP), 1, 
                        MPI_COMM_WORLD, &status
                    );
    }

   /// Set 2: HORIZONTAL Communication 
    if (has_neighbour(border_position::LEFT) && has_neighbour(border_position::RIGHT)) {

        // 1. Send to LEFT, Receive from RIGHT
        MPI_Sendrecv(   &p_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        &p_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        MPI_COMM_WORLD, &status
                    );
        // 2. Send to RIGHT , Receive from LEFT
        MPI_Sendrecv(   &p_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 10, 
                        &p_buffer.left_recv, buf_size_x, MPI_DOUBLE, neighbour(border_position::LEFT), 10, 
                        MPI_COMM_WORLD, &status
    }
    // Send to LEFT, Receive from LEFT (for Right-most row of processors)
    else if (has_neighbour(border_position::LEFT) && !has_neighbour(border_position::RIGHT)) {
         // Send to Left, Receive from Left
        MPI_Sendrecv(   &p_buffer.left_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        &p_buffer.left_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::LEFT), 9, 
                        MPI_COMM_WORLD, &status
                    );
    }
    // Send to RIGHT, Receive from RIGHT (for Left-most row of processors)
    else if (has_neighbour(border_position::RIGHT) && !has_neighbour(border_position::LEFT)) {
        // Send to Right, Receive from Right
        MPI_Sendrecv(   &p_buffer.right_send, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        &p_buffer.right_recv, buf_size_y, MPI_DOUBLE, neighbour(border_position::RIGHT), 9, 
                        MPI_COMM_WORLD, &status
                    );
    }


    // Copying data from recv buffers to fields natrices
    p_buffer.buffer_to_halo(grid, field.p_matrix());
}


Field_buffer::Field_buffer(Grid& grid, Matrix<double>& m) {
/// TODO: Deal with border cases
    
    // no of elements to send/recv
    auto n_x = grid.domain().size_x; // TODO: inside domain there is a domain_size_x and a size_x. none of them seem to be local size.
    auto n_y = grid.domain().size_y;

    // Allocate C-style arrays
    top_send    = new double[n_x];
    top_recv    = new double[n_x];
    bottom_recv = new double[n_x];
    bottom_send = new double[n_x]; 
    left_send   = new double[n_y];
    left_recv   = new double[n_y];
    right_send  = new double[n_y];
    right_recv  = new double[n_y];

    int i=0, j=0;

    // copy bottom elements from m to array
    j = 1;
    for ( i = 1; i < n_x + 1; i++)
        bottom_send[i-1] = m(i,j);

    // copy top elements from m to array
    j = n_y;
    for ( i = 1; i < n_x + 1; i++)
        top_send[i-1] = m(i,j);
        
    // copy left elements from m to array
    i = 1;
    for ( j = 1; i < n_y + 1; j++)
        left_send[j-1] = m(i,j);
    
    // copy right elements from m to array
    i = n_x;
    for ( j = 1; j < n_y + 1; j++)
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

    auto n_x = grid.domain().size_x; /// TODO: inside domain there is a domain_size_x and a size_x. Semantically none of them seem to be local size.
    auto n_y = grid.domain().size_y;

    int i=0, j=0;

    // copy bottom elements from array to m
    j = 0;
    for ( i = 1; i < n_x + 1 ; i++)
        m(i,j) = bottom_recv[i-1];

    // copy top elements from array to m
    j = n_y+1 ;
    for ( i = 1; i < n_x + 1; i++)
        m(i,j) = top_recv[i-1];
        
    // copy left elements from array to m
    i = 0;
    for ( j = 1; i < n_y + 1; j++)
        m(i,j) = left_recv[j-1];
    
    // copy right elements from array to m
    i = n_x+1;
    for ( j = 1; j < n_y + 1; j++)
        m(i,j) = right_recv[j-1];
}

// We only require double for our application,hence this function is not templatized. 
double Processor::reduce_min(double param){
    double param_comm = param;
    MPI_Allreduce(&param,&param_comm,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD); 
    return param_comm;
}

// Likewise we only need a bool version 
double Processor::reduce_sum(double rloc){
    double res=0; 
    MPI_Allreduce(&rloc,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    return res; 
} 