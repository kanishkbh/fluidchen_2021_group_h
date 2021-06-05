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

void Processor::communicate(Grid& grid, Fields& field, bool pressure_only) {

    MPI_Status status;
    
    size_t buf_size_x = grid.domain().local_size_x - 2; // TODO: inside domain there is a domain_size_x and a size_x. none of them seem to be local size.
    size_t buf_size_y = grid.domain().local_size_y - 2;

    Field_buffer p_buffer(grid, field.p_matrix());
    Field_buffer u_buffer(grid, field.u_matrix());
    Field_buffer v_buffer(grid, field.v_matrix());
    Field_buffer t_buffer(grid, field.t_matrix());

    /* TODO
        1. if (top & bottom neighbours)
            S-R from top and bottom
            else if ( Top & !bottom)
            S-R from top
            else if ( Bottom & ! Top)
            S-R frrom Bottom
        2. Same for Left Right
        3. Assign the received buffer to outermost field cells (halo locations).
        4. Check whether dynamic memory allocation-deallocation logic is correct.
        5. Clarify domain.local_size_x/y with Rahul
    */

    // WIP: Send to top, Receive from Bottom, and then vice versa
    if ()
    MPI_Sendrecv(   (void *)p_buffer.top_send, buf_size_x, MPI_DOUBLE, , 1, 
                    (void *)p_buffer.bottom_recv, buf_size_x, MPI_DOUBLE, ,1, 
                    MPI_COMM_WORLD, &status
                )

}


Field_buffer::Field_buffer(Grid& grid, Matrix<double>& m) {
    
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
        bottom_send = m(i,j);

    // copy top elements from m to array
    j = n_y - 2;
    for ( i = 1; i < n_x - 1; i++)
        top_send = m(i,j);
        
    // copy left elements from m to array
    i = 1;
    for ( j = 1; i < n_y - 1; j++)
        left_send = m(i,j);
    
    // copy right elements from m to array
    i = n_x - 2;
    for ( j = 1; j < n_y - 1; j++)
        right_send = m(i,j);
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