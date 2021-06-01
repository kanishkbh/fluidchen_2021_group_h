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