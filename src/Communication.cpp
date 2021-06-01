#include "Communication.hpp"

#include <iostream>

Processor::Processor(int i, int j, int id): -ip(i), _jp(j), _id(id) {}

const Processor* Processor::neighbour(border_position position) {
    return _neighbours.at(static_cast<int>(position));
}

void Processor::set_neighbour(Processor *processor, border_position position) {
    _neighbours.at(static_cast<int>(position)) = processor;
}

bool Processor::has_neighbour(border_position position) {
    return has_neigbour.at(static_cast<int>(position));
}

int Processor::ip() const {
    return _ip;
}

int Processor::jp() const {
    return _jp;
}

int Processor::id() const {
    return _id;
}