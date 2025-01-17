#include "Communication.hpp"

int Communication::_bottom_neighbor_rank = MPI_PROC_NULL;
int Communication::_top_neighbor_rank = MPI_PROC_NULL;
int Communication::_left_neighbor_rank = MPI_PROC_NULL;
int Communication::_right_neighbor_rank = MPI_PROC_NULL;
int Communication::_world_size = 1;
int Communication::_rank = 0;


void Communication::communicate_right(Matrix<double>& x, int tag) {
    if (_left_neighbor_rank == MPI_PROC_NULL && _right_neighbor_rank == MPI_PROC_NULL)
        return; //No horizontal communication

    std::vector<double> out = x.get_col(x.imax()-2);
    std::vector<double> in(x.jmax());

    if (_left_neighbor_rank == MPI_PROC_NULL) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_right_neighbor_rank == MPI_PROC_NULL)
        MPI_Recv(in.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag,
                     in.data(), x.jmax(), MPI_DOUBLE,
                     _left_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the left)
    for (int j = 0; j < x.jmax(); ++j)
        x(0, j) = in[j];
}

void Communication::communicate_left(Matrix<double>& x, int tag) {
    if (_left_neighbor_rank == MPI_PROC_NULL && _right_neighbor_rank == MPI_PROC_NULL)
        return; //No horizontal communication

    std::vector<double> out = x.get_col(1);
    std::vector<double> in(x.jmax());

    if (_right_neighbor_rank == MPI_PROC_NULL) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_left_neighbor_rank == MPI_PROC_NULL)
        MPI_Recv(in.data(), x.jmax(), MPI_DOUBLE, _right_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.jmax(), MPI_DOUBLE, _left_neighbor_rank, tag,
                     in.data(), x.jmax(), MPI_DOUBLE,
                     _right_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the right)
    auto i = x.imax()-1;
    for (int j = 0; j < x.jmax(); ++j)
        x(i, j) = in[j];
}

void Communication::communicate_top(Matrix<double>& x, int tag) {
    if (_top_neighbor_rank == MPI_PROC_NULL && _bottom_neighbor_rank == MPI_PROC_NULL)
        return; //No vertical communication

    std::vector<double> out = x.get_row(x.jmax()-2);
    std::vector<double> in(x.imax());

    if (_bottom_neighbor_rank == MPI_PROC_NULL) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_top_neighbor_rank == MPI_PROC_NULL)
        MPI_Recv(in.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag,
                     in.data(), x.imax(), MPI_DOUBLE,
                     _bottom_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the bottom)
    for (int i = 0; i < x.imax(); ++i)
        x(i, 0) = in[i];
}

void Communication::communicate_bottom(Matrix<double>& x, int tag) {

    if (_top_neighbor_rank == MPI_PROC_NULL && _bottom_neighbor_rank == MPI_PROC_NULL)
        return; //No vertical communication

    std::vector<double> out = x.get_row(1);
    std::vector<double> in(x.imax());

    if (_top_neighbor_rank == MPI_PROC_NULL) {
        //Nothing to read, just send then return
        MPI_Send(out.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag, MPI_COMM_WORLD);
        return;
    }

    // Read and, if necessary, send
    if (_bottom_neighbor_rank == MPI_PROC_NULL)
        MPI_Recv(in.data(), x.imax(), MPI_DOUBLE, _top_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    else {
        MPI_Sendrecv(out.data(), x.imax(), MPI_DOUBLE, _bottom_neighbor_rank, tag,
                     in.data(), x.imax(), MPI_DOUBLE,
                     _top_neighbor_rank, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    //Update x (from the top)
    auto j = x.jmax() - 1;
    for (int i = 0; i < x.imax(); ++i)
        x(i, j) = in[i];
}

void Communication::communicate_all(Matrix<double>& x, int tag) {
    communicate_right(x, tag);
    communicate_top(x, tag);
    communicate_left(x, tag);
    communicate_bottom(x, tag);
}

void Communication::communicate_sum_double(double* src, double* target) {
    MPI_Allreduce(src, target, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


void Communication::communicate_sum_int(int* src, int* target) {
    MPI_Allreduce(src, target, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
}

double Communication::reduce_min(double x) {
    double res;
    MPI_Allreduce(&x, &res, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return res;

}