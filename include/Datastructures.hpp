#pragma once

#include <vector>
#include <cmath>        // std::abs
#include <algorithm>
#include <iostream>
#include <iomanip>

/**
 * @brief General 2D data structure around std::vector, in column
 * major format.
 *
 */
template <typename T> class Matrix {

  public:
    Matrix<T>() = default;

    /**
     * @brief Move constructor
     * @param[in] Matrix from which data must be moved
     * */
    
    Matrix<T>(Matrix<T>&& rhs) : _imax(rhs._imax), _jmax(rhs._jmax) {
        _container = std::move(rhs._container);
    }

    /**
     * @brief Copy constructor
     * @param[in] Matrix from which data must be copied
     * */
    
    Matrix<T>(const Matrix<T>& rhs) : _imax(rhs._imax), _jmax(rhs._jmax), _container(rhs._container) {
    }

    /**
     * @brief Copy assignment
     * @param[in] Matrix from which data must be copied
     * */
    
    Matrix<T>& operator=(const Matrix<T>& rhs) {
        _imax = rhs._imax;
        _jmax = rhs._jmax;
        _container = rhs._container;


        return *this;
    }

    /**
     * @brief Move assignment
     * @param[in] Matrix from which data must be moved
     * */
    
    Matrix<T>& operator=(Matrix<T>&& rhs) {
        _imax = rhs._imax;
        _jmax = rhs._jmax;
        _container = std::move(rhs._container);


        return *this;
    }

    /**
     * @brief Constructor with initial value
     *
     * @param[in] number of elements in x direction
     * @param[in] number of elements in y direction
     * @param[in] initial value for the elements
     *
     */
    Matrix<T>(int i_max, int j_max, double init_val) : _imax(i_max), _jmax(j_max) {
        _container.resize(i_max * j_max);
        std::fill(_container.begin(), _container.end(), init_val);
    }

    /**
     * @brief Constructor without an initial value.
     *
     * @param[in] number of elements in x direction
     * @param[in] number of elements in y direction
     *
     */
    Matrix<T>(int i_max, int j_max) : _imax(i_max), _jmax(j_max) { _container.resize(i_max * j_max); }

    /**
     * @brief Element access and modify using index
     *
     * @param[in] x index
     * @param[in] y index
     * @param[out] reference to the value
     */
    T &operator()(int i, int j) { return _container.at(_imax * j + i); }

    /**
     * @brief Element access using index
     *
     * @param[in] x index
     * @param[in] y index
     * @param[out] value of the element
     */
    T operator()(int i, int j) const { return _container.at(_imax * j + i); }

    /**
     * @brief Pointer representation of underlying data
     *
     * @param[out] pointer to the beginning of the vector
     */
    const T *data() const { return _container.data(); }

    /**
     * @brief Access of the size of the structure
     *
     * @param[out] size of the data structure
     */
    int size() const { return _container.size(); }

    /// get the given row of the matrix
    std::vector<double> get_row(int row) {
        std::vector<T> row_data(_imax, -1);
        for (int i = 0; i < _imax; ++i) {
            row_data.at(i) = _container.at(i + _imax * row);
        }
        return row_data;
    }

    /// get the given column of the matrix
    std::vector<double> get_col(int col) {
        std::vector<T> col_data(_jmax, -1);
        for (int i = 0; i < _jmax; ++i) {
            col_data.at(i) = _container.at(col + i * _imax);
        }
        return col_data;
    }

    /// set the given column of matrix to given vector
    void set_col(const std::vector<double> &vec, int col) {
        for (int i = 0; i < _jmax; ++i) {
            _container.at(col + i * _imax) = vec.at(i);
        }
    }

    /// set the given row of matrix to given vector
    void set_row(const std::vector<double> &vec, int row) {
        for (int i = 0; i < _imax; ++i) {
            _container.at(i + row * _imax) = vec.at(i);
        }
    }

    /// get the number of elements in x direction
    int imax() const { return _imax; }

    /// get the number of elements in y direction
    int jmax() const { return _jmax; }

    /// Returns the absolute maximum from the fluid domain.
    T max() const {
        return std::fabs(*std::max_element(_container.cbegin(), _container.cend(), 
            [](const T a, const T b) {
            return (std::fabs(a) < std::fabs(b));
            }));
    }

    // Utility for debugging
    void pretty_print(std::ostream& os) const {
        auto f = os.flags();
        os << "-------------------------------" << std::endl;
        os << std::setiosflags(std::ios::fixed);
        os << std::setprecision(8);
        for (int j = jmax() - 1; j >= 0; --j) {
            for (int i = 0; i < imax(); ++i) {
                os << (*this)(i, j) << ' ';
            }
            os << std::endl;
        }
        os << "-------------------------------" << std::endl;
        os.flags(f);
    }

  private:
    /// Number of elements in x direction
    int _imax;
    /// Number of elements in y direction
    int _jmax;

    /// Data container
    std::vector<T> _container;
};
