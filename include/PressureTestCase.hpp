#ifndef DEF_PRESSURE_TEST_CASE
#define DEF_PRESSURE_TEST_CASE

#include <vector>
#include "Case.hpp"

class PressureTestCase : public Case {

    public:
    /**
     * @brief Parallel constructor for the Case.
     *
     * Reads input file, creates Fields, Grid, Boundary, Solver and sets
     * Discretization parameters Creates output directory
     *
     * @param[in] Input file name
     */
    PressureTestCase(std::string file_name, int argn, char **args);

    /**
     * @brief Do N iterations of the pressure solver and return a vector containing 
     * residuals at each solver step
     * */
    std::vector<double> pressure_solve(unsigned N);
    
};

#endif