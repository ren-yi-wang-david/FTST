#include <iostream>
#include <vector>
#include "mkl_solver.hpp" // thermal_solver definition

int main()
{
    thermal_solver TS;

    // A 3x3 matrix
    std::vector<int> ArowNEntry = {2, 3, 2};

    // row 0: (0,4), (1,-1)
    // row 1: (0,-1), (1,4), (2,-1)
    // row 2: (1,-1), (2,4)

    std::vector<std::vector<int>> ArowCols = {
        {0, 1},
        {0, 1, 2},
        {1, 2}};

    std::vector<std::vector<double>> ArowVals = {
        {4.0, -1.0},
        {-1.0, 4.0, -1.0},
        {-1.0, 4.0}};

    std::vector<double> b = {2.0, 6.0, 2.0};
    std::vector<double> x;

    TS.mkl_solver_csr(ArowNEntry, ArowCols, ArowVals, x, b);

    std::cout << "Solution x: ";
    for (double v : x)
        std::cout << v << " ";
    std::cout << std::endl;

    return 0;
}