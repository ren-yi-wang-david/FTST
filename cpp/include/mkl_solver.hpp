#pragma once
#include <vector>

class thermal_solver
{
public:
    thermal_solver() {}

    void mkl_solver_csr(
        std::vector<int> &ArowNEntry,
        std::vector<std::vector<int>> &ArowCols,
        std::vector<std::vector<double>> &ArowVals,
        std::vector<double> &x,
        std::vector<double> &b);
};
