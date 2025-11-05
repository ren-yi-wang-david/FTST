#include "thermal_analysis.hpp"
#include <vector>

namespace ftst {

std::vector<double> run_steady_heat(
    int Nx, int Ny, double dx, double dy, double k,
    const std::vector<double>& power,
    const std::vector<double>& bc)
{
    ThermalAnalysis solver(Nx, Ny, dx, dy, k);
    solver.set_power(power);
    solver.set_boundary(bc);
    solver.solve_steady();
    return solver.get_temperature();
}

} // namespace ftst
