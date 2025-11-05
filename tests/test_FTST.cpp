#include "../cpp/include/thermal_analysis.hpp"
#include <cassert>

int main() {
    ftst::ThermalAnalysis solver(5,5,1.0,1.0,1.0);
    solver.solve_steady(10);
    auto T = solver.get_temperature();
    assert(T.size() == 25);
    return 0;
}
