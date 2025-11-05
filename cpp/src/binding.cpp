#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "thermal_analysis.hpp"
#include "FTST.cpp"

namespace py = pybind11;

py::array_t<double> solve_py(int Nx, int Ny, double dx, double dy, double k,
                             py::array_t<double> power, py::array_t<double> bc)
{
    auto p_buf = power.unchecked<1>();
    auto b_buf = bc.unchecked<1>();

    std::vector<double> p_vec(p_buf.size()), b_vec(b_buf.size());
    for (ssize_t i = 0; i < p_buf.size(); ++i) p_vec[i] = p_buf(i);
    for (ssize_t i = 0; i < b_buf.size(); ++i) b_vec[i] = b_buf(i);

    auto T = ftst::run_steady_heat(Nx, Ny, dx, dy, k, p_vec, b_vec);

    py::array_t<double> out(T.size());
    auto out_buf = out.mutable_unchecked<1>();
    for (ssize_t i = 0; i < T.size(); ++i) out_buf(i) = T[i];
    return out;
}

PYBIND11_MODULE(FTST, m) {
    m.doc() = "FTST steady-state thermal simulator";
    m.def("solve", &solve_py, "Solve 2D steady heat distribution");
}
