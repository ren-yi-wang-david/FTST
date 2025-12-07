#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <array>
#include <algorithm>
#include <cctype>
#include <cmath>
#include <sstream>

#include "solver.hpp"

namespace py = pybind11;

namespace
{

struct SliceInfo
{
    py::array data;
    double x_min;
    double x_max;
    double y_min;
    double y_max;
    std::string xlabel;
    std::string ylabel;
    std::string title;
};

void set_rhs_from_array(
    Dense3DSolver &solver,
    const py::array_t<double, py::array::c_style | py::array::forcecast> &array)
{
    auto buf = array.request();
    if (buf.ndim != 3)
    {
        throw py::value_error("RHS array must be 3-dimensional");
    }

    if (static_cast<size_t>(buf.shape[0]) != solver.nx() ||
        static_cast<size_t>(buf.shape[1]) != solver.ny() ||
        static_cast<size_t>(buf.shape[2]) != solver.nz())
    {
        throw py::value_error("RHS array shape does not match solver grid size");
    }

    Matrix3D<double> rhs(solver.nx(), solver.ny(), solver.nz());
    auto *data = static_cast<const double *>(buf.ptr);

    size_t index = 0;
    for (size_t i = 0; i < solver.nx(); ++i)
        for (size_t j = 0; j < solver.ny(); ++j)
            for (size_t k = 0; k < solver.nz(); ++k)
                rhs(i, j, k) = data[index++];

    solver.setRHS(rhs);
}

py::array_t<double> solution_to_array(const Dense3DSolver &solver)
{
    const auto &phi = solver.solution();
    py::array_t<double> result(
        {solver.nx(), solver.ny(), solver.nz()});

    auto buf = result.mutable_unchecked<3>();
    for (size_t i = 0; i < solver.nx(); ++i)
        for (size_t j = 0; j < solver.ny(); ++j)
            for (size_t k = 0; k < solver.nz(); ++k)
                buf(i, j, k) = phi(i, j, k);

    return result;
}

SliceInfo make_slice(Dense3DSolver &solver,
                     const std::string &axis,
                     size_t index)
{
    const auto &phi = solver.solution();
    const double dx = solver.dx();
    const double lx = (solver.nx() > 1) ? (solver.nx() - 1) * dx : 0.0;
    const double ly = (solver.ny() > 1) ? (solver.ny() - 1) * dx : 0.0;
    const double lz = (solver.nz() > 1) ? (solver.nz() - 1) * dx : 0.0;

    if (axis == "z")
    {
        if (index >= solver.nz())
            throw py::value_error("slice index out of range for z-axis");

        py::array_t<double> slice({solver.ny(), solver.nx()});
        auto buf = slice.mutable_unchecked<2>();
        for (size_t j = 0; j < solver.ny(); ++j)
            for (size_t i = 0; i < solver.nx(); ++i)
                buf(j, i) = phi(i, j, index);

        std::ostringstream title;
        title << "Temperature Cut-plane at z=" << index * dx << " m";

        return {std::move(slice),
                0.0,
                lx,
                0.0,
                ly,
                "x (m)",
                "y (m)",
                title.str()};
    }

    if (axis == "y")
    {
        if (index >= solver.ny())
            throw py::value_error("slice index out of range for y-axis");

        py::array_t<double> slice({solver.nz(), solver.nx()});
        auto buf = slice.mutable_unchecked<2>();
        for (size_t k = 0; k < solver.nz(); ++k)
            for (size_t i = 0; i < solver.nx(); ++i)
                buf(k, i) = phi(i, index, k);

        std::ostringstream title;
        title << "Temperature Cut-plane at y=" << index * dx << " m";

        return {std::move(slice),
                0.0,
                lx,
                0.0,
                lz,
                "x (m)",
                "z (m)",
                title.str()};
    }

    if (axis == "x")
    {
        if (index >= solver.nx())
            throw py::value_error("slice index out of range for x-axis");

        py::array_t<double> slice({solver.ny(), solver.nz()});
        auto buf = slice.mutable_unchecked<2>();
        for (size_t j = 0; j < solver.ny(); ++j)
            for (size_t k = 0; k < solver.nz(); ++k)
                buf(j, k) = phi(index, j, k);

        std::ostringstream title;
        title << "Temperature Cut-plane at x=" << index * dx << " m";

        return {std::move(slice),
                0.0,
                lz,
                0.0,
                ly,
                "z (m)",
                "y (m)",
                title.str()};
    }

    throw py::value_error("axis must be 'x', 'y', or 'z'");
}

void plot_temperature_slice(Dense3DSolver &solver,
                            const std::string &axis,
                            py::object position_obj,
                            py::object index_obj,
                            const std::string &cmap)
{
    std::string normalized = axis;
    std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                   [](unsigned char c)
                   { return static_cast<char>(std::tolower(c)); });

    size_t axis_len = 0;
    if (normalized == "z")
        axis_len = solver.nz();
    else if (normalized == "y")
        axis_len = solver.ny();
    else if (normalized == "x")
        axis_len = solver.nx();
    else
        throw py::value_error("axis must be 'x', 'y', or 'z'");

    if (axis_len == 0)
        throw py::value_error("solver contains no points along the requested axis");

    const double dx = solver.dx();
    size_t idx = axis_len / 2;

    if (!position_obj.is_none())
    {
        if (!index_obj.is_none())
            throw py::value_error("Specify either position or index, not both");

        const double pos = position_obj.cast<double>();
        if (pos < 0.0)
            throw py::value_error("position must be non-negative");
        idx = static_cast<size_t>(std::llround(pos / dx));
        if (idx >= axis_len)
            throw py::value_error("position is outside the solver domain");
    }
    else if (!index_obj.is_none())
    {
        idx = index_obj.cast<size_t>();
        if (idx >= axis_len)
            throw py::value_error("slice index out of range");
    }

    auto info = make_slice(solver, normalized, idx);

    auto plt = py::module_::import("matplotlib.pyplot");
    plt.attr("figure")(py::arg("figsize") = py::make_tuple(7.0, 6.0));

    auto im = plt.attr("imshow")(info.data,
                                   py::arg("origin") = "lower",
                                   py::arg("extent") = py::make_tuple(info.x_min, info.x_max, info.y_min, info.y_max),
                                   py::arg("aspect") = "auto",
                                   py::arg("cmap") = cmap);
    plt.attr("colorbar")(im, py::arg("label") = "Temperature (°C)");
    plt.attr("xlabel")(info.xlabel);
    plt.attr("ylabel")(info.ylabel);
    plt.attr("title")(info.title);
    plt.attr("tight_layout")();
    plt.attr("show")();
}

} // namespace

PYBIND11_MODULE(ftst_dense, m)
{
    m.doc() = "Pybind11 bindings for the FTST Dense3DSolver";

    py::class_<Dense3DSolver>(m, "Dense3DSolver")
        .def(py::init<size_t, size_t, size_t, double>(),
             py::arg("nx"), py::arg("ny"), py::arg("nz"), py::arg("dx"),
             "Create a solver from explicit grid counts and spacing (meters)")
        .def(py::init<double, double, double, double>(),
             py::arg("lx"), py::arg("ly"), py::arg("lz"), py::arg("dx"),
             "Create a solver from physical lengths (meters) and spacing")
        .def_property_readonly("nx", &Dense3DSolver::nx)
        .def_property_readonly("ny", &Dense3DSolver::ny)
        .def_property_readonly("nz", &Dense3DSolver::nz)
        .def_property_readonly("dx", &Dense3DSolver::dx)
        .def("set_rhs", &set_rhs_from_array, py::arg("array"),
             "Provide the RHS field as a NumPy array with shape (nx, ny, nz)")
        .def("set_boundary",
             [](Dense3DSolver &self, const std::array<double, 6> &bc)
             {
                 self.setBoundary(bc);
             },
             py::arg("values"),
             "Set Dirichlet boundary conditions in the order [x-, x+, y-, y+, z-, z+]")
        .def("solve", &Dense3DSolver::solve,
             py::arg("max_iter"), py::arg("tolerance"))
        .def("solution_array", &solution_to_array,
             "Return the computed solution as a NumPy array");

    m.def("plot_temperature_slice", &plot_temperature_slice,
          py::arg("solver"),
          py::arg("axis") = "z",
          py::arg("position") = py::none(),
          py::arg("index") = py::none(),
          py::arg("cmap") = "inferno",
          "Plot a temperature slice using matplotlib.\n"
          "axis: 'x', 'y', or 'z' (default 'z').\n"
          "position: optional physical location in meters; overrides index.\n"
          "index: optional slice index; defaults to the mid-plane.");
}
