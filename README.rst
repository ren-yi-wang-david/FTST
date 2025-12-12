==============================================================
FTST
==============================================================

Basic Information
=================
FDM Thermal Simulator Tool

GitHub repository: https://github.com/ren-yi-wang-david/FTST

Problem to Solve
================
steady state heat transfer,thermal distributionetc.

Prospective Users
=================
this tool is developed for CAE,VLSI,CFD simulation user providing verification ,early stage prediction,thermal-aware co optimization.

Requirements
============
- **Python**: 3.10 or newer with `pip`.
- **Python packages** (install via `pip install -r requirements.txt`):
  `numpy`, `matplotlib`, `pandas`, `scipy`, `pytest`.
- **C++ toolchain**: g++ (C++17), OpenMP support, AVX2-capable flags (already configured in `Makefile`).
- **Build utilities**: GNU Make, pybind11 headers (`/usr/include/pybind11` by default), Python dev headers (`python3-config` must succeed).
- **Optional**: Google Test sources for C++ unit tests (path configurable via `GTEST_ROOT` in the Makefile).

System Architecture
===================
- **Discretization**: Finite Difference Method (FDM) converts the steady-state heat equation into a structured 3D grid.
- **Solver core (C++)**: Dense3DSolver implements a red-black Successive Over-Relaxation (SOR) algorithm optimized with AVX2 SIMD gather operations and multithreaded tile processing, enabling high throughput on modern CPUs.
- **Memory tracking**: Custom allocators track peak usage for profiling large grids.
- **Python binding**: `pybind11` exposes the C++ solver (`ftst_dense`), so Python drives configuration, execution, and data extraction.
- **Visualization & comparison**: Python utilities (`FTST.ThermalSolver`, `FTST.compare`) handle CSV exports, heatmaps, and Icepak comparisons for end users.
- **Workflow**: Users set geometry/boundaries → C++ solver iterates (SIMD + threads) → results stream back to Python → plots/CSV/validation artifacts are saved under `output/`.

API Description
===============

The following APIs provide the minimal end-to-end workflow for running the FTST
thermal solver, exporting results, generating cut-plane visualizations, and
optionally comparing against Ansys Icepak reference data.

.. code:: python

    from pathlib import Path
    from FTST import ThermalSolver

    ROOT = Path(__file__).resolve().parent

    def main() -> None:
        """Minimal end-to-end FTST workflow without prints."""

        # 1. Create solver instance
        #    Lx, Ly, Lz define the physical dimensions (meters),
        #    dx defines the uniform grid spacing.
        solver = ThermalSolver(Lx=0.01, Ly=0.01, Lz=0.01, dx=1e-4)

        # 2. Set fixed-temperature boundary conditions
        #    Six values correspond to: x-, x+, y-, y+, z-, z+
        solver.set_boundary([25.0, 25.0, 25.0, 25.0, 25.0, 60.0])

        # 3. Run the iterative SOR solver
        #    max_iter: maximum number of iterations
        #    tolerance: L-infinity residual threshold
        solver.solve(max_iter=5000, tolerance=1e-8)

        # 4. Export the 3D temperature field to CSV
        #    Output format: x_m, y_m, z_m, phi
        out_dir = ROOT / "output"
        out_dir.mkdir(exist_ok=True)
        csv_path = out_dir / "result_from_main.csv"
        solver.temp_csv(csv_path.as_posix())

        # 5. Generate a temperature cut-plane heatmap
        #    axis ∈ {'x','y','z'}, index selects the slice,
        #    filename=None → auto-generate heatmap_<axis>_<index>.png
        cx = solver.nx // 2
        solver.cutplane_plot("x", cx, None)

        # 6. Optional comparison with Icepak reference data
        #    If the Icepak file exists, FTST automatically writes the
        #    latest field and performs error/visual comparison.
        ice_path = ROOT / "tests" / "resources" / "ice.txt"
        if ice_path.exists():
            solver.compare(ice_path.as_posix(), csv_path.as_posix())

    if __name__ == "__main__":
        main()


How to Build FTST with Make
===========================

Compile the FTST C++ core and Python bindings:

.. code:: bash

    make

Clean all build artifacts:

.. code:: bash

    make clean

Install the FTST Python module into your environment:

.. code:: bash

    make install

After building, import the module in Python:

.. code:: python

    from FTST import ThermalSolver


main.py is the example of using FTST
Engineering Infrastructure
==========================
license:BSD 3-Clause License

Testing
==========
- **C++ / Google Test (`make test`)**
  - `Matrix3DTest.*` confirms allocator-backed storage and copy/move semantics.
  - `Dense3DSolverTest.*` exercises constant boundary convergence, zero-boundary stability, and medium-grid stress cases.
  - `GoldenTest.Analytic3D` compares the numeric field against the closed-form `sin(pi x) sin(pi y) sinh(pi z)/sinh(pi L)` solution as a regression check.
- **Python / pytest (`make pytest`)**
  - `tests/test_basic.py` mirrors the core solver use cases from Python to ensure the bindings preserve numerical behavior.
  - `tests/test_boundary.py` verifies each Dirichlet face matches the requested temperature (excluding the inevitable corner overlaps).
  - `tests/test_compare.py` creates synthetic FTST/Icepak CSVs and asserts that `compare_ice_to_ftst` reports zero error when the datasets align and emits the expected artifacts in `output/`.

Compiler: g++, python

Build System: make

Version Control: git

Automatic build system: GNU make

Documentation: README.md

Programming languages: Python (visualization), C++ (computation)

golden:ansys icepak

Schedule(completed)
========    

* Project Schedule (10 weeks from 10/6 to 12/8):

* Week 1 (10/6): Initialize repository and Make/CMake build system; set up project structure; study heat transfer modeling.

* Week 2 (10/13): Complete key notes on heat transfer modeling; draft FDM discretization derivations (documentation); set up baseline Google Test framework.

* Week 3 (10/20): Implement 2D FDM discretization for steady-state heat equation; add Dirichlet boundary condition support; run initial unit tests.

* Week 4 (10/27): Implement LU decomposition solver in C++; validate with small synthetic power maps; perform preliminary validation against Icepak.

* Week 5 (11/3): Expose solver to Python using pybind11; enable NumPy array interop; test Python API with simple cases; document mixed workflow.

* Week 6 (11/10): Re-check results against Icepak golden function; run micro-benchmarks; consolidate error analysis report.

* Week 7 (11/17): Performance analysis and optimize: time, memory, and error.

* Week 8 (11/24): Improve documentation (README, tutorials, API reference); draft v0.1 API freeze.

* Week 9 (12/1): Finalize release v0.1.0; polish code, demos, and pytest coverage.

* Week 10 (12/8): Prepare final report; demo presentation.

References
==========

HeatandMassTransfer7thEdition-Incropera-dewitt
Mittal, S. (2014). A study of successive over-relaxation method parallelisation over modern HPC languages. International journal of high performance computing and networking, 7(4), 292-298.
