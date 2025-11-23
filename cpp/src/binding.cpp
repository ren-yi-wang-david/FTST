#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <fstream>
#include <sstream>
#include <nlohmann/json.hpp>
#include "FTST.hpp"

namespace py = pybind11;
using json = nlohmann::json;

namespace ftst {

// ==============================================
// 全域暫存資料
// ==============================================
Geometry geometry;
Material material;
Boundary boundary;
PowerMap power_map;

// ==============================================
// JSON & CSV 載入
// ==============================================
void load_geometry(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("❌ geometry.json not found: " + path);
    json data; f >> data;

    geometry.name = data.value("name", "default_geometry");
    geometry.Nx = data.value("Nx", 10);
    geometry.Ny = data.value("Ny", 10);
    geometry.Nz = data.value("Nz", 5);
    geometry.dx = data.value("dx", 1.0);
    geometry.dy = data.value("dy", 1.0);
    geometry.dz = data.value("dz", 1.0);

    std::cout << "✅ Loaded geometry: " << geometry.name
              << " (" << geometry.Nx << "x" << geometry.Ny << "x" << geometry.Nz << ")\n";
}

void load_material(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("❌ materials.json not found: " + path);
    json data; f >> data;

    material.density = data["density"];
    material.specific_heat = data["specific_heat"];
    material.kx = data["conductivity"][0];
    material.ky = data["conductivity"][1];
    material.kz = data["conductivity"][2];

    std::cout << "✅ Loaded material: k=(" << material.kx << ", " << material.ky << ", " << material.kz << ")\n";
}

void load_boundary(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("❌ boundary.json not found: " + path);
    json data; f >> data;

    boundary.htc_top = data["htc"]["top"];
    boundary.htc_bottom = data["htc"]["bottom"];
    boundary.htc_side = data["htc"]["side"];
    boundary.ambient_T = data["ambient_temperature"];

    std::cout << "✅ Loaded boundary: HTC_top=" << boundary.htc_top
              << " Ambient=" << boundary.ambient_T << "°C\n";
}

void load_power(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) throw std::runtime_error("❌ power_map.csv not found: " + path);
    power_map.values.clear();

    std::string line;
    while (std::getline(f, line)) {
        std::stringstream ss(line);
        std::string val;
        while (std::getline(ss, val, ',')) {
            power_map.values.push_back(std::stod(val));
        }
    }

    std::cout << "✅ Loaded power map, size=" << power_map.values.size() << "\n";
}

// ==============================================
// Mock solver (生成虛擬溫度場)
// ==============================================
std::vector<double> solve() {
    int Nx = geometry.Nx, Ny = geometry.Ny, Nz = geometry.Nz;
    std::vector<double> T(Nx * Ny * Nz);

    double base = boundary.ambient_T;
    double scale = material.kx / 10.0;

    for (int k = 0; k < Nz; ++k) {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                int idx = k * Nx * Ny + j * Nx + i;
                T[idx] = base + scale * std::sin((double)i / Nx * 3.14) *
                                  std::sin((double)j / Ny * 3.14) *
                                  std::exp(-k * 0.2);
            }
        }
    }
    std::cout << "🔥 Mock solver finished.\n";
    return T;
}

// ==============================================
// 結果輸出
// ==============================================
void save_result(const std::string& path) {
    std::ofstream f(path);
    if (!f.is_open()) throw std::runtime_error("❌ Cannot write file: " + path);
    f << "Mock result file.\n";
    f.close();
    std::cout << "💾 Saved result to " << path << "\n";
}

}  // namespace ftst

// ==============================================
// pybind11 module
// ==============================================
PYBIND11_MODULE(FTST, m) {
    m.doc() = "FTST thermal simulator mock API";

    m.def("load_geometry", &ftst::load_geometry);
    m.def("load_material", &ftst::load_material);
    m.def("load_boundary", &ftst::load_boundary);
    m.def("load_power", &ftst::load_power);
    m.def("solve", &ftst::solve);
    m.def("save_result", &ftst::save_result);
}
