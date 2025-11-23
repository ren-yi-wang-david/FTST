#include <iostream>
#include "mkl_solver.hpp"
#include "parser.hpp"
#include "material.hpp"
using namespace std;
#include <chrono>
#include <fstream>
#include "thermal_analysis.hpp"

int main(int argc, char *argv[])
{



  std::string file_path = "../../input/";
  std::string design_name = "sim";
  std::string output_file_path = "../../output/";
  thermal_analysis TA;
  thermal_solver solver;
  analysis_area AA;
  
  // ---- 1. 設定建模空間起點 ----
  AA.org_point.org_x = 0;
  AA.org_point.org_y = 0;
  AA.org_point.org_z = 0;

  // ---- 2. 設定建模空間尺寸 ----
  AA.dim.dim_x = 2000; // μm
  AA.dim.dim_y = 1500; // μm
  AA.dim.dim_z = 500;  // μm

  // ---- 3. 設定均勻解析度 ----
  AA.global_res.res_x = 50; // X 每格 50 μm
  AA.global_res.res_y = 50; // Y 每格 50 μm
  AA.global_res.res_z = 10; // Z 每格 10 μm

  AA.init_grids();
  TA.Thermal_Analysis();
  std::ofstream fout1(output_file_path + "b_result.txt");
  for (int i = 0; i < TA.thermal_b_matrix.size(); i++)
  {
    if (TA.thermal_b_matrix[i] != 0)
      fout1 << i << " " << TA.thermal_b_matrix[i] << std::endl;
  }
  solver.mkl_solver_csr(
      TA.thermal_MNA_matrix.no_zero_term,
      TA.thermal_MNA_matrix.no_zero_term_column,
      TA.thermal_MNA_matrix.no_zero_term_value,
      TA.thermal_result,
      TA.thermal_b_matrix);
  TA.assign_tempture_to_grids();
  std::ofstream fout(output_file_path + "answer_result.txt");
  double max = -100000;
  double min = 1000000;
  for (int i = 0; i < TA.thermal_result.size(); i++)
  {
    if (TA.thermal_result[i] > max)
      max = TA.thermal_result[i];
    if (TA.thermal_result[i] < min)
      min = TA.thermal_result[i];
    fout << TA.thermal_result[i] << std::endl;
  }
  std::cout << "max tempture: " << max << std::endl;
  std::cout << "min tempture: " << min << std::endl;
  return 0;
}
