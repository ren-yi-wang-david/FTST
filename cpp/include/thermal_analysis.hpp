#pragma once
#include <iostream>
#include <vector>
#include <string>
#include<unordered_map>
#include<algorithm>
#include "parser.hpp"
class sparse_matrix
{
public:
  sparse_matrix() {};
  std::vector<int> no_zero_term;
  std::vector<std::vector<int>> no_zero_term_column;
  std::vector<std::vector<SADB>> no_zero_term_value;
};
class thermal_analysis
{
public:
  thermal_analysis()
  {
	total_nodes_number = -1;
	AnalArea = nullptr;
  };
  int total_nodes_number;
  analysis_area* AnalArea;
  sparse_matrix thermal_MNA_matrix;
  std::vector<SADB> thermal_b_matrix;
  std::vector<SADB> thermal_result;
  void stamp_thermal_MNA();
  void init_thermal_saprse_m();
  void stamp_resitance_to_sparsem
  (
	grid* Grid,
	grid* neighbor,
	SADB coefficient,
	int nonzero_index
  );
  void stamp_power_to_b_matrix();//have bugs
  void stamp_power_to_b_matrix_ver2();
  void stamp_power_to_b_matrix_ver3();//tiling+parallel+buffered write
  void Thermal_Analysis();
  void assign_tempture_to_grids();
  void stamp_power_to_b_matrix_ver4();//hash+tiling+parallel+buffered write
  void stamp_power_to_b_matrix_ver5();//tiling+parallel+buffered write
  void stamp_power_to_b_matrix_ver6();//tiling+parallel+buffered write
  
};


SADB caculate_thermal_resitance(
    SADB mat_value,
    int area_case, // 1: yz面積, 2: xz面積, 3: xy面積
    grid *Grid);