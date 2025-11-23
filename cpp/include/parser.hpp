#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include "material.hpp"
#include <unordered_set>
#include <iomanip>
typedef unsigned int SAINT;
typedef double SADB;
extern double unit; // 除上這個會變公尺
class point
{
public:
  point()
  {
    org_x = 0;
    org_y = 0;
    org_z = 0;
  }
  SAINT org_x;
  SAINT org_y;
  SAINT org_z;
};
class dimension
{
public:
  dimension()
  {
    dim_x = 0;
    dim_y = 0;
    dim_z = 0;
  }
  SAINT dim_x;
  SAINT dim_y;
  SAINT dim_z;
};
class resoultion
{
public:
  resoultion()
  {
    res_x = 0;
    res_y = 0;
    res_z = 0;
  }
  SAINT res_x;
  SAINT res_y;
  SAINT res_z;
};

class thermal_constraint
{
public:
  thermal_constraint()
  {
    constant_tempture = -1;
    heat_flux = -1;
    convection = -1;
  }
  double constant_tempture;
  double heat_flux;
  double convection;
};
class grid
{
public:
  grid()
  {
    index = -1;
    gz_top = nullptr;
    gz_bottom = nullptr;
    gx_left = nullptr;
    gx_right = nullptr;
    gy_up = nullptr;
    gy_down = nullptr;
    boundary_type = -1;
    boundary_constraint = 0;

    tempture = 0.0;
    thermal_bc = nullptr;
    thermal_bc_constraint = 0;
  }
  int index;
  grid *gz_top;
  grid *gz_bottom;
  grid *gx_right;
  grid *gx_left;
  grid *gy_up;
  grid *gy_down;
  dimension dim;
  int boundary_type;       // 邊界平面: -1;not B, 0;xL, 1;xR, 2;yD, 3;yU, 4;zB, 5;zT
  int boundary_constraint; // 邊界條件種類: 0:no, 1:constant displacement, 2:constant force, 3: mixed

  material *mat;
  point org;
  SADB tempture;
  /*********thermal**************/
  thermal_constraint *thermal_bc;
  int thermal_bc_constraint; // 邊界條件種類: 0:no, 1:constant tempture, 2:constant heat flux, 3: convection
};
class power_block
{
public:
  power_block()
  {
    power_value = 0;
  }
  point org;
  dimension dim;
  double power_value;
};
class die
{
public:
  die() {};
  std::string name;
  point org;
  dimension dim;
  resoultion res;
  material mat;
  std::vector<std::vector<power_block *>> power_bls;
  //*****存die包含的grid在那些間隔中*****//
  std::pair<int, int> x_level;
  std::pair<int, int> y_level;
  std::pair<int, int> z_level;
};
class analysis_area
{
public:
  analysis_area()
  {
    Ambient_T = 25; // 一開始室溫假設25度
  };
  double Ambient_T;
  point org_point;
  dimension dim;
  material mat;
  resoultion global_res;
  resoultion min_gap;
  std::vector<die *> dies;
  std::vector<std::vector<std::vector<grid *>>> inside_grids;
  std::vector<std::vector<std::vector<grid *>>> boundary_grids; // 邊界平面: 0;xL, 1;xR, 2;yD, 3;yU, 4;zB, 5;zT
  std::map<SAINT, SAINT> x_cut_map;
  std::map<SAINT, SAINT> y_cut_map;
  std::map<SAINT, SAINT> z_cut_map;
  std::unordered_map<std::string, material> mat_map;
  //================funtion===================================//
  //*********grid construct**********//
  void init_grids(); // 在此funtion建立grid及賦值參數給grid
  void construct_uniform_cuts() ;
  void generate_grids();
  void grids_find_neighbor();
  void assign_material_to_grid();
  //************boundary condition setting**************//
  void setting_bc(
      int bc_case,      // 0:LEFT, 1:RIGHT, 2:BACK, 3:FRONT, 4:BOTTOM, 5:TOP
      std::string type, 
      std::vector<double> bc_value);
  //************power setting**************//
  void assign_bc_to_boundary_grids(
      int bc_case,      // 0:LEFT, 1:RIGHT, 2:BACK, 3:FRONT, 4:BOTTOM, 5:TOP
      std::string type, // Displacement, Force, Stress
      std::vector<double> bc_value);
  //************power setting**************//
  // void parse_CSV_file(
  //     std::string file_path,
  //     std::string design_name,
  //   );
  /****************/
  //========================================================//
};
