#include"parser.hpp"
#include<unordered_set>
#include<fstream>
#include <chrono>
#include<algorithm>
void analysis_area::construct_uniform_cuts()
{
	std::cout << "Construct uniform grid..." << std::endl;

	x_cut_map.clear();
	y_cut_map.clear();
	z_cut_map.clear();

	//=============================
	// 1. 產生 X cut
	//=============================
	{
		SAINT start = org_point.org_x;
		SAINT end = org_point.org_x + dim.dim_x;
		SAINT step = global_res.res_x;

		int idx = 0;
		for (SAINT x = start; x <= end; x += step)
		{
			x_cut_map[idx++] = x;
		}

		// 若最後一刀沒切到終點 → 強制補齊
		if (x_cut_map.rbegin()->second != end)
		{
			x_cut_map[idx] = end;
		}
	}

	//=============================
	// 2. 產生 Y cut
	//=============================
	{

		SAINT start = org_point.org_y;
		SAINT end = org_point.org_y + dim.dim_y;
		SAINT step = global_res.res_y;
		int idx = 0;
		for (SAINT y = start; y <= end; y += step)
		{
			y_cut_map[idx++] = y;
		}

		if (y_cut_map.rbegin()->second != end)
		{
			y_cut_map[idx] = end;
		}
	}

	//=============================
	// 3. 產生 Z cut
	//=============================
	{
		SAINT start = org_point.org_z;
		SAINT end = org_point.org_z + dim.dim_z;
		SAINT step = global_res.res_z;

		int idx = 0;
		for (SAINT z = start; z <= end; z += step)
		{
			z_cut_map[idx++] = z;
		}

		if (z_cut_map.rbegin()->second != end)
		{
			z_cut_map[idx] = end;
		}
	}

	std::cout << "Uniform grid constructed: "
			  << x_cut_map.size() << " x, "
			  << y_cut_map.size() << " y, "
			  << z_cut_map.size() << " z cuts." << std::endl;
}

//==============================================//
void analysis_area::generate_grids()
{
  std::vector<SAINT> x_cut;
  std::vector<SAINT> y_cut;
  std::vector<SAINT> z_cut;
  x_cut.reserve(x_cut_map.size());
  y_cut.reserve(y_cut_map.size());
  z_cut.reserve(z_cut_map.size());
  for(auto x = x_cut_map.begin(); x != x_cut_map.end(); x++) x_cut.push_back(x->first);
  for(auto y = y_cut_map.begin(); y != y_cut_map.end(); y++) y_cut.push_back(y->first);
  for(auto z = z_cut_map.begin(); z != z_cut_map.end(); z++) z_cut.push_back(z->first);
  std::cout << (x_cut.size() - 1) << std::endl;
  std::cout << (y_cut.size() - 1) << std::endl;
  std::cout << (z_cut.size() - 1) << std::endl;
  std::vector<std::vector<std::vector<int>>> inside_grids_index_vector;//為了要讓生成grids可平行化，預先存下所有網格該在的index
  std::vector<std::vector<std::vector<int>>> boundary_grids_index_vector;//為了要讓生成grids可平行化，預先存下所有網格該在的index
  /*******************************/
  //先切內部的
  /*******************************/
  inside_grids.resize(z_cut.size() - 1);
  inside_grids_index_vector.resize(z_cut.size() - 1);
  for (int k = 0; k < z_cut.size() - 1; k++)
  {
	inside_grids[k].resize(y_cut.size() - 1);
	inside_grids_index_vector[k].resize(y_cut.size() - 1);
	for (int j = 0; j < y_cut.size() - 1;  j++)
	{
	  inside_grids[k][j].resize(x_cut.size() - 1);
	  inside_grids_index_vector[k][j].resize(x_cut.size() - 1);
	}
  }
  //*****為了要讓生成grids可平行化，事先存好index*******//
  int index = 0;
  for (int k = 0; k < (z_cut.size()) - 1; k++)
  {
	for (int j = 0; j < (y_cut.size()) - 1; j++)
	{
	  for (int i = 0; i < (x_cut.size()) - 1; i++)
	  {
		inside_grids_index_vector[k][j][i] = index;
		index++;
	  }
	}
  }
  //***********************************************//
 
  for (int k = 0; k < (z_cut.size()) - 1; k++)
  {
	for (int j = 0; j < (y_cut.size()) - 1; j++)
	{
	  for (int i = 0; i < (x_cut.size()) - 1; i++)
	  {
		grid* Grid = new grid();
		Grid->dim.dim_x = x_cut[i+1] - x_cut[i];
		Grid->dim.dim_y = y_cut[j+1] - y_cut[j];
		Grid->dim.dim_z = z_cut[k+1] - z_cut[k];
		Grid->org.org_x = x_cut[i];
		Grid->org.org_y = y_cut[j];
		Grid->org.org_z = z_cut[k];
		Grid->index = inside_grids_index_vector[k][j][i];
		Grid->mat = &mat_map["mold"];
		inside_grids[k][j][i] = Grid;
	  }
	}
  }
/*******************************/
//再切邊界面上的
/*******************************/
  boundary_grids.resize(6);
  boundary_grids_index_vector.resize(6);
  for (int bondary = 0; bondary < 6; bondary++)
  {
	if (bondary == 0 || bondary == 1)
	{
	  boundary_grids[bondary].resize(z_cut.size() - 1);
	  boundary_grids_index_vector[bondary].resize(z_cut.size() - 1);
	  for (size_t k = 0; k < z_cut.size() - 1; k++)
	  {
		boundary_grids[bondary][k].resize(y_cut.size() - 1);
		boundary_grids_index_vector[bondary][k].resize(y_cut.size() - 1);
	  }
	}
	else if (bondary == 2 || bondary == 3)
	{
	  boundary_grids[bondary].resize(z_cut.size() - 1);
	  boundary_grids_index_vector[bondary].resize(z_cut.size() - 1);
	  for (size_t k = 0; k < z_cut.size() - 1; k++)
	  {
		boundary_grids[bondary][k].resize(x_cut.size() - 1);
		boundary_grids_index_vector[bondary][k].resize(x_cut.size() - 1);
	  }
	}
	else if (bondary == 4 || bondary == 5)
	{
	  boundary_grids[bondary].resize(y_cut.size() - 1);
	  boundary_grids_index_vector[bondary].resize(y_cut.size() - 1);
	  for (size_t j = 0; j < y_cut.size() - 1; j++)
	  {
		boundary_grids[bondary][j].resize(x_cut.size() - 1);
		boundary_grids_index_vector[bondary][j].resize(x_cut.size() - 1);
	  }
	}
  }
  //*****為了要讓生成grids可平行化，事先存好index*******//
  for (int boundary = 0; boundary < 6; boundary++)
  {
	for ( int j = 0;j < boundary_grids_index_vector[boundary].size();j++)
	{
	  for (int i = 0; i < boundary_grids_index_vector[boundary][j].size(); i++)
	  {
		boundary_grids_index_vector[boundary][j][i] = index;
		index++;
	  }
	}
  }
  //***********************************************//
  for (int boundary = 0; boundary < 6; boundary++)
  {
	if (boundary == 0 || boundary ==1)
	{

	  for (int k = 0; k < z_cut.size() - 1; k++)
	  {
		for (int j = 0; j < y_cut.size() - 1; j++)
		{
		  grid* Grid = new grid();
		  Grid->dim.dim_y = y_cut[j + 1] - y_cut[j];
		  Grid->dim.dim_z = z_cut[k + 1] - z_cut[k];
		  Grid->org.org_y = y_cut[j];
		  Grid->org.org_z = z_cut[k];
		  Grid->index = boundary_grids_index_vector[boundary][k][j];
		  Grid->mat = &mat_map["mold"];
		  boundary_grids[boundary][k][j] = Grid;
		  if (boundary == 0)
		  {
			Grid->org.org_x = org_point.org_x;
			Grid->boundary_type = 0;
		  }
		  else
		  {
			Grid->org.org_x = org_point.org_x + dim.dim_x;
			Grid->boundary_type = 1;
		  }
		}
	  }
	}
	else if (boundary == 2 || boundary == 3)
	{

	  for (int k = 0; k < z_cut.size() - 1; k++)
	  {
		for (int i = 0; i < x_cut.size() - 1; i++)
		{
		  grid* Grid = new grid();
		  Grid->dim.dim_x = x_cut[i + 1] - x_cut[i];
		  Grid->dim.dim_z = z_cut[k + 1] - z_cut[k];
		  Grid->org.org_x = x_cut[i];
		  Grid->org.org_z = z_cut[k];
		  Grid->index = boundary_grids_index_vector[boundary][k][i];
		  Grid->mat = &mat_map["mold"];
		  boundary_grids[boundary][k][i] = Grid;
		  if (boundary == 2)
		  {
			Grid->org.org_y = org_point.org_y;
			Grid->boundary_type = 2;
		  }
		  else
		  {
			Grid->org.org_y = org_point.org_y + dim.dim_y;
			Grid->boundary_type = 3;
		  }
		}
	  }
	}
	else if (boundary == 4 || boundary == 5)
	{

	  for (int j = 0; j < y_cut.size() - 1; j++)
	  {
		for (size_t i = 0; i < x_cut.size() - 1; i++)
		{
		  grid* Grid = new grid();
		  Grid->dim.dim_x = x_cut[i + 1] - x_cut[i];
		  Grid->dim.dim_y = y_cut[j + 1] - y_cut[j];
		  Grid->org.org_x = x_cut[i];
		  Grid->org.org_y = y_cut[j];
		  Grid->index = boundary_grids_index_vector[boundary][j][i];
		  Grid->mat = &mat_map["mold"];
		  boundary_grids[boundary][j][i] = Grid;
		  if (boundary == 4)
		  {
			Grid->org.org_z = org_point.org_z;
			Grid->boundary_type = 4;
		  }
		  else
		  {
			Grid->org.org_z = org_point.org_z + dim.dim_z;
			Grid->boundary_type = 5;
		  }
		}
	  }
	}
  }
 
/*******************************/
}
//===========================================//
void analysis_area::grids_find_neighbor()
{
 
  size_t x_number = inside_grids[0][0].size();
  size_t y_number = inside_grids[0].size();
  size_t z_number = inside_grids.size();
  /*********************************************/
  //先處理內部grid的鄰居
  /*********************************************/

  for (int layer = 0; layer < z_number; layer++)
  {
	for (int j = 0; j < y_number; j++)
	{
	  for (int i = 0; i < x_number; i++)
	  {
		if (layer == 0)
		  inside_grids[layer][j][i]->gz_top = inside_grids[layer + 1][j][i];
		else if (layer == z_number - 1)
		  inside_grids[layer][j][i]->gz_bottom = inside_grids[layer - 1][j][i];
		else
		{
		  inside_grids[layer][j][i]->gz_top = inside_grids[layer + 1][j][i];
		  inside_grids[layer][j][i]->gz_bottom = inside_grids[layer - 1][j][i];
		}

		if (j == 0)
		  inside_grids[layer][j][i]->gy_up = inside_grids[layer][j + 1][i];
		else if (j == y_number - 1)
		  inside_grids[layer][j][i]->gy_down = inside_grids[layer][j - 1][i];
		else
		{
		  inside_grids[layer][j][i]->gy_up = inside_grids[layer][j + 1][i];
		  inside_grids[layer][j][i]->gy_down = inside_grids[layer][j - 1][i];
		}

		if (i == 0)
		  inside_grids[layer][j][i]->gx_right = inside_grids[layer][j][i + 1];
		else if (i == x_number - 1)
		  inside_grids[layer][j][i]->gx_left = inside_grids[layer][j][i - 1];
		else
		{
		  inside_grids[layer][j][i]->gx_right = inside_grids[layer][j][i + 1];
		  inside_grids[layer][j][i]->gx_left = inside_grids[layer][j][i - 1];
		}
	  }
	}
  }
  /*********************************************/
  //再處理外部grid的鄰居
  /*********************************************/
  for (int boundary = 0; boundary < 6; boundary++)
  {
	if (boundary == 0 || boundary == 1)
	{
	  
	  for (int k = 0; k < z_number; k++)
	  {
		for (int j = 0; j < y_number; j++)
		{
		  if (k == 0)
			boundary_grids[boundary][k][j]->gz_top = boundary_grids[boundary][k + 1][j];
		  else if (k == z_number - 1)
			boundary_grids[boundary][k][j]->gz_bottom = boundary_grids[boundary][k - 1][j];
		  else
		  {
			boundary_grids[boundary][k][j]->gz_top = boundary_grids[boundary][k + 1][j];
			boundary_grids[boundary][k][j]->gz_bottom = boundary_grids[boundary][k - 1][j];
		  }
		  if (j == 0)
			boundary_grids[boundary][k][j]->gy_up = boundary_grids[boundary][k][j + 1];
		  else if (j == y_number - 1)
			boundary_grids[boundary][k][j]->gy_down = boundary_grids[boundary][k][j - 1];
		  else
		  {
			boundary_grids[boundary][k][j]->gy_up = boundary_grids[boundary][k][j + 1];
			boundary_grids[boundary][k][j]->gy_down = boundary_grids[boundary][k][j - 1];
		  }

		  if (boundary == 0)
		  {
			boundary_grids[boundary][k][j]->gx_right = inside_grids[k][j][0];
			inside_grids[k][j][0]->gx_left = boundary_grids[boundary][k][j];
		  }
		  else if (boundary == 1)
		  {
			boundary_grids[boundary][k][j]->gx_left = inside_grids[k][j][x_number - 1];
			inside_grids[k][j][x_number - 1]->gx_right = boundary_grids[boundary][k][j];
		  }
		}
	  }
	}
	else if (boundary == 2 || boundary == 3)
	{
	  
	  for (int k = 0; k < z_number; k++)
	  {
		for (int i = 0; i < x_number; i++)
		{
		  if (k == 0)
			boundary_grids[boundary][k][i]->gz_top = boundary_grids[boundary][k + 1][i];
		  else if (k == z_number - 1)
			boundary_grids[boundary][k][i]->gz_bottom = boundary_grids[boundary][k - 1][i];
		  else
		  {
			boundary_grids[boundary][k][i]->gz_top = boundary_grids[boundary][k + 1][i];
			boundary_grids[boundary][k][i]->gz_bottom = boundary_grids[boundary][k - 1][i];
		  }
		  if (i == 0)
			boundary_grids[boundary][k][i]->gx_left = boundary_grids[boundary][k][i + 1];
		  else if (i == x_number - 1)
			boundary_grids[boundary][k][i]->gx_right = boundary_grids[boundary][k][i - 1];
		  else
		  {
			boundary_grids[boundary][k][i]->gx_left = boundary_grids[boundary][k][i + 1];
			boundary_grids[boundary][k][i]->gx_right = boundary_grids[boundary][k][i - 1];
		  }

		  if (boundary == 2)
		  {
			boundary_grids[boundary][k][i]->gy_up = inside_grids[k][0][i];
			inside_grids[k][0][i]->gy_down = boundary_grids[boundary][k][i];
		  }
		  else if (boundary == 3)
		  {
			boundary_grids[boundary][k][i]->gy_down = inside_grids[k][y_number - 1][i];
			inside_grids[k][y_number - 1][i]->gy_up = boundary_grids[boundary][k][i];
		  }
		}
	  }
	}
	else if (boundary == 4 || boundary == 5)
	{
	  
	  for (int j = 0; j < y_number; j++)
	  {
		for (int i = 0; i < x_number; i++)
		{
		  if (i == 0)
			boundary_grids[boundary][j][i]->gx_left = boundary_grids[boundary][j][i + 1];
		  else if (i == x_number - 1)
			boundary_grids[boundary][j][i]->gx_right = boundary_grids[boundary][j][i - 1];
		  else
		  {
			boundary_grids[boundary][j][i]->gx_left = boundary_grids[boundary][j][i + 1];
			boundary_grids[boundary][j][i]->gx_right = boundary_grids[boundary][j][i - 1];
		  }
		  if (j == 0)
			boundary_grids[boundary][j][i]->gy_up = boundary_grids[boundary][j + 1][i];
		  else if (j == y_number - 1)
			boundary_grids[boundary][j][i]->gy_down = boundary_grids[boundary][j - 1][i];
		  else
		  {
			boundary_grids[boundary][j][i]->gy_up = boundary_grids[boundary][j + 1][i];
			boundary_grids[boundary][j][i]->gy_down = boundary_grids[boundary][j - 1][i];
		  }
		  if (boundary == 4)
		  {
			boundary_grids[boundary][j][i]->gz_top = inside_grids[0][j][i];
			inside_grids[0][j][i]->gz_bottom = boundary_grids[boundary][j][i];
		  }
		  else if (boundary == 5)
		  {
			boundary_grids[boundary][j][i]->gz_bottom = inside_grids[z_number - 1][j][i];
			inside_grids[z_number - 1][j][i]->gz_top = boundary_grids[boundary][j][i];
		  }
		}
	  }
	}
  }
}


//======================================================//
void analysis_area::init_grids()
{
	construct_uniform_cuts();
	generate_grids();
	grids_find_neighbor();
	std::cout << "finish initial grids" << std::endl;
	setting_bc();
	std::cout << "finish set boundary condition " << std::endl;
}
//=====================================================//
void analysis_area::setting_bc(
	int bc_case,
	std::string type, // Displacement, Force, Stress
	std::vector<double> bc_value)
{
  for(int i = 0; i < 6; i++){
	  assign_bc_to_boundary_grids(i, "htc", {0.0});
  }    
}
//====================================================================//
void analysis_area::assign_bc_to_boundary_grids(
	int bc_case, // 0:LEFT, 1:RIGHT, 2:BACK, 3:FRONT, 4:BOTTOM, 5:TOP
	std::string type,
	std::vector<double> bc_value)
{
  if (type == "Tempture")
  {
    thermal_constraint* thermal_bc = new thermal_constraint();
  
    for (int i = 0; i < boundary_grids[bc_case].size(); i++)
    {
      for (int j = 0; j < boundary_grids[bc_case][i].size(); j++)
      {
        boundary_grids[bc_case][i][j]->thermal_bc_constraint = 1;
        boundary_grids[bc_case][i][j]->thermal_bc = thermal_bc;
        boundary_grids[bc_case][i][j]->thermal_bc->constant_tempture = bc_value[0];
      }
    }
  }
  else if (type == "Heatflux")
  {
    thermal_constraint* thermal_bc = new thermal_constraint();

    for (int i = 0; i < boundary_grids[bc_case].size(); i++)
    {
      for (int j = 0; j < boundary_grids[bc_case][i].size(); j++)
      {
        boundary_grids[bc_case][i][j]->thermal_bc_constraint = 2;
        boundary_grids[bc_case][i][j]->thermal_bc = thermal_bc;
        boundary_grids[bc_case][i][j]->thermal_bc->heat_flux = bc_value[0];
      }
    }
  }
  else if (type == "htc")
  {
    thermal_constraint* thermal_bc = new thermal_constraint();
 
    for (int i = 0; i < boundary_grids[bc_case].size(); i++)
    {
      for (int j = 0; j < boundary_grids[bc_case][i].size(); j++)
      {
        boundary_grids[bc_case][i][j]->thermal_bc_constraint = 3;
        boundary_grids[bc_case][i][j]->thermal_bc = thermal_bc;
        boundary_grids[bc_case][i][j]->thermal_bc->convection = bc_value[0];
      }
    }
  }
  }
}