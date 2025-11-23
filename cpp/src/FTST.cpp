#include "thermal_analysis.hpp"
#include <chrono>
#include <tuple>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include "parser.hpp"
//===================================================================================//
double caculate_area_overlapratio(
	power_block *pbl,
	grid *Grid)
{
	double ratio = 1;
	double Ax = pbl->org.org_x;
	double Ay = pbl->org.org_y;
	double Aw = pbl->dim.dim_x;
	double Ah = pbl->dim.dim_y;
	double bx = Grid->org.org_x;
	double by = Grid->org.org_y;
	double bw = Grid->dim.dim_x;
	double bh = Grid->dim.dim_y;
	double overlap_width = std::max(0.0, std::min(Ax + Aw, bx + bw) - std::max(Ax, bx));
	double overlap_height = std::max(0.0, std::min(Ay + Ah, by + bh) - std::max(Ay, by));
	double overlap_area = overlap_width * overlap_height;
	double area_A = Aw * Ah;
	ratio = overlap_area / area_A;
	return ratio;
}
//===================================================================================//
void thermal_analysis::init_thermal_saprse_m()
{

	int inside_nodes_number = (AnalArea->x_cut_map.size() - 1) * (AnalArea->y_cut_map.size() - 1) * (AnalArea->z_cut_map.size() - 1);
	int x_boundary_nodes_number = (AnalArea->y_cut_map.size() - 1) * (AnalArea->z_cut_map.size() - 1);
	int y_boundary_nodes_number = (AnalArea->x_cut_map.size() - 1) * (AnalArea->z_cut_map.size() - 1);
	int z_boundary_nodes_number = (AnalArea->x_cut_map.size() - 1) * (AnalArea->y_cut_map.size() - 1);
	int total_nodes_number = inside_nodes_number + 2 * x_boundary_nodes_number + 2 * y_boundary_nodes_number + 2 * z_boundary_nodes_number; // 包含inside_grids和boundary_grids的總數量
	std::cout << "total inside nodes: " << inside_nodes_number << std::endl;
	std::cout << "total boundary nodes: " << 2 * x_boundary_nodes_number + 2 * y_boundary_nodes_number + 2 * z_boundary_nodes_number << std::endl;
	std::cout << "total nodes: " << total_nodes_number << std::endl;
	this->total_nodes_number = total_nodes_number;
	thermal_MNA_matrix.no_zero_term.resize(total_nodes_number, 7); // 先預設每個grid都是7個非0項(6個鄰居+自己)
	thermal_MNA_matrix.no_zero_term_column.resize(total_nodes_number);
	thermal_MNA_matrix.no_zero_term_value.resize(total_nodes_number);
	thermal_b_matrix.resize(total_nodes_number);
	for (int k = 0; k < AnalArea->boundary_grids.size(); k++)
	{

		for (int j = 0; j < AnalArea->boundary_grids[k].size(); j++)
		{
			for (int i = 0; i < AnalArea->boundary_grids[k][j].size(); i++)
			{
				auto Grid = AnalArea->boundary_grids[k][j][i];
				thermal_MNA_matrix.no_zero_term[Grid->index] = 2;
			}
		}
	}

	for (int i = 0; i < total_nodes_number; i++)
	{
		thermal_MNA_matrix.no_zero_term_column[i].resize(thermal_MNA_matrix.no_zero_term[i], -1);
		thermal_MNA_matrix.no_zero_term_value[i].resize(thermal_MNA_matrix.no_zero_term[i], 0);
	}

	std::cout << "init_thermal_saprse_m: " << std::endl;
}
//===================================================================================//
void thermal_analysis::stamp_thermal_MNA()
{
	auto start = std::chrono::high_resolution_clock::now();

	for (int k = 0; k < AnalArea->inside_grids.size(); k++)
	{
		for (int j = 0; j < AnalArea->inside_grids[k].size(); j++)
		{
			for (int i = 0; i < AnalArea->inside_grids[k][j].size(); i++)
			{
				auto Grid = AnalArea->inside_grids[k][j][i];
				int no_zero_index = 1; // 第0個都固定給對角線元素
				thermal_MNA_matrix.no_zero_term_column[Grid->index][0] = Grid->index;
				/******************************gx_left**********************************/
				if (Grid->gx_left->thermal_bc_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_x, 1, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gx_left->mat->k.k_x, 1, Grid->gx_left);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gx_left, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_x, 1, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gx_left, coefficient, no_zero_index);
					no_zero_index++;
				}
				/******************************gx_right**********************************/
				if (Grid->gx_right->thermal_bc_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_x, 1, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gx_right->mat->k.k_x, 1, Grid->gx_right);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gx_right, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_x, 1, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gx_right, coefficient, no_zero_index);
					no_zero_index++;
				}
				/******************************gy_down**********************************/
				if (Grid->gy_down->thermal_bc_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_y, 2, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gy_down->mat->k.k_y, 2, Grid->gy_down);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gy_down, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_y, 2, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gy_down, coefficient, no_zero_index);
					no_zero_index++;
				}
				/******************************gy_up**********************************/
				if (Grid->gy_up->thermal_bc_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_y, 2, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gy_up->mat->k.k_y, 2, Grid->gy_up);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gy_up, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_y, 2, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gy_up, coefficient, no_zero_index);
					no_zero_index++;
				}
				/******************************gz_bottom**********************************/
				if (Grid->gz_bottom->thermal_bc_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_z, 3, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gz_bottom->mat->k.k_z, 3, Grid->gz_bottom);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gz_bottom, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_z, 3, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gz_bottom, coefficient, no_zero_index);
					no_zero_index++;
				}
				/******************************gz_top**********************************/
				if (Grid->gz_top->boundary_constraint == 0)
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_z, 3, Grid);
					SADB R_neighbor = caculate_thermal_resitance(Grid->gz_top->mat->k.k_z, 3, Grid->gz_top);
					SADB coefficient = 1 / (R_itself + R_neighbor);
					stamp_resitance_to_sparsem(Grid, Grid->gz_top, coefficient, no_zero_index);
					no_zero_index++;
				}
				else
				{
					SADB R_itself = caculate_thermal_resitance(Grid->mat->k.k_z, 3, Grid);
					SADB coefficient = 1 / R_itself;
					stamp_resitance_to_sparsem(Grid, Grid->gz_top, coefficient, no_zero_index);
					no_zero_index++;
				}
			}
		}
	}
	for (int k = 0; k < AnalArea->boundary_grids.size(); k++)
	{

		for (int j = 0; j < AnalArea->boundary_grids[k].size(); j++)
		{
			for (int i = 0; i < AnalArea->boundary_grids[k][j].size(); i++)
			{
				auto Grid = AnalArea->boundary_grids[k][j][i];
				if (Grid->thermal_bc_constraint == 1)
					continue;
				double area = 0;
				if (k == 0 || k == 1)
					area = Grid->dim.dim_y * Grid->dim.dim_z * unit * unit;
				else if (k == 2 || k == 3)
					area = Grid->dim.dim_x * Grid->dim.dim_z * unit * unit;
				else if (k == 4 || k == 5)
					area = Grid->dim.dim_x * Grid->dim.dim_y * unit * unit;
				if (Grid->thermal_bc_constraint == 2)
				{
					thermal_b_matrix[Grid->index] += Grid->thermal_bc->heat_flux * area;
				}
				else
				{
					thermal_MNA_matrix.no_zero_term_value[Grid->index][0] += Grid->thermal_bc->convection * area;
					thermal_b_matrix[Grid->index] += Grid->thermal_bc->convection * area * AnalArea->Ambient_T;
				}
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "stamp_thermal_MNA執行時間: " << elapsed.count() << " 毫秒" << std::endl;
}
void thermal_analysis::stamp_power_to_b_matrix_ver2() // 笨方法，但能work
{
	std::ofstream fout;
	std::string file_path = "C:\\Users\\jenyi\\Desktop\\input_file_test\\answer_result.txt";
	fout.open(file_path);
	auto start = std::chrono::high_resolution_clock::now();
	for (int die_number = 0; die_number < AnalArea->dies.size(); die_number++)
	{
		auto die = AnalArea->dies[die_number];
		fout << die->name << std::endl;
		if (die->power_bls.size() == 0)
			continue;
		std::vector<power_block *> pb_sets;
		for (int i = 0; i < die->power_bls.size(); i++)
		{
			for (int j = 0; j < die->power_bls[i].size(); j++)
			{
				pb_sets.push_back(die->power_bls[i][j]);
			}
		}
#pragma omp parallel for num_threads(numThreads)
		for (int k = die->z_level.first; k < die->z_level.second; k++)
		{
			for (int j = die->y_level.first; j < die->y_level.second; j++)
			{
				for (int i = die->x_level.first; i < die->x_level.second; i++)
				{
					auto Grid = AnalArea->inside_grids[k][j][i];
					std::vector<power_block *> overlap_plbs;
					for (auto overlap_bl : pb_sets)
					{
						if (((overlap_bl->org.org_x <= Grid->org.org_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x > Grid->org.org_x) || (overlap_bl->org.org_x < Grid->org.org_x + Grid->dim.dim_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x >= Grid->org.org_x + Grid->dim.dim_x) || (overlap_bl->org.org_x >= Grid->org.org_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x <= Grid->org.org_x + Grid->dim.dim_x)) && ((overlap_bl->org.org_y <= Grid->org.org_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y > Grid->org.org_y) || (overlap_bl->org.org_y < Grid->org.org_y + Grid->dim.dim_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y >= Grid->org.org_y + Grid->dim.dim_y) || (overlap_bl->org.org_y >= Grid->org.org_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y <= Grid->org.org_y + Grid->dim.dim_y)))
						{
							overlap_plbs.push_back(overlap_bl);
						}
					}
					for (int id = 0; id < overlap_plbs.size(); id++)
					{
						double ratio = caculate_area_overlapratio(overlap_plbs[id], Grid);
						thermal_b_matrix[Grid->index] += (ratio * overlap_plbs[id]->power_value);
					}
					fout << Grid->index << " " << Grid->org.org_x + Grid->dim.dim_x << " " << Grid->org.org_y + Grid->dim.dim_y << " " << thermal_b_matrix[Grid->index] << std::endl;
				}
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "stamp_power_to_b_matrix執行時間: " << elapsed.count() << " 毫秒" << std::endl;
}
//===================================================================================//
void thermal_analysis::stamp_power_to_b_matrix_ver2() // 笨方法，但能work
{
	std::ofstream fout;
	std::string file_path = "../../output/answer_result.txt";
	fout.open(file_path);
	auto start = std::chrono::high_resolution_clock::now();
	for (int die_number = 0; die_number < AnalArea->dies.size(); die_number++)
	{
		auto die = AnalArea->dies[die_number];
		fout << die->name << std::endl;
		if (die->power_bls.size() == 0)
			continue;
		std::vector<power_block *> pb_sets;
		for (int i = 0; i < die->power_bls.size(); i++)
		{
			for (int j = 0; j < die->power_bls[i].size(); j++)
			{
				pb_sets.push_back(die->power_bls[i][j]);
			}
		}
// #pragma omp parallel for num_threads(numThreads)
		for (int k = die->z_level.first; k < die->z_level.second; k++)
		{
			for (int j = die->y_level.first; j < die->y_level.second; j++)
			{
				for (int i = die->x_level.first; i < die->x_level.second; i++)
				{
					auto Grid = AnalArea->inside_grids[k][j][i];
					std::vector<power_block *> overlap_plbs;
					for (auto overlap_bl : pb_sets)
					{
						if (((overlap_bl->org.org_x <= Grid->org.org_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x > Grid->org.org_x) || (overlap_bl->org.org_x < Grid->org.org_x + Grid->dim.dim_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x >= Grid->org.org_x + Grid->dim.dim_x) || (overlap_bl->org.org_x >= Grid->org.org_x && overlap_bl->org.org_x + overlap_bl->dim.dim_x <= Grid->org.org_x + Grid->dim.dim_x)) && ((overlap_bl->org.org_y <= Grid->org.org_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y > Grid->org.org_y) || (overlap_bl->org.org_y < Grid->org.org_y + Grid->dim.dim_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y >= Grid->org.org_y + Grid->dim.dim_y) || (overlap_bl->org.org_y >= Grid->org.org_y && overlap_bl->org.org_y + overlap_bl->dim.dim_y <= Grid->org.org_y + Grid->dim.dim_y)))
						{
							overlap_plbs.push_back(overlap_bl);
						}
					}
					for (int id = 0; id < overlap_plbs.size(); id++)
					{
						double ratio = caculate_area_overlapratio(overlap_plbs[id], Grid);
						thermal_b_matrix[Grid->index] += (ratio * overlap_plbs[id]->power_value);
					}
					fout << Grid->index << " " << Grid->org.org_x + Grid->dim.dim_x << " " << Grid->org.org_y + Grid->dim.dim_y << " " << thermal_b_matrix[Grid->index] << std::endl;
				}
			}
		}
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> elapsed = end - start;
	std::cout << "stamp_power_to_b_matrix執行時間: " << elapsed.count() << " 毫秒" << std::endl;
}
// void thermal_analysis::stamp_power_to_b_matrix_ver5()
// {
// 	auto start = std::chrono::high_resolution_clock::now();

// 	std::ofstream fout("C:\\Users\\jenyi\\Desktop\\input_file_test\\answer_result_tiled.txt");
// 	if (!fout.is_open())
// 	{
// 		std::cerr << "Cannot open output file!\n";
// 		return;
// 	}

// 	for (int die_number = 0; die_number < AnalArea->dies.size(); die_number++)
// 	{
// 		auto die = AnalArea->dies[die_number];
// 		fout << die->name << std::endl;

// 		if (die->power_bls.empty())
// 			continue;

// 		// === 攤平成一維 power block 集合 ===
// 		std::vector<power_block *> pb_sets;
// 		for (auto &row : die->power_bls)
// 			for (auto *pb : row)
// 				pb_sets.push_back(pb);

// 		// === 取出非均勻網格對應表 ===
// 		auto &x_cut_map = AnalArea->x_cut_map; // key: 實際座標, value: grid index
// 		auto &y_cut_map = AnalArea->y_cut_map;

// 		for (int k = die->z_level.first; k < die->z_level.second; k++)
// 		{
// 			for (auto *p : pb_sets)
// 			{
// 				double llx = p->org.org_x;
// 				double lly = p->org.org_y;
// 				double urx = p->org.org_x + p->dim.dim_x;
// 				double ury = p->org.org_y + p->dim.dim_y;

// 				// === 找出覆蓋範圍的左下、右上 index ===
// 				// upper_bound -> 第一個比 key 大的元素
// 				auto it_x_start = x_cut_map.upper_bound(llx);
// 				auto it_x_end = x_cut_map.upper_bound(urx);
// 				auto it_y_start = y_cut_map.upper_bound(lly);
// 				auto it_y_end = y_cut_map.upper_bound(ury);

// 				int start_x = 0;
// 				int end_x = die->x_level.second;
// 				int start_y = 0;
// 				int end_y = die->y_level.second;

// 				if (it_x_start != x_cut_map.begin())
// 					start_x = std::prev(it_x_start)->second;
// 				if (it_x_end != x_cut_map.end())
// 					end_x = it_x_end->second;
// 				if (it_y_start != y_cut_map.begin())
// 					start_y = std::prev(it_y_start)->second;
// 				if (it_y_end != y_cut_map.end())
// 					end_y = it_y_end->second;

// 				// === 掃描整個矩形範圍（左上→右下） ===
// 				for (int j = start_y; j < end_y; ++j)
// 				{
// 					for (int i = start_x; i < end_x; ++i)
// 					{
// 						auto Grid = AnalArea->inside_grids[k][j][i];

// 						// 檢查是否與 power block 有實際重疊
// 						if (Grid->org.org_x + Grid->dim.dim_x <= llx)
// 							continue;
// 						if (Grid->org.org_y + Grid->dim.dim_y <= lly)
// 							continue;
// 						if (Grid->org.org_x >= urx)
// 							continue;
// 						if (Grid->org.org_y >= ury)
// 							continue;

// 						double ratio = caculate_area_overlapratio(p, Grid);
// 						if (ratio <= 0.0)
// 							continue;

// 						double local_power = ratio * p->power_value;

// 						thermal_b_matrix[Grid->index] += local_power;

// 						{
// 							fout << Grid->index << " "
// 								 << Grid->org.org_x + Grid->dim.dim_x << " "
// 								 << Grid->org.org_y + Grid->dim.dim_y << " "
// 								 << thermal_b_matrix[Grid->index] << "\n";
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}

// 	fout.close();

// 	auto end = std::chrono::high_resolution_clock::now();
// 	std::chrono::duration<double, std::milli> elapsed = end - start;
// 	std::cout << "stamp_power_to_b_matrix_ver5 : " << elapsed.count() << " ms\n";
// }

// //===================================================================================//
// void thermal_analysis::stamp_power_to_b_matrix_ver6()
// {

// 	if (!fout.is_open())
// 	{
// 		std::cerr << "Cannot open output file!\n";
// 		return;
// 	}

// 	// === 掃所有 die ===
// 	for (int die_number = 0; die_number < AnalArea->dies.size(); die_number++)
// 	{
// 		auto die = AnalArea->dies[die_number];
// 		fout << die->name << std::endl;

// 		if (die->power_bls.empty())
// 			continue;

// 		// === 攤平成一維 power block 集合 ===
// 		std::vector<power_block *> pb_sets;
// 		for (auto &row : die->power_bls)
// 			for (auto *pb : row)
// 				pb_sets.push_back(pb);

// 		// === 非均勻網格映射表 ===
// 		auto &x_cut_map = AnalArea->x_cut_map; // key: 實際座標, value: grid index
// 		auto &y_cut_map = AnalArea->y_cut_map;

// 		const int TILE_X = 32; // tile 大小可調 (8~16 效果較佳)
// 		const int TILE_Y = 32;

// 		// === 平行處理每個 z 層 ===

// 		for (int k = die->z_level.first; k < die->z_level.second; k++)
// 		{
// 			for (auto *p : pb_sets)
// 			{
// 				double llx = p->org.org_x;
// 				double lly = p->org.org_y;
// 				double urx = p->org.org_x + p->dim.dim_x;
// 				double ury = p->org.org_y + p->dim.dim_y;

// 				// === 利用 map 查 grid 範圍 ===
// 				auto it_x_start = x_cut_map.upper_bound(llx);
// 				auto it_x_end = x_cut_map.upper_bound(urx);
// 				auto it_y_start = y_cut_map.upper_bound(lly);
// 				auto it_y_end = y_cut_map.upper_bound(ury);

// 				int start_x = 0, end_x = die->x_level.second;
// 				int start_y = 0, end_y = die->y_level.second;

// 				if (it_x_start != x_cut_map.begin())
// 					start_x = std::prev(it_x_start)->second;
// 				if (it_x_end != x_cut_map.end())
// 					end_x = it_x_end->second;
// 				if (it_y_start != y_cut_map.begin())
// 					start_y = std::prev(it_y_start)->second;
// 				if (it_y_end != y_cut_map.end())
// 					end_y = it_y_end->second;

// 				// === Tiling 掃描整個矩形區域 ===
// 				for (int ty = start_y; ty < end_y; ty += TILE_Y)
// 				{
// 					for (int tx = start_x; tx < end_x; tx += TILE_X)
// 					{
// 						int y_end_tile = min(ty + TILE_Y, end_y);
// 						int x_end_tile = min(tx + TILE_X, end_x);

// 						// === 掃 tile 內的所有 grid ===
// 						for (int j = ty; j < y_end_tile; ++j)
// 						{
// 							for (int i = tx; i < x_end_tile; ++i)
// 							{
// 								auto Grid = AnalArea->inside_grids[k][j][i];

// 								// 快速排除沒重疊的格子
// 								if (Grid->org.org_x + Grid->dim.dim_x <= llx)
// 									continue;
// 								if (Grid->org.org_y + Grid->dim.dim_y <= lly)
// 									continue;
// 								if (Grid->org.org_x >= urx)
// 									continue;
// 								if (Grid->org.org_y >= ury)
// 									continue;

// 								// 計算重疊面積比例
// 								double ratio = caculate_area_overlapratio(p, Grid);
// 								if (ratio <= 0.0)
// 									continue;

// 								double local_power = ratio * p->power_value;

// #pragma omp atomic
// 								thermal_b_matrix[Grid->index] += local_power;

// 								// Debug 輸出 (保留格式與 ver2 一致)
// #pragma omp critical
// 								{
// 									fout << Grid->index << " "
// 										 << Grid->org.org_x + Grid->dim.dim_x << " "
// 										 << Grid->org.org_y + Grid->dim.dim_y << " "
// 										 << thermal_b_matrix[Grid->index] << "\n";
// 								}
// 							}
// 						}
// 					}
// 				}
// 			}
// 		}
// 	}

// 	fout.close();

// }

//===================================================================================//
void thermal_analysis::stamp_resitance_to_sparsem(
	grid *Grid,
	grid *neighbor,
	SADB coefficient,
	int nonzero_index)
{
	thermal_MNA_matrix.no_zero_term_column[Grid->index][nonzero_index] = neighbor->index;
	thermal_MNA_matrix.no_zero_term_value[Grid->index][0] += coefficient;
	if (neighbor->thermal_bc_constraint != 1)
	{
		thermal_MNA_matrix.no_zero_term_value[Grid->index][nonzero_index] = -coefficient;
		if (neighbor->thermal_bc_constraint != 0) // 為非定溫的邊界點
		{
			thermal_MNA_matrix.no_zero_term_column[neighbor->index][0] = neighbor->index;
			thermal_MNA_matrix.no_zero_term_column[neighbor->index][1] = Grid->index;
			thermal_MNA_matrix.no_zero_term_value[neighbor->index][0] += coefficient;
			thermal_MNA_matrix.no_zero_term_value[neighbor->index][1] = -coefficient;
		}
	}
	else // 為定溫的邊界點
	{
		thermal_MNA_matrix.no_zero_term_column[neighbor->index][0] = neighbor->index;
		thermal_MNA_matrix.no_zero_term_column[neighbor->index][1] = Grid->index;
		thermal_MNA_matrix.no_zero_term_value[neighbor->index][0] = 1;
		thermal_b_matrix[Grid->index] += coefficient * neighbor->thermal_bc->constant_tempture;
		thermal_b_matrix[neighbor->index] += neighbor->thermal_bc->constant_tempture;
	}
}
//===================================================================================//
void thermal_analysis::Thermal_Analysis()
{
	init_thermal_saprse_m();
	stamp_power_to_b_matrix_ver2();
	stamp_thermal_MNA();
}
//=============================================================================//
void thermal_analysis::assign_tempture_to_grids()
{

	for (int k = 0; k < AnalArea->inside_grids.size(); k++)
	{
		for (int j = 0; j < AnalArea->inside_grids[k].size(); j++)
		{
			for (int i = 0; i < AnalArea->inside_grids[k][j].size(); i++)
			{
				auto Grid = AnalArea->inside_grids[k][j][i];
				Grid->tempture = thermal_result[Grid->index];
			}
		}
	}
	for (int k = 0; k < 6; k++)
	{

		for (int j = 0; j < AnalArea->boundary_grids[k].size(); j++)
		{
			for (int i = 0; i < AnalArea->boundary_grids[k][j].size(); i++)
			{
				auto Grid = AnalArea->boundary_grids[k][j][i];
				Grid->tempture = thermal_result[Grid->index];
			}
		}
	}
}
//=============================================================================//
SADB caculate_thermal_resitance(
	SADB mat_value,
	int area_case, // 1: yz面積, 2: xz面積, 3: xy面積
	grid *Grid)
{
	if (area_case == 1)
		return (Grid->dim.dim_x / 2.0) / ((Grid->dim.dim_y * Grid->dim.dim_z * unit) * (mat_value));
	else if (area_case == 2)
		return (Grid->dim.dim_y / 2.0) / ((Grid->dim.dim_x * Grid->dim.dim_z * unit) * (mat_value));
	else if (area_case == 3)
		return (Grid->dim.dim_z / 2.0) / ((Grid->dim.dim_x * Grid->dim.dim_y * unit) * (mat_value));
}
//========================================================//