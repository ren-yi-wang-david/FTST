// #include "parser.hpp"
// #include <unordered_set>
// #include <fstream>
// #include <chrono>
// #include <algorithm>

// void analysis_area::parse_CSV_file(die *die, std::ifstream &PI_file)
// {
//     std::string line;

//     // name → power
//     std::unordered_map<std::string, double> power_table;

//     // name → geometry (xs, ys, xe, ye)
//     struct Geo
//     {
//         double xs, ys, xe, ye;
//     };
//     std::unordered_map<std::string, Geo> geo_table;

//     // =========================================================
//     // STEP 1：讀 Power Table
//     // =========================================================

//     PI_file.clear();
//     PI_file.seekg(0, std::ios::beg);

//     while (std::getline(PI_file, line))
//     {
//         if (line.find("name,power") != std::string::npos)
//             break;
//     }

//     while (std::getline(PI_file, line))
//     {
//         if (line.size() == 0)
//             continue;
//         if (line[0] == '#')
//             continue; // 跳註解

//         std::stringstream ss(line);
//         std::string token;
//         std::vector<std::string> t;

//         while (std::getline(ss, token, ','))
//             t.push_back(token);
//         if (t.size() < 2)
//             continue;

//         std::string name = t[0];
//         if (name == "")
//             continue;
//         if (!std::isdigit(name[0]) && !std::isalpha(name[0]))
//             continue;

//         double power = std::stod(t[1]);
//         power_table[name] = power;
//     }

//     // =========================================================
//     // STEP 2：讀 Blocks Prism Geometry
//     // =========================================================

//     PI_file.clear();
//     PI_file.seekg(0, std::ios::beg);

//     while (std::getline(PI_file, line))
//     {
//         if (line.find("name,xs,ys,zs") != std::string::npos)
//             break;
//     }

//     while (std::getline(PI_file, line))
//     {
//         if (line.size() == 0)
//             continue;

//         std::stringstream ss(line);
//         std::string token;
//         std::vector<std::string> t;

//         while (std::getline(ss, token, ','))
//             t.push_back(token);
//         if (t.size() < 7)
//             continue;

//         std::string name = t[0];
//         if (name == "")
//             continue;
//         if (power_table.find(name) == power_table.end())
//             continue;

//         double xs = std::stod(t[1]);
//         double ys = std::stod(t[2]);
//         double xe = std::stod(t[4]);
//         double ye = std::stod(t[5]);

//         geo_table[name] = {xs, ys, xe, ye};
//     }

//     // =========================================================
//     // STEP 3：把兩個 table 合併建 power_block
//     // =========================================================

//     std::vector<power_block *> temp;

//     for (auto &kv : power_table)
//     {
//         const std::string &name = kv.first;

//         if (geo_table.find(name) == geo_table.end())
//             continue; // 無 geometry → 無法建 block

//         Geo g = geo_table[name];
//         double power = kv.second;

//         power_block *pb = new power_block();

//         pb->org.org_x = g.xs + die->org.org_x;
//         pb->org.org_y = g.ys + die->org.org_y;
//         pb->dim.dim_x = g.xe - g.xs;
//         pb->dim.dim_y = g.ye - g.ys;

//         pb->power_value = power;

//         temp.push_back(pb);
//     }

//     // =========================================================
//     // STEP 4：依照你原本邏輯切成 rows
//     // =========================================================

//     int row_number = 0;
//     for (auto *pb : temp)
//     {
//         if (pb->org.org_x + pb->dim.dim_x ==
//             die->org.org_x + die->dim.dim_x)
//         {
//             row_number++;
//         }
//     }

//     die->power_bls.resize(row_number);

//     int j = 0;
//     for (int i = 0; i < temp.size(); i++)
//     {
//         die->power_bls[j].push_back(temp[i]);
//         if (temp[i]->org.org_x + temp[i]->dim.dim_x >=
//             die->org.org_x + die->dim.dim_x)
//         {
//             j++;
//         }
//     }
// }
