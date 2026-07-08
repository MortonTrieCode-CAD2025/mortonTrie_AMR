// /**
//  * @file D3Q19_BC_Manager.cpp
//  * @brief 
//  * @date 2023-07-18
//  */

#include "BC/D3Q19_BC_Manager.hpp"
// #include "Lat_Manager.h"
// // #include "Collision_Model.h"
// #include "LBM_Manager.hpp"
// #include "Morton_Assist.h"

#ifdef BC_NO_GHOST
#if (C_DIMS == 3)
 #if (C_Q == 19)
constexpr std::array<std::array<int, 14>,26> D3Q19_BC_Strategy::bdry_dirs;
 #elif (C_Q == 27)
constexpr std::array<std::array<int, 18>,C_Q-1> D3Q19_BC_Strategy::bdry_dirs;
 #endif
#endif
#endif

// void D3Q19_BC_Strategy::initial_Feq(D_mapLat* bc_ptr)
// {
//     D_Phy_DDF feq_temp[C_Q];
//     for (auto i_bc_lat : *bc_ptr)
//     {
//         D_morton lat_code = i_bc_lat.first;
//         calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
//         for (D_int i_q = 0; i_q < C_Q; ++i_q)
//         {
//             df_ptr->f[i_q].at(lat_code)  = feq_temp[i_q];
//             df_ptr->fcol[i_q].at(lat_code)  = feq_temp[i_q];
//         }
        
//     }
// }