/**
 * @file NonEquilibrumExtrapolation.cpp
 * @brief Non-equilibrum Extrapolation Method
 * @date 2023-07-27
 */

#include "BC/NonEquilibrumExtrapolation.hpp"
#include "Lat_Manager.hpp"
#include "Morton_Assist.h"
#include "LBM_Manager.hpp"
#include <tuple>

constexpr std::array<uint,9> nonEquilibrumExtrapolation_West::lat_bc_idx;
constexpr std::array<uint,9> nonEquilibrumExtrapolation_East::lat_bc_idx;
constexpr std::array<uint,3> nonEquilibrumExtrapolation_South::lat_bc_idx;
constexpr std::array<uint,3> nonEquilibrumExtrapolation_North::lat_bc_idx;
constexpr std::array<uint,1> nonEquilibrumExtrapolation_Top::lat_bc_idx;
constexpr std::array<uint,1> nonEquilibrumExtrapolation_Bot::lat_bc_idx;

void nonEquilibrumExtrapolation_West::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
    initial_Feq(westBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_West::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_West::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
    for (auto i_bc_lat : *westBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_x1(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), {0.,0.,0.}, feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_West::lat_bc_idx)
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif


}

void nonEquilibrumExtrapolation_East::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
    initial_Feq(eastBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_East::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_East::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
    for (auto i_bc_lat : *eastBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_x0(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), velocity_ptr->at(nbr_code), feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_East::lat_bc_idx) 
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif
}

void nonEquilibrumExtrapolation_South::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
    initial_Feq(southBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_South::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_South::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
    for (auto i_bc_lat : *southBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_y1(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), {0.,0.,0.}, feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_South::lat_bc_idx) 
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif
}

void nonEquilibrumExtrapolation_North::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
    initial_Feq(northBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_North::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_North::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
    for (auto i_bc_lat : *northBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_y0(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), {0.,0.,0.}, feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_North::lat_bc_idx) 
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif
}

void nonEquilibrumExtrapolation_Bot::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
    initial_Feq(botBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_Bot::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_Bot::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
    for (auto i_bc_lat : *botBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_z1(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), {0.,0.,0.}, feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_Bot::lat_bc_idx) 
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif
}

void nonEquilibrumExtrapolation_Top::initialBCStrategy()
{
#ifndef BC_NO_GHOST
    D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
    initial_Feq(topBc_ptr);
#else
    for (uint bc_idx : nonEquilibrumExtrapolation_Top::lat_bc_idx) {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        initial_Feq(bc_ptr);
    }
#endif
}

void nonEquilibrumExtrapolation_Top::applyBCStrategy()
{
#ifndef BC_NO_GHOST
 #ifndef eflow_stable_bc
    D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
    for (auto i_bc_lat : *topBc_ptr) {
        D_morton lat_code = i_bc_lat.first;
        D_morton nbr_code = Morton_Assist::find_z0(lat_code, 0);
        D_Phy_DDF feq[C_Q], nbr_feq[C_Q];
        calculateFeq_HeLuo(density_ptr->at(lat_code), {0.,0.,0.}, feq);
        calculateFeq_HeLuo(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), nbr_feq);
        for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            f_ptr[i_q].at(lat_code) = feq[i_q] + f_ptr[i_q].at(lat_code) - nbr_feq[i_q];
        }
    }
 #endif
#else
 #ifdef eflow_stable_bc
    for (auto bc_idx : nonEquilibrumExtrapolation_Top::lat_bc_idx) 
    {
        D_setLat* bc_ptr = &(Lat_Manager::pointer_me->lat_bc.at(bc_idx));
        auto dirs = D3Q19_BC_Strategy::bdry_dirs.at(bc_idx);
        for (D_morton lat_code : *bc_ptr) 
        {
            for (auto i_q : dirs) 
            {
                if (i_q < 0) break;
                D_morton nbr_code = Morton_Assist::find_neighbor(lat_code, 0, ex[i_q], ey[i_q], ez[i_q]);
                D_Phy_DDF feq_w = calculateFeq(density_ptr->at(lat_code), bc_velocity, i_q);
                D_Phy_DDF feq_n = calculateFeq(density_ptr->at(nbr_code), velocity_ptr->at(nbr_code), i_q);
                f_ptr[i_q].at(lat_code) = feq_w + f_ptr[i_q].at(nbr_code) - feq_n;
            } // for i_q
        } // for bc_lat
    }
 #endif
#endif
}
