/**
 * @file D3Q19_BC_Class.h
 * @brief 
 * @date 2023-07-18
 */
#pragma once

#include "user.h"
#include <map>
#include <memory>
#include "Lat_Manager.hpp"
#include "LBM_Manager.hpp"
#include <tuple>

class D3Q19_BC_Strategy
{
protected:
    D_DDFSlot* f_ptr;
    _s_DDF* df_ptr;
    D_map_define<D_Phy_Rho>* density_ptr;
    D_map_define<D_uvw>* velocity_ptr;
    D_Phy_real CF_U_;
    // User* user_;
    D_Phy_real tau_;
    enum Bdry_Type : uint {
        BDRY_FACE_WEST  = 0,
        BDRY_FACE_EAST  = 1,
        BDRY_FACE_SOUTH = 2,
        BDRY_FACE_NORTH = 3,
        BDRY_FACE_BOT   = 4,
        BDRY_FACE_TOP   = 5,

        BDRY_EDGE_EN = 6,
        BDRY_EDGE_WN = 7,
        BDRY_EDGE_WS = 8,
        BDRY_EDGE_ES = 9,
        BDRY_EDGE_NT = 10,
        BDRY_EDGE_ST = 11,
        BDRY_EDGE_NB = 12,
        BDRY_EDGE_SB = 13,
        BDRY_EDGE_ET = 14,
        BDRY_EDGE_WT = 15,
        BDRY_EDGE_EB = 16,
        BDRY_EDGE_WB = 17,

        BDRY_CORNER_ENT = 18,
        BDRY_CORNER_WNT = 19,
        BDRY_CORNER_WST = 20,
        BDRY_CORNER_EST = 21,
        BDRY_CORNER_ENB = 22,
        BDRY_CORNER_WNB = 23,
        BDRY_CORNER_WSB = 24,
        BDRY_CORNER_ESB = 25 };

#ifdef BC_NO_GHOST
// The directions for the BC lattices,
//  the reserved directions pointing to existing neighbor lattices
#if (C_DIMS == 3)
 #if (C_Q == 19)
    static constexpr std::array<std::array<int, 14>,26>
    // std::array<std::array<int, 14>,26>  (outer = 6 faces + 12 edges + 8 corners = 26, NOT C_Q-1=18)
     bdry_dirs = {{{0,1,3,4,5,6,7,9,11,13,15,16,17,18}, //FACE_WEST[14]
                    {0,2,3,4,5,6,8,10,12,14,15,16,17,18}, //FACE_EAST[14]
                     {0,1,2,3,5,6,7,10,11,12,13,14,15,17}, //FACE_SOUTH[14]
                      {0,1,2,4,5,6,8,9,11,12,13,14,16,18}, //FACE_NORTH[14]
                       {0,1,2,3,4,5,7,8,9,10,11,14,15,18}, //FACE_BOT[14]
                        {0,1,2,3,4,6,7,8,9,10,12,13,16,17}, //FACE_TOP[14]

                   {0,2,4,5,6,8,12,14,16,18,-1,-1,-1,-1}, //EDGE_NE[10]
                    {0,1,4,5,6,9,11,13,16,18,-1,-1,-1,-1}, //EDGE_NW[10]
                     {0,1,3,5,6,7,11,13,15,17,-1,-1,-1,-1}, //EDGE_SW[10]
                      {0,2,3,5,6,10,12,14,15,17,-1,-1,-1,-1}, //EDGE_SE[10]
                       {0,1,2,4,6,8,9,12,13,16,-1,-1,-1,-1}, //EDGE_TN[10]
                        {0,1,2,3,6,7,10,12,13,17,-1,-1,-1,-1}, //EDGE_TS[10]
                         {0,1,2,4,5,8,9,11,14,18,-1,-1,-1,-1}, //EDGE_BN[10]
                          {0,1,2,3,5,7,10,11,14,15,-1,-1,-1,-1}, //EDGE_BS[10]
                           {0,2,3,4,6,8,10,12,16,17,-1,-1,-1,-1}, //EDGE_TE[10]
                            {0,1,3,4,6,7,9,13,16,17,-1,-1,-1,-1}, //EDGE_TW[10]
                             {0,2,3,4,5,8,10,14,15,18,-1,-1,-1,-1}, //EDGE_BE[10]
                              {0,1,3,4,5,7,9,11,15,18,-1,-1,-1,-1}, //EDGE_BW[10]

                   {0,2,4,6,8,12,16,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NET[7]
                    {0,1,4,6,9,13,16,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NWT[7]
                     {0,1,3,6,7,13,17,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SWT[7]
                      {0,2,3,6,10,12,17,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SET[7]
                       {0,2,4,5,8,14,18,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NEB[7]
                        {0,1,4,5,9,11,18,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NWB[7]
                         {0,1,3,5,7,11,15,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SWB[7]
                          {0,2,3,5,10,14,15,-1,-1,-1,-1,-1,-1,-1} //CORNER_SEB[7]
     }};
 #elif (C_Q == 27)
    static constexpr std::array<std::array<int, 18>,C_Q-1>
    // std::array<std::array<int, 18>,C_Q-1>
     bdry_dirs = {{{0,1,3,4,5,6,7,9,11,13,15,16,17,18,19,21,24,26}, //FACE_WEST[18]
                    {0,2,3,4,5,6,8,10,12,14,15,16,17,18,20,22,23,25}, //FACE_EAST[18] 
                     {0,1,2,3,5,6,7,10,11,12,13,14,15,17,19,22,24,25}, //FACE_SOUTH[18]
                      {0,1,2,4,5,6,8,9,11,12,13,14,16,18,20,21,23,26}, //FACE_NORTH[18]
                       {0,1,2,3,4,5,7,8,9,10,11,14,15,18,19,21,23,25}, //FACE_BOT[18]
                        {0,1,2,3,4,6,7,8,9,10,12,13,16,17,20,22,24,26}, //FACE_TOP[18]

                  {0,2,4,5,6,8,12,14,16,18,20,23,-1,-1,-1,-1,-1,-1}, //EDGE_NE[12]
                   {0,1,4,5,6,9,11,13,16,18,21,26,-1,-1,-1,-1,-1,-1}, //EDGE_NW[12]
                    {0,1,3,5,6,7,11,13,15,17,19,24,-1,-1,-1,-1,-1,-1}, //EDGE_SW[12]
                     {0,2,3,5,6,10,12,14,15,17,22,25,-1,-1,-1,-1,-1,-1}, //EDGE_SE[12]
                      {0,1,2,4,6,8,9,12,13,16,20,26,-1,-1,-1,-1,-1,-1}, //EDGE_TN[12]
                       {0,1,2,3,6,7,10,12,13,17,22,24,-1,-1,-1,-1,-1,-1}, //EDGE_TS[12]
                        {0,1,2,4,5,8,9,11,14,18,21,23,-1,-1,-1,-1,-1,-1}, //EDGE_BN[12]
                         {0,1,2,3,5,7,10,11,14,15,19,25,-1,-1,-1,-1,-1,-1}, //EDGE_BS[12]
                          {0,2,3,4,6,8,10,12,16,17,20,22,-1,-1,-1,-1,-1,-1}, //EDGE_TE[12]
                           {0,1,3,4,6,7,9,13,16,17,24,26,-1,-1,-1,-1,-1,-1}, //EDGE_TW[12]
                            {0,2,3,4,5,8,10,14,15,18,23,25,-1,-1,-1,-1,-1,-1}, //EDGE_BE[12]
                             {0,1,3,4,5,7,9,11,15,18,19,21,-1,-1,-1,-1,-1,-1}, //EDGE_BW[12]

                  {0,2,4,6,8,12,16,20,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NET[8]
                   {0,1,4,6,9,13,16,26,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NWT[8]
                    {0,1,3,6,7,13,17,24,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SWT[8]
                     {0,2,3,6,10,12,17,22,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SET[8]
                      {0,2,4,5,8,14,18,23,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NEB[8]
                       {0,1,4,5,9,11,18,21,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_NWB[8]
                        {0,1,3,5,7,11,15,19,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}, //CORNER_SWB[8]
                         {0,2,3,5,10,14,15,25,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1} //CORNER_SEB[8]
                }};  
 #endif
#endif
#endif

    // void calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& velocity);
    void calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq) { 
        LBM_Manager::calculateFeq(rho, velocity, feq); 
    }

    void calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq) { 
        LBM_Manager::calculateFeq_HeLuo(rho, velocity, feq); 
    }

    D_Phy_DDF calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, const int dir) {
        return LBM_Manager::calculateFeq(rho, velocity, dir);
    }

    D_Phy_DDF calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, const int dir) {
        return LBM_Manager::calculateFeq_HeLuo(rho, velocity, dir); 
    }
#ifndef BC_NO_GHOST
    void initial_Feq(D_mapLat* bc_ptr) {
        for (auto i_bc_lat : *bc_ptr)
        {
            D_morton lat_code = i_bc_lat.first;
            for (D_int i_q = 0; i_q < C_Q; ++i_q)
            {
                D_Phy_DDF feq_temp = calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), i_q);
                df_ptr->f[i_q].at(lat_code)  = feq_temp;
                df_ptr->fcol[i_q].at(lat_code)  = feq_temp;
            }
        }
    }
#else
    void initial_Feq(D_setLat* bc_ptr) {
        for (D_morton lat_code : *bc_ptr)
        {
            for (D_int i_q = 0; i_q < C_Q; ++i_q)
            {
                D_Phy_DDF feq_temp = calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), i_q);
                df_ptr->f[i_q].at(lat_code)  = feq_temp;
                df_ptr->fcol[i_q].at(lat_code)  = feq_temp;
            }
            
        }
    }
#endif

public:
    D3Q19_BC_Strategy() 
     : df_ptr(&(LBM_Manager::pointer_me->user[0].df)),
       density_ptr(&(LBM_Manager::pointer_me->user[0].density)),
       velocity_ptr(&(LBM_Manager::pointer_me->user[0].velocity)),
       tau_(LBM_Manager::pointer_me->tau[0]),
       CF_U_(LBM_Manager::pointer_me->CF_U) {}

    virtual void applyBCStrategy() = 0;
    virtual void initialBCStrategy() = 0;
    // void initial_BCstrategy_Ptrs(User* user) {
    // void initial_BCstrategy_Ptrs() {
    //     // df_ptr = &(user_[0].df);
    //     density_ptr = &(user_[0].density);
    //     velocity_ptr = &(user_[0].velocity);
    // };
    // void initial_BCstrategy_ddfPtrs(User* user, const int position) {
    void initial_BCstrategy_ddfPtrs(const int position) {
        if (position == 1)
            // After stream
            // f_ptr = user_[0].df.f;
            f_ptr = df_ptr->f;
        else if (position == 2)
            // After collision
            // f_ptr = user_[0].df.fcol;
            f_ptr = df_ptr->fcol;
    }
    virtual ~D3Q19_BC_Strategy() = default;
};

class D3Q19_BoundaryCondtion_Manager
{
private:
    std::map<std::string, std::unique_ptr<D3Q19_BC_Strategy>> boundaryConditionStrategies_;
    // User* user_;
    // D_Phy_real tau_;
    // D_map_define<D_Phy_Rho>* density_ptr;
    // D_map_define<D_uvw>* velocity_ptr;

public:
    // D3Q19_BoundaryCondtion_Manager(User *user, const std::array<D_Phy_real,C_max_level+1> &tau)
    //     :user_(user), tau_(tau[0]), density_ptr(&(user[0].density)), velocity_ptr(&(user[0].velocity)) {};

    template<typename strateType>
    void setup(const std::string& envelope, std::unique_ptr<strateType> strategy) {
        boundaryConditionStrategies_[envelope] = std::move(strategy);
        // boundaryConditionStrategies_[envelope]->initial_BCstrategy_Ptrs(user_);
    };

    void applyBoundaryCondition(const std::string &surface, const std::string &apply_position) {
        auto bc = boundaryConditionStrategies_.find(surface);

        if (bc != boundaryConditionStrategies_.end()) {
            if (apply_position == "AfterStream")
                // bc->second->initial_BCstrategy_ddfPtrs(user_, 1);
                bc->second->initial_BCstrategy_ddfPtrs(1);
            else if (apply_position == "AfterCollision")
                // bc->second->initial_BCstrategy_ddfPtrs(user_, 2);
                bc->second->initial_BCstrategy_ddfPtrs(2);
            else {
                std::stringstream error_log;
                error_log << "[applyBoundaryCondition] Apply " << surface << " BC at unknown position: " << apply_position << std::endl;
                error_log << "[applyBoundaryCondition] Built in BC: AfterStream, AfterCollision " << std::endl;
                log_error(error_log.str(), Log_function::logfile);
            }

            bc->second->applyBCStrategy();
        }
        else {
            std::stringstream error_log;
            error_log << "[applyBoundaryCondition]" << surface << " not contains in built-in BC: West, East, North, SOUTH, TOP, BOT" << std::endl;
            log_error(error_log.str(), Log_function::logfile);
        }
    }

    void initialBoundaryCondition(const std::string &surface) {
        auto bc = boundaryConditionStrategies_.find(surface);
        if (bc != boundaryConditionStrategies_.end()) {
            bc->second->initialBCStrategy();
        }
        else {
            std::stringstream error_log;
            error_log << "[initialBoundaryCondition]" << surface << " not contains in built-in BC: West, East, North, SOUTH, TOP, BOT" << std::endl;
            log_error(error_log.str(), Log_function::logfile);
        }
    }

    // void lbgkcollision();

    // void stream();

    virtual ~D3Q19_BoundaryCondtion_Manager() = default;
};