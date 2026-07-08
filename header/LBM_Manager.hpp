/**
 * @file LBM_Manager.h
 * @brief Define physics manager class
 * @version 0.1
 * @date 2023-03-27
 * 
 */

#pragma once

#include "Constants.h"
#include "user.h"
// #include "Collision_Manager.hpp"
// #include "Stream_Manager.hpp"
// #include "ComplexBoundary_Manager.hpp"
// #include "BC/D3Q19_BC_Manager.hpp"
#include "tecplot.h"

class Collision_Model;
class Stream_Model;
class D3Q19_BC_Strategy;
class ComplexBoundary_Model;
class D3Q19_BoundaryCondtion_Manager;

class LBM_Manager
{
friend class D3Q19_BC_Strategy;
friend class Collision_Model;
friend class SinglePointInterpolation;
// friend void LBM_Partition::calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& velocity);
// friend void IO_TECPLOTH::write(D_int iter);
// friend class nonEquilibrumBounceBack_Velocity_BC;

private:
    
    Collision_Model* collideModel; 
    Stream_Model* streamModel; 
    ComplexBoundary_Model* solidBCModel;
    D3Q19_BoundaryCondtion_Manager *bc_manager_;

    // User<D_Phy_DDF, D_Phy_Velocity, D_Phy_Rho> *user;
    User* user;
    std::vector<unsigned int> run_order;   /// order of runing simulation at different refinement level

    std::array<D_Phy_real, C_max_level+1> tau;  ///> Lattice relaxation time for levels

    // Conversion Factor for LBM
    std::array<D_Phy_real, C_max_level+1> CF_L;   ///> [m]      Conversion Factor of length for levels
    std::array<D_Phy_real, C_max_level+1> CF_T;   ///> [sec]    Conversion Factor of time for levels
    D_Phy_real CF_U;                              ///> [m/sec]  Conversion factor of velocity for levels
    std::array<D_Phy_real, C_max_level+1> CF_Nu;  ///> [m2/sec] Conversion factor of kinetic viscosoity for levels
    D_Phy_real CF_Rho; ///> [kg/m3]  Conversion factor of density for levels
#ifdef DDFTRANSFACTOR
    std::array<D_Phy_real[2], C_max_level> transFactor_C2F; ///> f(fine) = transFactor_C2F[0] * f(coarse) + transFactor_C2F[1] * feq
    std::array<D_Phy_real[2], C_max_level> transFactor_F2C; ///> f(coarse) = transFactor_F2C[0] * f(fine) + transFactor_F2C[1] * feq
#endif

    void LBMkernel(D_int refine_level = 0);

    // allocate memory
    DEPRECATED void allocMacros();
    DEPRECATED void allocDDF();
    void allocMacrosAndDDF();
    void checkMach(D_Phy_Velocity u0);
    void setUnit(std::string collision_Model);

    // initialize flow field
    void initialBC();
    void initialFeq();

    // functions for physics process
    static void calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& v);
    void calculateMacros(D_int level);
    void calculateMacros();

    static void calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq);
    static void calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, D_Phy_DDF *feq);
    static D_Phy_DDF calculateFeq(const D_Phy_Rho rho, const D_uvw velocity, const int dir);
    static D_Phy_DDF calculateFeq_HeLuo(const D_Phy_Rho rho, const D_uvw velocity, const int dir);

    void AMR_transDDF_coalescence(D_int level);
    void AMR_transDDF_explosion(D_int level);

    void defineRunOrder();

    void load_ddf(D_Phy_DDF* ddf, D_int level, D_morton code);
    void store_ddf(const D_Phy_DDF* ddf, D_int level, D_morton code);
    
    // utilities for LBM simulation
    void infoPrint(D_int iter_no);

    // friend void D3Q19_BC_Strategy::calculateFeq(D_Phy_Rho rho, D_uvw velocity, D_Phy_DDF *feq);
    // friend void D3Q19_BC_Strategy::calculateFeq_HeLuo(D_Phy_Rho rho, D_uvw velocity, D_Phy_DDF *feq);
    // friend D_Phy_DDF D3Q19_BC_Strategy::calculateFeq(D_Phy_Rho rho, D_uvw velocity, int dir);
    // friend D_Phy_DDF D3Q19_BC_Strategy::calculateFeq_HeLuo(D_Phy_Rho rho, D_uvw velocity, int dir);

    // friend void Collision_Model::calculateFeq(D_Phy_Rho rho, D_uvw velocity, D_Phy_DDF *feq);
    // friend void Collision_Model::calculateFeq_HeLuo(D_Phy_Rho rho, D_uvw velocity, D_Phy_DDF *feq);
    // friend void D3Q19_BC_Strategy::calculateMacros(const D_Phy_DDF* ddf, D_Phy_Rho& rho, D_uvw& velocity);

    friend void IO_TECPLOTH::write(D_int iter);

    std::array<uint, C_max_level+1> run_col_times{0,0,0};
    std::array<uint, C_max_level+1> run_stm_times{0,0,0};

    // IB force via momentum-exchange; Cd = 2*F_lb/(rho_lb*U_lb^2*pi*d_lb^2/4)
    D_Phy_real m_ibFx_sum{0}, m_ibFy_sum{0}, m_ibFz_sum{0};
    D_int m_ibForce_count{0};
    D_Phy_real m_cd_sum{0};
    D_int m_cd_count{0};
    void accumIBForce(D_Phy_real fx, D_Phy_real fy, D_Phy_real fz);

public:
    LBM_Manager();
    static LBM_Manager* pointer_me;
    void fluidSimulate();
    void diagL2MaxF(const char* tag);
    void diagL2RhoMin(const char* tag);
};