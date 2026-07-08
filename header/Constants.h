/**
* @file
* @brief Constants used in the solver.
* @note define constants include: 1) lattice paramters; 2) phycical parameters; 3) numerical parameters; and 4) other parameters.
*/
#pragma once

#include "General.h"
#define C_DEBUG 0
#define C_CHECK_MORTON_BOUNDARY 0
#define C_SPEED_UP 0

#define C_MAP_TYPE 1 // 1=unordered_map, 2=map, 3=MortonTrie

// Derive the backend macros consumed by Grid_Manager.h::D_map_define and
// user.h::_s_DDF from the single C_MAP_TYPE switch above, so the two cannot
// drift out of sync. settings.hpp keeps the same names commented out for the
// rare case where one wants to force a backend independently.
#if (C_MAP_TYPE == 1) && !defined(USE_HASHMAP) && !defined(USE_RBTREE) && !defined(USE_MORTONTRIE)
  #define USE_HASHMAP
#elif (C_MAP_TYPE == 2) && !defined(USE_HASHMAP) && !defined(USE_RBTREE) && !defined(USE_MORTONTRIE)
  #define USE_RBTREE
#elif (C_MAP_TYPE == 3) && !defined(USE_HASHMAP) && !defined(USE_RBTREE) && !defined(USE_MORTONTRIE)
  #define USE_MORTONTRIE
#endif
 
#define C_SEARCH_METHOD 1 // method to identify if nodes need to be removed, when set as 3, use icount_refine in a region;
                                                                          // when set as 2, use two way search and use icout_refine on sreaching boundary
                                                                          // when set as 1, use two way search
                                                                          // 1 and 2 does not impletement 2D, and does not check numerical boundary when find nodes on the searching boundary

// Grid informaion
#define C_FSI_INTERFACE 0               /// if C_SOLID_BOUNDARY == 1, using IBM, an addtional map will be used to store number of points near the node;
#define C_SOLID_BOUNDARY 1              ///< Choose method to implement solid boundary condition, 1 static boundary, 2 moving boundary;
                                        //  if C_SOLID_BOUNDARY == 2: update nodes at each finest time step

namespace airplane
{
    const std::string F_model_path = "/home/AMR-framework/stl/airplane.txt";
    const std::string OUTPUT_NAME = "airplane";

    // Region settings
    constexpr D_real C_xb = 20000;  ///< Background boundary distance in x direction
    constexpr D_real C_yb = 20000;  ///< Background boundary distance in y direction
    constexpr D_real C_zb = 20000;  ///< Background boundary distance in z direction

    // Legacy domain boundary variables for compatibility
    constexpr D_real xb_domain = C_xb;
    constexpr D_real yb_domain = C_yb;
    constexpr D_real zb_domain = C_zb;

    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {0, 0, 0,
                               C_xb, C_yb, C_zb};

    // Solid model origin (absolute coordinate). STL local (0,0,0) is placed
    // here. Edit these three values to manually control where the solid sits
    // inside the domain (e.g., sphere center for sphere-flow test).
    constexpr D_real C_solid_origin[3] = {(C_domain[0] + C_domain[3]) / 2,
                                          (C_domain[1] + C_domain[4]) / 2,
                                          (C_domain[2] + C_domain[5]) / 2};

    // Mesh settings
    const D_real C_dx = 256;    ///< Grid space of the background mesh
    const int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    const int C_y0b_offset = 0;
    const int C_z0b_offset = 0;
    const unsigned int C_max_level = 4;    ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
    const unsigned int C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)
    const unsigned int C_extend_inner = 4;       ///< Number of nodes extended from solid points in the finest block
    const unsigned int C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    const unsigned int C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    const unsigned int C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    const unsigned int C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    const unsigned int C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    const unsigned int C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)
    const unsigned int C_extend = 4;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)
    const unsigned int C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    const unsigned int C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    const unsigned int C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    const unsigned int C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    const unsigned int C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    const unsigned int C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)
    const unsigned int C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid

    // LBM physical parameters
    static D_real U_dt = 0.001f;               ///< [s]
    constexpr D_real U_u0 = 1.0;               ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;                ///< initial velocity in y direction
    constexpr D_real U_w0 = 0.;                ///< initial velocity in z direction
    constexpr D_real U_rho0 = 1.;              ///< initial density
    constexpr D_real U_kineViscosity = 0.002;  ///< kinetic viscosity
    constexpr uint U_iters = 100;              ///< number of LBM iterations

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;

}

namespace hm
{
    const std::string F_model_path = "/home/AMR-framework/stl/ln_j20.stl";

    #define C_DIMS 3                      ///< Number of dimensions
    #define C_Q 19                        ///< Number of discrete velocities (global: last #define wins; all set to 19 for D3Q19 sphere validation)

    // Region settings - Unified domain boundary definition
    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {-10, -50, -60,
                               400, 150, 60};

    // Solid model origin (absolute coordinate). STL local (0,0,0) is placed
    // here. Edit to manually control where the solid sits inside the domain.
    constexpr D_real C_solid_origin[3] = {(C_domain[0] + C_domain[3]) / 2,
                                          (C_domain[1] + C_domain[4]) / 2,
                                          (C_domain[2] + C_domain[5]) / 2};

    // Mesh settings
    constexpr D_real C_dx = 20;      ///< Grid space of the background mesh
    constexpr int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    constexpr int C_y0b_offset = 0;
    constexpr int C_z0b_offset = 0;
    constexpr uint C_max_level = 4;    ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
    constexpr uint C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)
    constexpr uint C_extend_inner = 4;       ///< Number of nodes extended from solid points in the finest block
    constexpr uint C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    constexpr uint C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    constexpr uint C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    constexpr uint C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    constexpr uint C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    constexpr uint C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)
    constexpr uint C_extend = 4;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)
    constexpr uint C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    constexpr uint C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    constexpr uint C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    constexpr uint C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    constexpr uint C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    constexpr uint C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)
    constexpr uint C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid

    // Simulation setting
    constexpr uint U_iters = 100000;
    constexpr uint U_outputIters = 1000;

    // Numerical parameters
    static D_real U_dt = 2.0;               ///< [s]

    // Initial physical variables
    constexpr D_real U_u0 = 10.0;  ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_w0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_rho0 = 1.;  ///< initial density
    constexpr D_real U_kineViscosity = 0.00001;	   ///< kinetic viscosoity
    constexpr D_real U_dynaViscosity = U_kineViscosity * U_rho0;

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;
}  

namespace sphere
{
    const std::string F_model_path = "stl/sphere.stl";
    const std::string OUTPUT_NAME = "sphere";

    #define C_DIMS 3                      ///< Number of dimensions
    #define C_Q 19                        ///< D3Q19

    // Sphere flow benchmark: Re = U*d/nu = 100, target Cd ≈ 1.09
    // sphere.stl center at STL (10,10,10), placed at domain (4,4,4)

    // Region settings
    constexpr D_real C_xb = 12;  ///< Background boundary distance in x direction (domain length, streamwise)
    constexpr D_real C_yb = 8;   ///< Background boundary distance in y direction
    constexpr D_real C_zb = 8;   ///< Background boundary distance in z direction

    // Legacy domain boundary variables for compatibility (absolute upper bound, since C_domain starts at 0)
    constexpr D_real xb_domain = C_xb;
    constexpr D_real yb_domain = C_yb;
    constexpr D_real zb_domain = C_zb;

    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {0, 0, 0,
                               C_xb, C_yb, C_zb};

    // Solid model origin: STL (10,10,10) → domain (4,4,4)
    constexpr D_real C_solid_origin[3] = {-6, -6, -6};

    // Mesh settings
    constexpr D_real C_dx = 0.1;      ///< background grid spacing
    constexpr int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    constexpr int C_y0b_offset = 0;
    constexpr int C_z0b_offset = 0;
    constexpr uint C_max_level = 2;
    constexpr bool C_IB_use_FH = true;
    constexpr uint C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)
    constexpr uint C_extend_inner = 2;       ///< Number of nodes extended from solid points in the finest block
    constexpr uint C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    constexpr uint C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    constexpr uint C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    constexpr uint C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    constexpr uint C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    constexpr uint C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)
    constexpr uint C_extend = 4;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)
    constexpr uint C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    constexpr uint C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    constexpr uint C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    constexpr uint C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    constexpr uint C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    constexpr uint C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)
    constexpr uint C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid

    // Simulation setting
    constexpr uint U_iters = 1200;
    constexpr uint U_outputIters = 1000;

    // Numerical parameters
    static D_real U_dt = 0.0058f;

    // Initial physical variables
    constexpr D_real U_u0 = 1.0;   ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;    ///< initial velocity in y direction
    constexpr D_real U_w0 = 0.;    ///< initial velocity in z direction
    constexpr D_real U_rho0 = 1.;  ///< initial density
    constexpr D_real U_kineViscosity = 0.01;	   ///< Re = 100
    constexpr D_real U_dynaViscosity = U_kineViscosity * U_rho0;

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;
}

namespace dji
{
    const std::string F_model_path = "/home/lb-amr/stl/dji_air3.stl";
    const std::string OUTPUT_NAME = "dji";


    #define C_DIMS 3                      ///< Number of dimensions
    #define C_Q 19                        ///< Number of discrete velocities (global: last #define wins; all set to 19 for D3Q19 sphere validation)



    // Region settings - Unified domain boundary definition
    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {-100, -100, -100,
                               300, 300, 300};

    // Solid model origin (absolute coordinate). STL local (0,0,0) is placed
    // here. Edit to manually control where the solid sits inside the domain.
    constexpr D_real C_solid_origin[3] = {(C_domain[0] + C_domain[3]) / 2,
                                          (C_domain[1] + C_domain[4]) / 2,
                                          (C_domain[2] + C_domain[5]) / 2};

    // Mesh settings
    constexpr D_real C_dx = 0.25;      ///< Grid space of the background mesh

    constexpr int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    constexpr int C_y0b_offset = 0;
    constexpr int C_z0b_offset = 0;

    constexpr uint C_max_level = 3;    ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
    const unsigned int C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)

    const unsigned int C_extend_inner = 2;       ///< Number of nodes extended from solid points in the finest block
    const unsigned int C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    const unsigned int C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    const unsigned int C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    const unsigned int C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    const unsigned int C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    const unsigned int C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)

    const unsigned int C_extend = 2;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)

    const unsigned int C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    const unsigned int C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    const unsigned int C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    const unsigned int C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    const unsigned int C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    const unsigned int C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)

    const unsigned int C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid



    // Simulation setting
    constexpr uint U_iters = 100;
    constexpr uint U_outputIters = 100;

    // Numerical parameters
    static D_real U_dt = 0.001f;   ///< [s]

    // Initial physical variables
    constexpr D_real U_u0 = 1.0;   ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_w0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_rho0 = 1.;  ///< initial density
    constexpr D_real U_kineViscosity = 0.002;	   ///< kinetic viscosoity
    constexpr D_real U_dynaViscosity = U_kineViscosity * U_rho0;

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;

    
}

namespace paris
{
    const std::string F_model_path = "/home/lb-amr/stl/EiffelTowerTALL.stl";
    const std::string OUTPUT_NAME = "paris";


    #define C_DIMS 3                      ///< Number of dimensions
    #define C_Q 19                        ///< Number of discrete velocities (global: last #define wins; all set to 19 for D3Q19 sphere validation)



    // Region settings - Unified domain boundary definition
    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {-200, -200, -100,
                               600, 600, 900};

    // Solid model origin (absolute coordinate). STL local (0,0,0) is placed
    // here. Edit to manually control where the solid sits inside the domain.
    constexpr D_real C_solid_origin[3] = {(C_domain[0] + C_domain[3]) / 2,
                                          (C_domain[1] + C_domain[4]) / 2,
                                          (C_domain[2] + C_domain[5]) / 2};

    // Mesh settings
    constexpr D_real C_dx = 128;      ///< Grid space of the background mesh

    constexpr int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    constexpr int C_y0b_offset = 0;
    constexpr int C_z0b_offset = 0;

    constexpr uint C_max_level = 3;    ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
    const unsigned int C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)

    const unsigned int C_extend_inner = 2;       ///< Number of nodes extended from solid points in the finest block
    const unsigned int C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    const unsigned int C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    const unsigned int C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    const unsigned int C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    const unsigned int C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    const unsigned int C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)

    const unsigned int C_extend = 2;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)

    const unsigned int C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    const unsigned int C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    const unsigned int C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    const unsigned int C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    const unsigned int C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    const unsigned int C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)

    const unsigned int C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid



    // Simulation setting
    constexpr uint U_iters = 100;
    constexpr uint U_outputIters = 100;

    // Numerical parameters
    static D_real U_dt = 0.001f;   ///< [s]

    // Initial physical variables
    constexpr D_real U_u0 = 1.0;   ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_w0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_rho0 = 1.;  ///< initial density
    constexpr D_real U_kineViscosity = 0.002;	   ///< kinetic viscosoity
    constexpr D_real U_dynaViscosity = U_kineViscosity * U_rho0;

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;

    
}

namespace bunny
{
    const std::string F_model_path = "/home/lb-amr/stl/stanford_bunny.stl";
    const std::string OUTPUT_NAME = "bunny";


    #define C_DIMS 3                      ///< Number of dimensions
    #define C_Q 19                        ///< Number of discrete velocities (global: last #define wins; all set to 19 for D3Q19 sphere validation)



    // Region settings - Unified domain boundary definition
    // Domain boundaries: {x_min, y_min, z_min, x_max, y_max, z_max}
    constexpr D_real C_domain[6] = {-400, -400, -400,
                               1200, 1200, 1200};

    // Solid model origin (absolute coordinate). STL local (0,0,0) is placed
    // here. Edit to manually control where the solid sits inside the domain.
    constexpr D_real C_solid_origin[3] = {(C_domain[0] + C_domain[3]) / 2,
                                          (C_domain[1] + C_domain[4]) / 2,
                                          (C_domain[2] + C_domain[5]) / 2};

    // Mesh settings
    constexpr D_real C_dx = 5.5;      ///< Grid space of the background mesh

    constexpr int C_x0b_offset = 0;    ///< offset to avoid Morton_assit::find_x0 exceeding the boundary limit, the offset distance is (x0b_offset * C_dx)
    constexpr int C_y0b_offset = 0;
    constexpr int C_z0b_offset = 0;

    constexpr uint C_max_level = 3;    ///< Maximum refinement level, the background mesh is level 0, others are 1, 2, 3, ..., C_maxlevel
    const unsigned int C_overlap = 1;         ///< Number of overlapping grids (corresponding to coarser grid at ilevel - 1)

    const unsigned int C_extend_inner = 2;       ///< Number of nodes extended from solid points in the finest block
    const unsigned int C_extend_inner_x0 = 2;      ///< additional extension in -x direction (inner block)
    const unsigned int C_extend_inner_x1 = 2;      ///< additional extension in +x direction (inner block)
    const unsigned int C_extend_inner_y0 = 2;      ///< additional extension in -y direction (inner block)
    const unsigned int C_extend_inner_y1 = 2;      ///< additional extension in +y direction (inner block)
    const unsigned int C_extend_inner_z0 = 2;      ///< additional extension in -z direction (inner block)
    const unsigned int C_extend_inner_z1 = 2;      ///< additional extension in +z direction (inner block)

    const unsigned int C_extend = 2;             ///< Number of nodes extended from the inner block to outer blocks (iblock < C_max_level)

    const unsigned int C_extend_outer_x0 = 2;      ///< additional extension in -x direction (outer block)
    const unsigned int C_extend_outer_x1 = 2;      ///< additional extension in +x direction (outer block)
    const unsigned int C_extend_outer_y0 = 2;      ///< additional extension in -y direction (outer block)
    const unsigned int C_extend_outer_y1 = 2;      ///< additional extension in +y direction (outer block)
    const unsigned int C_extend_outer_z0 = 2;      ///< additional extension in -z direction (outer block)
    const unsigned int C_extend_outer_z1 = 2;      ///< additional extension in +z direction (outer block)

    const unsigned int C_extend_ghost = 0;      ///< Number of nodes extended from solid points inside the solid



    // Simulation setting
    constexpr uint U_iters = 100;
    constexpr uint U_outputIters = 100;

    // Numerical parameters
    static D_real U_dt = 0.001f;   ///< [s]

    // Initial physical variables
    constexpr D_real U_u0 = 1.0;   ///< initial velocity in x direction
    constexpr D_real U_v0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_w0 = 0.;    ///< initial velocity in x direction
    constexpr D_real U_rho0 = 1.;  ///< initial density
    constexpr D_real U_kineViscosity = 0.002;	   ///< kinetic viscosoity
    constexpr D_real U_dynaViscosity = U_kineViscosity * U_rho0;

    // Reference physical variables
    constexpr D_real R_velocity = U_u0;
    constexpr D_real R_length = 1.;
    constexpr D_real R_Re = R_velocity * R_length / U_kineViscosity;



}


using namespace sphere;

#if (C_FSI_INTERFACE == 1)
const unsigned int C_extend_IB = 1;      ///< Number of nodes extended from solid points for IBM
#endif

// Physical Constraint for stability
constexpr D_real MAX_MACH = 0.09;           ///< Max Mach number of problems that LB can solve
constexpr D_real MIN_TAU = 0.500001;        ///< Min lattice tau, for stability of LB simulation

#if (C_DIMS == 3)
 #if (C_Q == 19)
/* D3Q19             
                         ___________________              ________15_________  
                        /:                 /|            /:                 /|
                       / :      5         / |          14 :               11 |
  z                   /__:______|_3_____ /  |         /________18________ /  |
 /|\ y               |   :      |/      |   |         |  10              |   7 
  | /                | 2-:------0-------|-1 |         |   :      .       |   |
  |/                 |   :     /|       |   |         |   :              |   |
  o------->x         |   :----/-|-------|---|         8   :-------17-----9---|
                     |  /    4  |       |  /          |  /               |   /
                     | /        6       | /           | 12               | 13  
                     |/_________________|/            |/_______16________|/    
*/
constexpr int ex[C_Q] = {
    0,  
    1, -1,  0,  0,  0,  0,  
    1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0};
constexpr int ey[C_Q] = {
    0,  
    0,  0,  1, -1,  0,  0,  
    1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1};
constexpr int ez[C_Q] = {
    0,  
    0,  0,  0,  0,  1, -1,  
    0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1};

constexpr unsigned int e_inv[C_Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17};

 #elif (C_Q == 27)
/* D3Q27             
                             ___________________              ________15_________              25________________19 
                            /:                 /|            /:                 /|            /:                 /| 
                           / :      5         / |          14 :               11 |           / :                / | 
      z                   /__:______|_3_____ /  |         /________18________ /  |         23__:______________21  | 
     /|\ y               |   :      |/      |   |         |  10              |   7         |   :              |   | 
      | /                | 2-:------0-------|-1 |         |   :      .       |   |         |   :              |   | 
      |/                 |   :     /|       |   |         |   :              |   |         |   :              |   | 
      o------->x         |   :----/-|-------|---|         8   :-------17-----9---|         |   22-------------|--24 
                         |  /    4  |       |  /          |  /               |   /         |  /               |  /  
                         | /        6       | /           | 12               | 13          | /                | /   
                         |/_________________|/            |/_______16________|/            20_________________26    
*/
           
constexpr int ex[C_Q] = {
    0, 
    1, -1,  0,  0,  0,  0, 
    1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0, 
    1, -1,  1, -1, -1,  1, -1,  1};
constexpr int ey[C_Q] = {
    0, 
    0,  0,  1, -1, 0, 0, 
    1, -1, -1,  1,  0,  0,  0,  0,  1, -1,  1, -1,
    1, -1, -1,  1, -1,  1,  1, -1};
constexpr int ez[C_Q] = {
    0,  
    0,  0,  0,  0,  1, -1,  
    0,  0,  0,  0,  1, -1, -1,  1,  1, -1, -1,  1,
    1, -1,  1, -1,  1, -1,  1, -1};

constexpr unsigned int e_inv[C_Q] = {0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17, 20, 19, 22, 21, 24, 23, 26, 25};

 #endif
#endif