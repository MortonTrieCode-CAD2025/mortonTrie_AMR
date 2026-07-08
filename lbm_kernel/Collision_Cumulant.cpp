/**
 * @file Collision_Model.cpp
 * @brief the implementation of the cumulant collsion
 * @version 0.1
 * @date 2024-02-21
 * 
 */

#include "lbm_kernel/Collision_Cumulant.hpp"
#include "Lat_Manager.hpp"
#include "lbm_kernel/les.hpp"

// rho floor: prevent 1/rho & 1/rho^2 super-exponential blow-up when rho->0.
// bwd_cum_trans _b ~ u^6/rho^2 amplifies a local rho->0 into f×1e30/step (the NaN amplifier
// identified 2026-06-24). Clamp |rho|>=RHO_FLOOR so cumulants stay bounded; the root cause
// (who pushes rho->0) then surfaces as a bounded rho-anomaly instead of an opaque NaN.
#ifndef RHO_FLOOR
#define RHO_FLOOR 1e-6
#endif
#define SAFE_IRHO(r) (1.0 / (((r) > -RHO_FLOOR && (r) < RHO_FLOOR) ? (((r) < 0.0) ? -RHO_FLOOR : RHO_FLOOR) : (r)))

void Collision_Cumulant::collide(D_int i_level)
{
#if (C_Q != 27)
    log_error("The collision model is set as CUMULANT model. The discreted model is resetted to D3Q27.", Log_function::logfile);
#endif

    // D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(i_level));
    // D_Phy_real tau_alias = tau_[i_level];

    // auto collision_function = [i_level, &tau_alias] (D_mapLat* lattice_ptr, User* user, _s_DDF* f_ptr)
    auto collision_function = [i_level, this] (D_mapLat* lattice_ptr)
    {
        _s_DDF* f_ptr = &(user_[i_level].df);
        for (auto lat_iter : *lattice_ptr)
        {
            D_morton lat_code = lat_iter.first;
            // D_morton xxx_code("0000000000000000000000000000000000001010101001001111100001000010");
            
            D_Phy_real tau0 = tau_[i_level];
            #ifdef LES
                D_Phy_real stressNorm = Turbulent_LES::calculate_stress_tensor(ddf, feq);
                tau0 = Turbulent_LES::smag_subgrid_tau(omega, stressNorm);
            #endif

            D_Phy_DDF ddf[27];
#if (C_MAP_TYPE == 3)
            // Hot path: 1 trie walk + 27 contiguous double loads instead of 27
            // independent `f[i_q].at(code)` calls (each a full root-to-leaf walk).
            f_ptr->f[0].parent()->load_all(lat_code, refine_level, ddf);
#else
            for (D_int i_q = 0; i_q < 27; ++i_q) {
                ddf[i_q] = f_ptr->f[i_q].at(lat_code);
            }
#endif
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i_q = 0; i_q < C_Q; ++i_q) {
            //         std::cout << "  f[" << i_q << "] = " << ddf[i_q];
            //     }
            //     std::cout << std::endl;
            // }

            // declare the local array
            D_Phy_DDF k[27] = {0.}, C[27] = {0.};

            // Cache velocity/density for this cell: collision never updates them,
            // and caching collapses 5 root-to-leaf trie walks into 1 each.
            const D_uvw& vel = user_[i_level].velocity.at(lat_code);
            const D_Phy_Rho& rho = user_[i_level].density.at(lat_code);

            // initialize the k
            fwd_central_moment_trans(ddf, vel, k);
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i = 0; i < C_Q; ++i) {
            //         printf(" fwd_central_moment_trans u=%f v=%f w=%f k[%d] = %f ddf[%d] = %f\n",user_[i_level].velocity.at(lat_code).x,user_[i_level].velocity.at(lat_code).y,user_[i_level].velocity.at(lat_code).z, i, k[i], i, ddf[i]);
            //     }
            // }
            // for (int i = 0; i < 27; ++i) {
            //     printf(" fwd_central_moment_trans u=%f v=%f w=%f k[%d] = %f ddf[%d] = %f\n",user_[i_level].velocity.at(lat_code).x,user_[i_level].velocity.at(lat_code).y,user_[i_level].velocity.at(lat_code).z, i, k[i], i, ddf[i]);
            // }

            // forward cumulants transformation
            fwd_cum_trans(k, rho, C);
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i = 0; i < C_Q; ++i) {
            //         printf(" fwd_cum_trans k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            //     }
            // }
            // for (int i = 0; i < 27; ++i) {
            //     printf(" fwd_cum_trans k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            // }

            // cumulant collide kernel
            cum_collide_kernel(C, k, tau0, vel, rho);
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i = 0; i < C_Q; ++i) {
            //         printf(" cum_collide_kernel k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            //     }
            // }
            // for (int i = 0; i < 27; ++i) {
            //     printf(" cum_collide_kernel k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            // }

            // backward cumulants transformation
            bwd_cum_trans(C, k, rho);
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i = 0; i < C_Q; ++i) {
            //         printf(" bwd_cum_trans k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            //     }
            // }
            // for (int i = 0; i < 27; ++i) {
            //     printf(" bwd_cum_trans k[%d] = %f C[%d] = %f\n", i, k[i], i, C[i]);
            // }

            // store k to ddf
            bwd_central_moment_trans(k, ddf, vel);
            // if (i_level==4 && xxx_code==lat_code ) {
            //     for (D_int i = 0; i < C_Q; ++i) {
            //         printf(" bwd_central_moment_trans k[%d] = %f ddf[%d] = %f\n", i, k[i], i, ddf[i]);
            //     }
            // }
            // for (int i = 0; i < 27; ++i) {
            //     printf(" bwd_central_moment_trans k[%d] = %f ddf[%d] = %f\n", i, k[i], i, ddf[i]);
            // }

#if (C_MAP_TYPE == 3)
            // Hot path: 1 trie walk + 27 contiguous double stores instead of 27
            // independent `fcol[i_q].at(code) = ...` calls (each a full
            // root-to-leaf walk).
            f_ptr->fcol[0].parent()->store_all(lat_code, refine_level, ddf);
#else
            for (D_int i_q = 0; i_q < 27; ++i_q) {
                f_ptr->fcol[i_q].at(lat_code) = ddf[i_q];
            }
#endif
        }
    };
    
    collision_function(&(Lat_Manager::pointer_me->lat_f.at(i_level)));

    collision_function(&(Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level)));
    
    for (auto lat_iter : Lat_Manager::pointer_me->lat_overlap_C2F.at(i_level))
    {
        D_morton lat_code = lat_iter.first;
#if (C_MAP_TYPE == 3)
        // Bulk copy f -> fcol for overlap cells: 2 trie walks (1 load + 1 store)
        // instead of 54 (27 reads on f + 27 writes on fcol).
        D_Phy_DDF tmp[27];
        user_[i_level].df.f[0].parent()->load_all(lat_code, refine_level, tmp);
        user_[i_level].df.fcol[0].parent()->store_all(lat_code, refine_level, tmp);
#else
        for (D_int i_q = 0; i_q < 27; ++i_q)
        {
            user_[i_level].df.fcol[i_q].at(lat_code) = user_[i_level].df.f[i_q].at(lat_code);
        }
#endif
    }
}

// decode the array to variable
#define k000  (k[0])

#define k100  (k[1])
#define k010  (k[2])
#define k001  (k[3])

#define k200  (k[4])
#define k020  (k[5])
#define k002  (k[6])
#define k110  (k[7])
#define k101  (k[8])
#define k011  (k[9])

#define k210  (k[10])
#define k201  (k[11])
#define k021  (k[12])
#define k120  (k[13])
#define k102  (k[14])
#define k012  (k[15])
#define k111  (k[16])

#define k220  (k[17])
#define k202  (k[18])
#define k022  (k[19])
#define k211  (k[20])
#define k121  (k[21])
#define k112  (k[22])

#define k221  (k[23])
#define k212  (k[24])
#define k122  (k[25])

#define k222  (k[26])

#define C000  (C[0])

#define C100  (C[1])
#define C010  (C[2])
#define C001  (C[3])

#define C200  (C[4])
#define C020  (C[5])
#define C002  (C[6])
#define C110  (C[7])
#define C101  (C[8])
#define C011  (C[9])

#define C210  (C[10])
#define C201  (C[11])
#define C021  (C[12])
#define C120  (C[13])
#define C102  (C[14])
#define C012  (C[15])
#define C111  (C[16])

#define C220  (C[17])
#define C202  (C[18])
#define C022  (C[19])
#define C211  (C[20])
#define C121  (C[21])
#define C112  (C[22])

#define C221  (C[23])
#define C212  (C[24])
#define C122  (C[25])

#define C222  (C[26])

// macro for backward moment transformation
#define kN00 (k1[0 ])
#define kN01 (k1[1 ])
#define kN02 (k1[2 ])
#define kN10 (k1[3 ])
#define kN11 (k1[4 ])
#define kN12 (k1[5 ])
#define kN20 (k1[6 ])
#define kN21 (k1[7 ])
#define kN22 (k1[8 ])
#define kZ00 (k1[9 ])
#define kZ01 (k1[10])
#define kZ02 (k1[11])
#define kZ10 (k1[12])
#define kZ11 (k1[13])
#define kZ12 (k1[14])
#define kZ20 (k1[15])
#define kZ21 (k1[16])
#define kZ22 (k1[17])
#define kP00 (k1[18])
#define kP01 (k1[19])
#define kP02 (k1[20])
#define kP10 (k1[21])
#define kP11 (k1[22])
#define kP12 (k1[23])
#define kP20 (k1[24])
#define kP21 (k1[25])
#define kP22 (k1[26])

#define kNN0 (k2[0 ])
#define kNN1 (k2[1 ])
#define kNN2 (k2[2 ])
#define kNZ0 (k2[3 ])
#define kNZ1 (k2[4 ])
#define kNZ2 (k2[5 ])
#define kNP0 (k2[6 ])
#define kNP1 (k2[7 ])
#define kNP2 (k2[8 ])
#define kZN0 (k2[9 ])
#define kZN1 (k2[10])
#define kZN2 (k2[11])
#define kZZ0 (k2[12])
#define kZZ1 (k2[13])
#define kZZ2 (k2[14])
#define kZP0 (k2[15])
#define kZP1 (k2[16])
#define kZP2 (k2[17])
#define kPN0 (k2[18])
#define kPN1 (k2[19])
#define kPN2 (k2[20])
#define kPZ0 (k2[21])
#define kPZ1 (k2[22])
#define kPZ2 (k2[23])
#define kPP0 (k2[24])
#define kPP1 (k2[25])
#define kPP2 (k2[26])

#define kPZZ (k3[0 ])
#define kNZZ (k3[1 ])
#define kZPZ (k3[2 ])
#define kZNZ (k3[3 ])
#define kZZP (k3[4 ])
#define kZZN (k3[5 ])
#define kZPP (k3[6 ])
#define kZNN (k3[7 ])
#define kZPN (k3[8 ])
#define kZNP (k3[9 ])
#define kPZP (k3[10])
#define kNZN (k3[11])
#define kPZN (k3[12])
#define kNZP (k3[13])
#define kPPZ (k3[14])
#define kNNZ (k3[15])
#define kPNZ (k3[16])
#define kNPZ (k3[17])
#define kPPP (k3[18])
#define kNNN (k3[19])
#define kNNP (k3[20])
#define kPPN (k3[21])
#define kNPP (k3[22])
#define kPNN (k3[23])
#define kPNP (k3[24])
#define kNPN (k3[25])
#define kZZZ (k3[26])

/**
 * @brief step 1. initialize the k
 * @param[in] ddf
 * @param[out] k
*/
void Collision_Cumulant::fwd_central_moment_trans(const D_Phy_DDF* ddf, const D_uvw velocity, D_Phy_DDF* k)
{
    D_Phy_DDF cMux, cMuy, cMuz, f;
    for (int q = 0; q < 27; ++q) 
    {
        cMux = ex[q] - velocity.x;
        cMuy = ey[q] - velocity.y;
        cMuz = ez[q] - velocity.z;
        f = ddf[q];

        // order 0
        k000 += f;

        // order 1
        k100 += cMux * f;
        k010 += cMuy * f;
        k001 += cMuz * f;

        // Order 2
        k200 += cMux * cMux * f;
        k020 += cMuy * cMuy * f;
        k002 += cMuz * cMuz * f;
        k110 += cMux * cMuy * f;
        k101 += cMux * cMuz * f;
        k011 += cMuy * cMuz * f;
        // Order 3
        k210 += cMux * cMux * cMuy * f;
        k201 += cMux * cMux * cMuz * f;
        k021 += cMuy * cMuy * cMuz * f;
        k120 += cMux * cMuy * cMuy * f;
        k102 += cMux * cMuz * cMuz * f;
        k012 += cMuy * cMuz * cMuz * f;
        k111 += cMux * cMuy * cMuz * f;
        // Order 4
        k220 += cMux * cMux * cMuy * cMuy * f;
        k202 += cMux * cMux * cMuz * cMuz * f;
        k022 += cMuy * cMuy * cMuz * cMuz * f;
        k211 += cMux * cMux * cMuy * cMuz * f;
        k121 += cMux * cMuy * cMuy * cMuz * f;
        k112 += cMux * cMuy * cMuz * cMuz * f;
        // Order 5
        k221 += cMux * cMux * cMuy * cMuy * cMuz * f;
        k212 += cMux * cMux * cMuy * cMuz * cMuz * f;
        k122 += cMux * cMuy * cMuy * cMuz * cMuz * f;
        // Order 6
        k222 += cMux * cMux * cMuy * cMuy * cMuz * cMuz * f;
    }
}

/**
 * @brief step 2. forward cumulants transformation
 * @param[in] k
 * @param[out] C
*/
void Collision_Cumulant::fwd_cum_trans(const D_Phy_DDF *k, const D_Phy_Rho rho, D_Phy_DDF *C) 
{
    D_Phy_Rho irho = SAFE_IRHO(rho);

    C110 = k110, C011 = k011, C101 = k101;
    C200 = k200, C020 = k020, C002 = k002;
    C120 = k120, C012 = k012, C201 = k201;
    C111 = k111;

    C211 = k211 - (k200 * k011 + 2 * k110 * k101) * irho;
    C121 = k121 - (k020 * k101 + 2 * k011 * k110) * irho;
    C112 = k112 - (k002 * k110 + 2 * k101 * k011) * irho;

    C220 = k220 - (k200 * k020 + 2 * k110 * k110) * irho;
    C022 = k022 - (k020 * k002 + 2 * k011 * k011) * irho;
    C202 = k202 - (k002 * k200 + 2 * k101 * k101) * irho;

    C122 = k122 - (k002 * k120 + k020 * k102 + 4 * k011 * k111 + 2 * (k101 * k021 + k110 * k012)) * irho;
    C212 = k212 - (k200 * k012 + k002 * k210 + 4 * k101 * k111 + 2 * (k110 * k102 + k011 * k201)) * irho;
    C221 = k221 - (k020 * k201 + k200 * k021 + 4 * k110 * k111 + 2 * (k011 * k210 + k101 * k120)) * irho;

    D_Phy_real _a = (4 * k111 * k111 + k200 * k022 + k020 * k202 + k002 * k220 + 4 * (k011 * k211 + k101 * k121 + k110 * k112) + 2 * (k120 * k102 + k210 * k012 + k201 * k021)) * irho;
    D_Phy_real _b = (16 * k110 * k101 * k011 + 4 * (k101 * k101 * k020 + k011 * k011 * k200 + k110 * k110 * k002) + 2 * k200 * k020 * k002) * irho * irho;
    C222 = k222 - _a + _b;
}

/**
 * @brief step 3. cumulant collide kernel
 * @param[in] k
 * @param[out] C
*/
void Collision_Cumulant::cum_collide_kernel(D_Phy_DDF *C, const D_Phy_DDF *k, const D_Phy_real tau, const D_uvw velocity, const D_Phy_Rho rho) 
{
    // prepare some variables
    // D_Phy_Velocity _u = velocity.x, _v = velocity.y, _w = velocity.z;
    D_Phy_DDF _omega1 = 1./tau;
    D_Phy_DDF _omega2 = 1.;

    D_Phy_Rho irho = SAFE_IRHO(rho);
    D_Phy_DDF X1, X2, X3, Y1, Y2, Y3; // tmp

    // derivation
    D_Phy_DDF Dx_u = -(_omega1 * (2 * C200 - C020 - C002) + _omega2 * (C200 + C020 + C002 - k000)) * irho * 0.5;
    D_Phy_DDF Dy_v = Dx_u + 3 * _omega1 * (C200 - C020) * 0.5 * irho;
    D_Phy_DDF Dz_w = Dx_u + 3 * _omega1 * (C200 - C002) * 0.5 * irho;

    // 2nd order
    C110 *= (1 - _omega1);
    C101 *= (1 - _omega1);
    C011 *= (1 - _omega1);
    X1 = (1 - _omega1) * (C200 - C020) - 3 * rho * (1 - _omega1 * 0.5) * (SQ(velocity.x) * Dx_u - SQ(velocity.y) * Dy_v);  // C200* - C020*
    X2 = (1 - _omega1) * (C200 - C002) - 3 * rho * (1 - _omega1 * 0.5) * (SQ(velocity.x) * Dx_u - SQ(velocity.z) * Dz_w);  // C200* - C002*
    X3 = k000 * _omega2 + (1 - _omega2) * (C200 + C020 + C002) - 3 * rho * (1 - _omega2 * 0.5) * (SQ(velocity.x) * Dx_u + SQ(velocity.y) * Dy_v + SQ(velocity.z) * Dz_w);  // C200* + C020* + C002*
    C200 = (X1 + X2 + X3) / 3;
    C020 = C200 - X1;
    C002 = C200 - X2;

    C120 = 0, C102 = 0, C210 = 0;
    C012 = 0, C201 = 0, C021 = 0;
    C111 = 0;
    C220 = 0, C202 = 0, C022 = 0;
    C211 = 0, C121 = 0, C112 = 0;
    C221 = 0, C212 = 0, C122 = 0;
    C222 = 0;

}

/**
 * @brief step 4. backward cumulants transformation
 * @param[in] C
 * @param[out] k
*/
void Collision_Cumulant::bwd_cum_trans(const D_Phy_DDF *C, D_Phy_DDF *k, D_Phy_Rho rho) 
{
    D_Phy_Rho irho = SAFE_IRHO(rho);
    D_Phy_Rho drho = rho - 1.0;

    k110 = C110, k011 = C011, k101 = C101;
    k200 = C200, k020 = C020, k002 = C002;
    k120 = C120, k012 = C012, k201 = C201;
    k111 = C111;

    k211 = C211 + (k200 * k011 + 2 * k110 * k101) * irho;
    k121 = C121 + (k020 * k101 + 2 * k011 * k110) * irho;
    k112 = C112 + (k002 * k110 + 2 * k101 * k011) * irho;

    k220 = C220 + (k200 * k020 + 2 * k110 * k110) * irho;
    k022 = C022 + (k020 * k002 + 2 * k011 * k011) * irho;
    k202 = C202 + (k002 * k200 + 2 * k101 * k101) * irho;

    k122 = C122 + (k002 * k120 + k020 * k102 + 4 * k011 * k111 + 2 * (k101 * k021 + k110 * k012)) * irho;
    k212 = C212 + (k200 * k012 + k002 * k210 + 4 * k101 * k111 + 2 * (k110 * k102 + k011 * k201)) * irho;
    k221 = C221 + (k020 * k201 + k200 * k021 + 4 * k110 * k111 + 2 * (k011 * k210 + k101 * k120)) * irho;

    D_Phy_DDF _a = (4 * k111 * k111 + k200 * k022 + k020 * k202 + k002 * k220 + 4 * (k011 * k211 + k101 * k121 + k110 * k112) + 2 * (k120 * k102 + k210 * k012 + k201 * k021)) * irho;
    D_Phy_DDF _b = (16 * k110 * k101 * k011 + 4 * (k101 * k101 * k020 + k011 * k011 * k200 + k110 * k110 * k002) + 2 * k200 * k020 * k002) * irho * irho;
    k222 = C222 + _a - _b;
}

/**
 * @brief step 5. store k to ddf
*/
void Collision_Cumulant::bwd_central_moment_trans(const D_Phy_DDF *k, D_Phy_DDF *ddf, const D_uvw velocity) 
{
    D_Phy_Velocity _u = velocity.x, _v = velocity.y, _w = velocity.z;

    D_Phy_DDF k1[27] = {0.};
    D_Phy_DDF k2[27] = {0.};
    D_Phy_DDF k3[27] = {0.};

    D_Phy_Velocity u2 = SQ(_u), v2 = SQ(_v), w2 = SQ(_w);

    // beta = 0, gamma = 0
    kZ00 = k000 * (1 - u2) - k100 * 2 * _u - k200;
    kN00 = (k000 * (u2 - _u) + k100 * (2 * _u - 1) + k200) * 0.5;
    kP00 = (k000 * (u2 + _u) + k100 * (2 * _u + 1) + k200) * 0.5;

    // beta = 0, gamma = 1
    kZ01 = k001 * (1 - u2) - k101 * 2 * _u - k201;
    kN01 = (k001 * (u2 - _u) + k101 * (2 * _u - 1) + k201) * 0.5;
    kP01 = (k001 * (u2 + _u) + k101 * (2 * _u + 1) + k201) * 0.5;

    // beta = 0, gamma = 2
    kZ02 = k002 * (1 - u2) - k102 * 2 * _u - k202;
    kN02 = (k002 * (u2 - _u) + k102 * (2 * _u - 1) + k202) * 0.5;
    kP02 = (k002 * (u2 + _u) + k102 * (2 * _u + 1) + k202) * 0.5;

    // beta = 1, gamma = 0
    kZ10 = k010 * (1 - u2) - k110 * 2 * _u - k210;
    kN10 = (k010 * (u2 - _u) + k110 * (2 * _u - 1) + k210) * 0.5;
    kP10 = (k010 * (u2 + _u) + k110 * (2 * _u + 1) + k210) * 0.5;

    // beta = 1, gamma = 1
    kZ11 = k011 * (1 - u2) - k111 * 2 * _u - k211;
    kN11 = (k011 * (u2 - _u) + k111 * (2 * _u - 1) + k211) * 0.5;
    kP11 = (k011 * (u2 + _u) + k111 * (2 * _u + 1) + k211) * 0.5;

    // beta = 1, gamma = 2
    kZ12 = k012 * (1 - u2) - k112 * 2 * _u - k212;
    kN12 = (k012 * (u2 - _u) + k112 * (2 * _u - 1) + k212) * 0.5;
    kP12 = (k012 * (u2 + _u) + k112 * (2 * _u + 1) + k212) * 0.5;

    // beta = 2, gamma = 0
    kZ20 = k020 * (1 - u2) - k120 * 2 * _u - k220;
    kN20 = (k020 * (u2 - _u) + k120 * (2 * _u - 1) + k220) * 0.5;
    kP20 = (k020 * (u2 + _u) + k120 * (2 * _u + 1) + k220) * 0.5;

    // beta = 2, gamma = 1
    kZ21 = k021 * (1 - u2) - k121 * 2 * _u - k221;
    kN21 = (k021 * (u2 - _u) + k121 * (2 * _u - 1) + k221) * 0.5;
    kP21 = (k021 * (u2 + _u) + k121 * (2 * _u + 1) + k221) * 0.5;

    // beta = 2, gamma = 2
    kZ22 = k022 * (1 - u2) - k122 * 2 * _u - k222;
    kN22 = (k022 * (u2 - _u) + k122 * (2 * _u - 1) + k222) * 0.5;
    kP22 = (k022 * (u2 + _u) + k122 * (2 * _u + 1) + k222) * 0.5;

    // ------- split line --------------

    // i = -1, gamma = 0
    kNZ0 = kN00 * (1 - v2) - kN10 * 2 * _v - kN20;
    kNN0 = (kN00 * (v2 - _v) + kN10 * (2 * _v - 1) + kN20) * 0.5;
    kNP0 = (kN00 * (v2 + _v) + kN10 * (2 * _v + 1) + kN20) * 0.5;

    // i = -1, gamma = 1
    kNZ1 = kN01 * (1 - v2) - kN11 * 2 * _v - kN21;
    kNN1 = (kN01 * (v2 - _v) + kN11 * (2 * _v - 1) + kN21) * 0.5;
    kNP1 = (kN01 * (v2 + _v) + kN11 * (2 * _v + 1) + kN21) * 0.5;

    // i = -1, gamma = 2
    kNZ2 = kN02 * (1 - v2) - kN12 * 2 * _v - kN22;
    kNN2 = (kN02 * (v2 - _v) + kN12 * (2 * _v - 1) + kN22) * 0.5;
    kNP2 = (kN02 * (v2 + _v) + kN12 * (2 * _v + 1) + kN22) * 0.5;

    // i = 0, gamma = 0
    kZZ0 = kZ00 * (1 - v2) - kZ10 * 2 * _v - kZ20;
    kZN0 = (kZ00 * (v2 - _v) + kZ10 * (2 * _v - 1) + kZ20) * 0.5;
    kZP0 = (kZ00 * (v2 + _v) + kZ10 * (2 * _v + 1) + kZ20) * 0.5;

    // i = 0, gamma = 1
    kZZ1 = kZ01 * (1 - v2) - kZ11 * 2 * _v - kZ21;
    kZN1 = (kZ01 * (v2 - _v) + kZ11 * (2 * _v - 1) + kZ21) * 0.5;
    kZP1 = (kZ01 * (v2 + _v) + kZ11 * (2 * _v + 1) + kZ21) * 0.5;

    // i = 0, gamma = 2
    kZZ2 = kZ02 * (1 - v2) - kZ12 * 2 * _v - kZ22;
    kZN2 = (kZ02 * (v2 - _v) + kZ12 * (2 * _v - 1) + kZ22) * 0.5;
    kZP2 = (kZ02 * (v2 + _v) + kZ12 * (2 * _v + 1) + kZ22) * 0.5;

    // i = 1, gamma = 0
    kPZ0 = kP00 * (1 - v2) - kP10 * 2 * _v - kP20;
    kPN0 = (kP00 * (v2 - _v) + kP10 * (2 * _v - 1) + kP20) * 0.5;
    kPP0 = (kP00 * (v2 + _v) + kP10 * (2 * _v + 1) + kP20) * 0.5;

    // i = 1, gamma = 1
    kPZ1 = kP01 * (1 - v2) - kP11 * 2 * _v - kP21;
    kPN1 = (kP01 * (v2 - _v) + kP11 * (2 * _v - 1) + kP21) * 0.5;
    kPP1 = (kP01 * (v2 + _v) + kP11 * (2 * _v + 1) + kP21) * 0.5;

    // i = 1, gamma = 2
    kPZ2 = kP02 * (1 - v2) - kP12 * 2 * _v - kP22;
    kPN2 = (kP02 * (v2 - _v) + kP12 * (2 * _v - 1) + kP22) * 0.5;
    kPP2 = (kP02 * (v2 + _v) + kP12 * (2 * _v + 1) + kP22) * 0.5;

    // ------- split line --------------

    // i = -1, j = -1
    kNNZ = kNN0 * (1 - w2) - kNN1 * 2 * _w - kNN2;
    kNNN = (kNN0 * (w2 - _w) + kNN1 * (2 * _w - 1) + kNN2) * 0.5;
    kNNP = (kNN0 * (w2 + _w) + kNN1 * (2 * _w + 1) + kNN2) * 0.5;

    // i = -1, j = 0
    kNZZ = kNZ0 * (1 - w2) - kNZ1 * 2 * _w - kNZ2;
    kNZN = (kNZ0 * (w2 - _w) + kNZ1 * (2 * _w - 1) + kNZ2) * 0.5;
    kNZP = (kNZ0 * (w2 + _w) + kNZ1 * (2 * _w + 1) + kNZ2) * 0.5;

    // i = -1, j = 1
    kNPZ = kNP0 * (1 - w2) - kNP1 * 2 * _w - kNP2;
    kNPN = (kNP0 * (w2 - _w) + kNP1 * (2 * _w - 1) + kNP2) * 0.5;
    kNPP = (kNP0 * (w2 + _w) + kNP1 * (2 * _w + 1) + kNP2) * 0.5;

    // i = 0, j = -1
    kZNZ = kZN0 * (1 - w2) - kZN1 * 2 * _w - kZN2;
    kZNN = (kZN0 * (w2 - _w) + kZN1 * (2 * _w - 1) + kZN2) * 0.5;
    kZNP = (kZN0 * (w2 + _w) + kZN1 * (2 * _w + 1) + kZN2) * 0.5;

    // i = 0, j = 0
    kZZZ = kZZ0 * (1 - w2) - kZZ1 * 2 * _w - kZZ2;
    kZZN = (kZZ0 * (w2 - _w) + kZZ1 * (2 * _w - 1) + kZZ2) * 0.5;
    kZZP = (kZZ0 * (w2 + _w) + kZZ1 * (2 * _w + 1) + kZZ2) * 0.5;

    // i = 0, j = 1
    kZPZ = kZP0 * (1 - w2) - kZP1 * 2 * _w - kZP2;
    kZPN = (kZP0 * (w2 - _w) + kZP1 * (2 * _w - 1) + kZP2) * 0.5;
    kZPP = (kZP0 * (w2 + _w) + kZP1 * (2 * _w + 1) + kZP2) * 0.5;

    // i = 1, j = -1
    kPNZ = kPN0 * (1 - w2) - kPN1 * 2 * _w - kPN2;
    kPNN = (kPN0 * (w2 - _w) + kPN1 * (2 * _w - 1) + kPN2) * 0.5;
    kPNP = (kPN0 * (w2 + _w) + kPN1 * (2 * _w + 1) + kPN2) * 0.5;

    // i = 1, j = 0
    kPZZ = kPZ0 * (1 - w2) - kPZ1 * 2 * _w - kPZ2;
    kPZN = (kPZ0 * (w2 - _w) + kPZ1 * (2 * _w - 1) + kPZ2) * 0.5;
    kPZP = (kPZ0 * (w2 + _w) + kPZ1 * (2 * _w + 1) + kPZ2) * 0.5;

    // i = 1, j = 1
    kPPZ = kPP0 * (1 - w2) - kPP1 * 2 * _w - kPP2;
    kPPN = (kPP0 * (w2 - _w) + kPP1 * (2 * _w - 1) + kPP2) * 0.5;
    kPPP = (kPP0 * (w2 + _w) + kPP1 * (2 * _w + 1) + kPP2) * 0.5;

    ddf[0] = kZZZ;
    ddf[1] = kPZZ;
    ddf[2] = kNZZ;
    ddf[3] = kZPZ;
    ddf[4] = kZNZ;
    ddf[5] = kZZP;
    ddf[6] = kZZN;
    ddf[7] = kPPZ;
    ddf[8] = kNNZ;
    ddf[9] = kPNZ;
    ddf[10] = kNPZ;
    ddf[11] = kPZP;
    ddf[12] = kNZN;
    ddf[13] = kPZN;
    ddf[14] = kNZP;
    ddf[15] = kZPP;
    ddf[16] = kZNN;
    ddf[17] = kZPN;
    ddf[18] = kZNP;
    ddf[19] = kPPP;
    ddf[20] = kNNN;
    ddf[21] = kPNP;
    ddf[22] = kNPN;
    ddf[23] = kNNP;
    ddf[24] = kPPN;
    ddf[25] = kNPP;
    ddf[26] = kPNN;
}