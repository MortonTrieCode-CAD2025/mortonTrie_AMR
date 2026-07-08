#include "LBM_Manager.hpp"
#include "lbm_kernel/Collision_Manager.hpp"
#include "lbm_kernel/Stream_Manager.hpp"
#include "ComplexBoundary_Manager.hpp"
#include "BC/D3Q19_BC_Manager.hpp"
#include "IO_Manager.h"

void LBM_Manager::fluidSimulate()
{
    initialBC();
    for (D_int i_iter = 0; i_iter < U_iters; ++i_iter) {
        LBMkernel(0);
        calculateMacros();
        infoPrint(i_iter);
        if (i_iter % U_outputIters == 0)
            IO_Manager::pointer_me->writeFlow(i_iter);
    }
}

void LBM_Manager::LBMkernel(D_int refine_level)
{
    // D_morton xxx_code("0000000000000000000000000000000000001010101001001110000110100001");

    Timer tmr;
    double t0 = tmr.elapsed();
    calculateMacros(refine_level);
    // std::cout << "[1] Before Collision A " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    collideModel->collide(refine_level);
    ++run_col_times[refine_level];
    // std::cout << "[2] After Collision A " << run_col_times[refine_level] << " time " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    double t1 = tmr.elapsed();
    printf("[level %d %d] collide[A]\t%f s.\n", refine_level, run_col_times[refine_level], t1-t0);
    for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
        D_morton lat_code = lat.first;
        double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
        double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
        if ( isnanf(ddf0_temp)
          || isnanf(ddf1_temp) ) {
            std::cout << "  - [A] Appear NaN ddf in fcol [" << refine_level << "] Code = " << lat_code  << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
        }
    }
    if (refine_level == C_max_level && run_col_times[0] <= 3) { LBM_Manager::pointer_me->diagL2MaxF("post-collideA"); LBM_Manager::pointer_me->diagL2RhoMin("post-collideA"); }

    if (refine_level < C_max_level)
        LBMkernel(refine_level+1);
    
    if (refine_level > 0) {
        // std::cout << "[3] Before AMR_transDDF_explosion A " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
        AMR_transDDF_explosion(refine_level);
        // std::cout << "[4] After AMR_transDDF_explosion A " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    }
    for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
        D_morton lat_code = lat.first;
        double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
        double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
        if ( isnanf(ddf0_temp)
          || isnanf(ddf1_temp) ) {
            std::cout << "  - [B] Appear NaN ddf in fcol [" << refine_level << "] Code = " << lat_code  << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
        }
    }

    if (refine_level == 0) {
        // NEQ extrapolation BC only needs AfterStream timing (applied to f after stream).
        // AfterCollision (applied to fcol) is intentionally omitted; add it here if ZouHe
        // or other BC strategies requiring pre-stream correction are introduced.
    }

    double t2 = tmr.elapsed();
    // std::cout << "[5] Before stream A " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    // if (refine_level == 4 && run_stm_times[refine_level]==12) {
    //     std::cout << "Wrong ddf" << std::endl;
    // }
    streamModel->stream(refine_level);
    ++run_stm_times[refine_level];
    // std::cout << "[6] After stream A " << run_stm_times[refine_level] << " time " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    double t3 = tmr.elapsed();
    printf("[level %d %d] stream[A]\t%f s.\n", refine_level, run_col_times[refine_level], t3-t2);
    for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
        D_morton lat_code = lat.first;
        double ddf0_temp = user[refine_level].df.f[0].at(lat_code);
        double ddf1_temp = user[refine_level].df.f[1].at(lat_code);
        if ( isnanf(ddf0_temp)
          || isnanf(ddf1_temp) ) {
            std::cout << "  - [C] Appear NaN ddf in f [" << refine_level << "] Code = " << lat_code  << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
        }
    }
    if (refine_level == C_max_level && run_col_times[0] <= 3) { LBM_Manager::pointer_me->diagL2MaxF("post-streamA"); LBM_Manager::pointer_me->diagL2RhoMin("post-streamA"); }

    if (refine_level == 0) {
        bc_manager_->applyBoundaryCondition("West", "AfterStream");
         bc_manager_->applyBoundaryCondition("East", "AfterStream");
          bc_manager_->applyBoundaryCondition("South", "AfterStream");
           bc_manager_->applyBoundaryCondition("North", "AfterStream");
            bc_manager_->applyBoundaryCondition("Bot", "AfterStream");
             bc_manager_->applyBoundaryCondition("Top", "AfterStream");
    }
    else if (refine_level == C_max_level) {
        solidBCModel->treat();
        // [D] NaN check on IB cells right after the IB correction (the key AMR x IB
        // coupling point). IB treat writes f[i_q] for every dis cell; if any goes
        // Coalescence NaN guard
        for (auto dis_ib : Lat_Manager::pointer_me->dis) {
            D_morton lat_code = dis_ib.first;
            double ddf0_temp = user[C_max_level].df.f[0].at(lat_code);
            double ddf1_temp = user[C_max_level].df.f[1].at(lat_code);
            if ( isnanf(ddf0_temp)
              || isnanf(ddf1_temp) ) {
                std::cout << "  - [D] Appear NaN ddf in IB f [" << C_max_level << "] Code = " << lat_code << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
                break;
            }
        }
    }
    if (refine_level == C_max_level && run_col_times[0] <= 3) { LBM_Manager::pointer_me->diagL2MaxF("post-IBtreatA"); LBM_Manager::pointer_me->diagL2RhoMin("post-IBtreatA"); }
    // for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
    //     D_morton lat_code = lat.first;
        
    //     if ( isnanf(user[refine_level].df.f->at(lat_code)) ) {
    //         std::cout << "  - [D] Appear NaN ddf in f [" << refine_level << "] Code = " << lat_code << std::endl;
    //     }
    // }

    if (refine_level < C_max_level){
        // std::cout << "[7] Before AMR_transDDF_coalescence B " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
        AMR_transDDF_coalescence(refine_level);
        // std::cout << "[8] After AMR_transDDF_coalescence B " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
        // [E] NaN check on the coarse overlap cells whose f was just rewritten by
        // /8.0 average or to bad fine-side f feeding it.
        for (auto lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(refine_level)) {
            D_morton lat_code = lat.first;
            double ddf0_temp = user[refine_level].df.f[0].at(lat_code);
            double ddf1_temp = user[refine_level].df.f[1].at(lat_code);
            if ( isnanf(ddf0_temp)
              || isnanf(ddf1_temp) ) {
                std::cout << "  - [E] Appear NaN ddf in overlap f [" << refine_level << "] Code = " << lat_code << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
                break;
            }
        }
    }
    // for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
    //     D_morton lat_code = lat.first;
    //     double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
    //     double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
    //     if ( isnanf(ddf0_temp) || ddf0_temp<0. || ddf0_temp>100. 
    //       || isnanf(ddf1_temp) || ddf1_temp<0. || ddf1_temp>100. ) {
    //         std::cout << "  - [E] Appear NaN ddf in f [" << refine_level << "] Code = " << lat_code << std::endl;
    //     }
    // }

    if (refine_level == 0)
        return;

    // std::cout << "[9] Before collision B " << " time " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    double t4 = tmr.elapsed();
    calculateMacros(refine_level);
    collideModel->collide(refine_level);
    // if (run_col_times[refine_level]==1) {
    //     std::cout << "  Access to falut" << std::endl;
    // }
    ++run_col_times[refine_level];
    double t5 = tmr.elapsed();
    std::cout << "[level " << refine_level << "] collide[B]\t" << t5-t4 << " s." << std::endl;
    if (refine_level == C_max_level && run_col_times[0] <= 3) { LBM_Manager::pointer_me->diagL2MaxF("post-collideB"); LBM_Manager::pointer_me->diagL2RhoMin("post-collideB"); }
    // std::cout << "[10] After collision B " << run_col_times[refine_level] << " time " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    // for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
    //     D_morton lat_code = lat.first;
    //     double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
    //     double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
    //     if ( isnanf(ddf0_temp) || ddf0_temp<0. || ddf0_temp>100. 
    //       || isnanf(ddf1_temp) || ddf1_temp<0. || ddf1_temp>100. ) {
    //         // [F] Appear NaN ddf in fcol [4] Code = 0000000000000000000000000000000000001010101001001111100001000010 ddf0 0.263572 ddf1 -0.0735331
    //         std::cout << "  - [F] Appear NaN ddf in fcol [" << refine_level << "] Code = " << lat_code << " ddf0 " << ddf0_temp << " ddf1 " << ddf1_temp << std::endl;
    //     }
    // }

    // After collision BC

    if (refine_level < C_max_level)
        LBMkernel(refine_level+1);

    
    if (refine_level == 0) {
        // NEQ extrapolation BC only needs AfterStream timing; AfterCollision intentionally omitted.
    }
    // std::cout << "[11] Before stream B " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    double t6 = tmr.elapsed();
    streamModel->stream(refine_level);
    ++run_stm_times[refine_level];
    // std::cout << "[12] After stream B " << run_stm_times[refine_level] << " time " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    double t7 = tmr.elapsed();
    std::cout << "[level " << refine_level << "] stream[B]\t" << t7-t6 << " s." << std::endl;
    if (refine_level == C_max_level && run_col_times[0] <= 3) { LBM_Manager::pointer_me->diagL2MaxF("post-streamB"); LBM_Manager::pointer_me->diagL2RhoMin("post-streamB"); }
    // for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
    //     D_morton lat_code = lat.first;
    //     double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
    //     double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
    //     if ( isnanf(ddf0_temp) || ddf0_temp<0. || ddf0_temp>100. 
    //       || isnanf(ddf1_temp) || ddf1_temp<0. || ddf1_temp>100. ) {
    //         std::cout << "  - [G] Appear NaN ddf in f [" << refine_level << "] Code = " << lat_code << std::endl;
    //     }
    // }
    
    if (refine_level == 0) {
        bc_manager_->applyBoundaryCondition("West", "AfterStream");
         bc_manager_->applyBoundaryCondition("East", "AfterStream");
          bc_manager_->applyBoundaryCondition("South", "AfterStream");
           bc_manager_->applyBoundaryCondition("North", "AfterStream");
            bc_manager_->applyBoundaryCondition("Bot", "AfterStream");
             bc_manager_->applyBoundaryCondition("Top", "AfterStream");
    }
    else if (refine_level == C_max_level) {
        solidBCModel->treat();
    }

    if (refine_level < C_max_level) {
        // std::cout << "[13] Before AMR_transDDF_coalescence B " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
        AMR_transDDF_coalescence(refine_level);
        // std::cout << "[14] After AMR_transDDF_coalescence B " << user[4].df.f[0].at(xxx_code) << " " << user[4].df.f[1].at(xxx_code) << std::endl;
    }
    // for (auto lat : Lat_Manager::pointer_me->lat_f.at(refine_level)) {
    //     D_morton lat_code = lat.first;
    //     double ddf0_temp = user[refine_level].df.fcol[0].at(lat_code);
    //     double ddf1_temp = user[refine_level].df.fcol[1].at(lat_code);
    //     if ( isnanf(ddf0_temp) || ddf0_temp<0. || ddf0_temp>100. 
    //       || isnanf(ddf1_temp) || ddf1_temp<0. || ddf1_temp>100. ) {
    //         std::cout << "  - [H] Appear NaN ddf in f [" << refine_level << "] Code = " << lat_code << std::endl;
    //     }
    // }
    
}

