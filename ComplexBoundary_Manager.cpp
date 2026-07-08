/**
 * @file SolidStream_Manager.h
 * @brief fluid stream aroung obstacle
 * @version 0.1
 * @date 2024-02-26
 *
 */

#include "ComplexBoundary_Manager.hpp"
#include "Lat_Manager.hpp"
#include "Morton_Assist.h"
#include "LBM_Manager.hpp"
#include <iomanip>

void UniformBoundary::treat()
{
    auto dis_ptr = &(Lat_Manager::pointer_me->dis);
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(C_max_level));

    for (auto dis_ib : *dis_ptr)
    {
        D_morton lat_code = dis_ib.first;

        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_real distance = dis_ib.second[e_inv[i_q]-1];
            if (distance != -1) {
                D_morton inv_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[e_inv[i_q]], ey[e_inv[i_q]], ez[e_inv[i_q]]);
                D_Phy_DDF fcol = user_[C_max_level].df.fcol[i_q].at(lat_code);
                D_Phy_DDF fcol_inv = user_[C_max_level].df.fcol[e_inv[i_q]].at(lat_code);
                D_Phy_DDF fcol_of_invId = user_[C_max_level].df.fcol[i_q].at(inv_code);
                user_[C_max_level].df.f[e_inv[i_q]].at(lat_code) = (distance * (fcol_inv + fcol) + (1. - distance) * fcol_of_invId) / (1. + distance);
            }
            else {
                D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[i_q], ey[i_q], ez[i_q]);
                user_[C_max_level].df.f[i_q].at(lat_code) = user_[C_max_level].df.fcol[i_q].at(neighbor_code);
            }
        }
    }

}


void SinglePointInterpolation::treat()
{
    auto dis_ptr = &(Lat_Manager::pointer_me->dis);
    D_mapLat* lat_ptr = &(Lat_Manager::pointer_me->lat_f.at(C_max_level));
    D_real dx_inner = Lat_Manager::pointer_me->dx.at(C_max_level);

    // IB hydrodynamic force (momentum-exchange)
    D_Phy_real fx_step = 0., fy_step = 0., fz_step = 0.;

    // Front/back hemisphere diagnostic
    const D_real sphere_cx = 4.0;  // sphere center x
    D_Phy_real front_fcol_sum = 0., back_fcol_sum = 0.;
    D_Phy_real front_rho_sum  = 0., back_rho_sum  = 0.;
    D_Phy_real front_u_sum    = 0., back_u_sum    = 0.;
    D_int      front_links = 0,    back_links = 0;
    D_int      front_cells = 0,    back_cells = 0;

    D_Phy_real postib_rhosum = 0., postib_usum = 0., postib_umax = 0.;
    D_int      postib_count  = 0;

    for (auto dis_ib : *dis_ptr)
    {
        D_morton lat_code = dis_ib.first;
        if (lat_ptr->find(lat_code) == lat_ptr->end()) continue;
        D_vec lat_center;
        Morton_Assist::compute_coordinate(lat_code, C_max_level, lat_center.x, lat_center.y, lat_center.z);
        lat_center += (D_vec){dx_inner/2., dx_inner/2., dx_inner/2.};
        D_real CF_U_ = LBM_Manager::pointer_me->CF_U;
        D_uvw u_b = user_[C_max_level].velocity.at(lat_code);
        D_Phy_Rho rho_b = user_[C_max_level].density.at(lat_code);

        user_[C_max_level].df.f[0].at(lat_code) = user_[C_max_level].df.fcol[0].at(lat_code);

        bool is_front = (lat_center.x < sphere_cx);
        if (is_front) { ++front_cells; front_rho_sum += rho_b; }
        else         { ++back_cells;  back_rho_sum  += rho_b; }

        for (D_int i_q = 1; i_q < C_Q; ++i_q) {
            D_real distance = dis_ib.second[e_inv[i_q]-1];

            if (distance != -1.) {

                D_vec vecR = (D_vec){0.,0.,0} - lat_center;
                D_vec angVel(0.,0.,0);
                D_vec linVel = cross_product(angVel, vecR);
                D_vec cartVel(0.,0.,0.);
                D_uvw uw((linVel.x + cartVel.x)/CF_U_, (linVel.y + cartVel.y)/CF_U_, (linVel.z + cartVel.z)/CF_U_);

                // IB wall-link interpolation: FH (Filippova-Hänel, JCP 147, 1998) or distance-weighted
                D_Phy_DDF fcol_inc = user_[C_max_level].df.fcol[e_inv[i_q]].at(lat_code);
                D_Phy_DDF f_new;
                {
                    D_morton inv_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[e_inv[i_q]], ey[e_inv[i_q]], ez[e_inv[i_q]]);
                    auto fcol_it = user_[C_max_level].df.fcol[i_q].find(inv_code);
                    D_Phy_DDF f_stream = (fcol_it != user_[C_max_level].df.fcol[i_q].end())
                                         ? fcol_it->second : fcol_inc;  // fallback to bounce-back

                    if (ib_method_ == IB_Method::FH) {
                        // FH: piecewise interpolation, continuous at d=0.5
                        constexpr D_Phy_real w[C_Q] = {
                            1./3.,
                            1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
                            1./36., 1./36., 1./36., 1./36., 1./36., 1./36.,
                            1./36., 1./36., 1./36., 1./36., 1./36., 1./36.
                        };
                        if (distance < 0.5) {
                            f_new = (1.0 - 2.0 * distance) * fcol_inc
                                  + 2.0 * distance * f_stream;
                        } else {
                            D_Phy_DDF f_eq = w[e_inv[i_q]] * rho_b;
                            D_Phy_real inv_2d = 1.0 / (2.0 * distance);
                            f_new = inv_2d * f_stream + (1.0 - inv_2d) * f_eq;
                        }
                    } else {
                        f_new = (1.0 - distance) * fcol_inc + distance * f_stream;
                    }
                }
                user_[C_max_level].df.f[i_q].at(lat_code) = f_new;
                D_Phy_DDF mterm = fcol_inc + f_new;
                fx_step += -mterm * ex[i_q];
                fy_step += -mterm * ey[i_q];
                fz_step += -mterm * ez[i_q];

                D_Phy_real umag = std::sqrt(u_b.x*u_b.x + u_b.y*u_b.y + u_b.z*u_b.z);
                if (is_front) {
                    front_fcol_sum += fcol_inc;
                    front_u_sum    += umag;
                    ++front_links;
                } else {
                    back_fcol_sum += fcol_inc;
                    back_u_sum    += umag;
                    ++back_links;
                }
            }
            else {
                // normal stream in pull mode
                D_morton neighbor_code = Morton_Assist::find_neighbor(lat_code, C_max_level, ex[i_q], ey[i_q], ez[i_q]);
                auto fcol_it = user_[C_max_level].df.fcol[i_q].find(neighbor_code);
                if (fcol_it == user_[C_max_level].df.fcol[i_q].end()) {
                    user_[C_max_level].df.f[i_q].at(lat_code) = user_[C_max_level].df.fcol[e_inv[i_q]].at(lat_code);
                } else {
                    user_[C_max_level].df.f[i_q].at(lat_code) = fcol_it->second;
                }
            }

        } // for i_q

        // Mass-conservation safety: clamp extreme ρ to [0.3, 3.0]
        D_Phy_Rho rho_after = 0.;
        for (D_int q = 0; q < C_Q; ++q)
            rho_after += user_[C_max_level].df.f[q].at(lat_code);

        if      (rho_after < 0.3)  user_[C_max_level].df.f[0].at(lat_code) += 1.0 * (0.3 - rho_after);
        else if (rho_after < 0.5)  user_[C_max_level].df.f[0].at(lat_code) += 0.5 * (0.5 - rho_after);
        else if (rho_after > 3.0)  user_[C_max_level].df.f[0].at(lat_code) += 1.0 * (3.0 - rho_after);
        else if (rho_after > 2.0)  user_[C_max_level].df.f[0].at(lat_code) += 0.5 * (2.0 - rho_after);

        // Post-IB |u|/ρ collection
        {
            D_Phy_Rho rho_pib = 0.; D_uvw u_pib = {0.,0.,0.};
            for (D_int q = 0; q < C_Q; ++q) {
                D_Phy_DDF fq = user_[C_max_level].df.f[q].at(lat_code);
                rho_pib += fq; u_pib.x += fq*ex[q]; u_pib.y += fq*ey[q]; u_pib.z += fq*ez[q];
            }
            if (rho_pib > 0.) {
                u_pib = u_pib / rho_pib;
                D_Phy_real um = std::sqrt(u_pib.x*u_pib.x + u_pib.y*u_pib.y + u_pib.z*u_pib.z);
                postib_rhosum += rho_pib; postib_usum += um;
                if (um > postib_umax) postib_umax = um;
                ++postib_count;
            }
        }

    } // for distance

    LBM_Manager::pointer_me->accumIBForce(fx_step, fy_step, fz_step);

    // MDF: Multi-Direct Forcing on wall-intersecting directions
    {
        constexpr int    MDF_ITERS = 8;
        constexpr D_Phy_real MDF_RELAX = 0.5;
        constexpr D_Phy_real w[C_Q] = {
            1./3.,
            1./18., 1./18., 1./18., 1./18., 1./18., 1./18.,
            1./36., 1./36., 1./36., 1./36., 1./36., 1./36.,
            1./36., 1./36., 1./36., 1./36., 1./36., 1./36.
        };
        for (int mdf_iter = 0; mdf_iter < MDF_ITERS; ++mdf_iter) {
            for (auto dis_ib : *dis_ptr) {
                D_morton lat_code = dis_ib.first;
                if (lat_ptr->find(lat_code) == lat_ptr->end()) continue;
                auto f0_it = user_[C_max_level].df.f[0].find(lat_code);
                if (f0_it == user_[C_max_level].df.f[0].end()) continue;

                // Compute current u from post-treatment f
                D_Phy_Rho rho = 0.; D_uvw u = {0.,0.,0.};
                for (D_int q = 0; q < C_Q; ++q) {
                    D_Phy_DDF fq = user_[C_max_level].df.f[q].at(lat_code);
                    rho += fq; u.x += fq*ex[q]; u.y += fq*ey[q]; u.z += fq*ez[q];
                }
                if (rho <= 0.) continue;
                u = u / rho;
                D_Phy_real um = std::sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
                if (um < 1e-10) continue;

                // Momentum deficit
                D_uvw dM = {-MDF_RELAX * rho * u.x,
                            -MDF_RELAX * rho * u.y,
                            -MDF_RELAX * rho * u.z};

                // Spread to wall-intersecting directions
                D_Phy_DDF df0_adj = 0.;
                for (D_int q = 1; q < C_Q; ++q) {
                    D_real dist = dis_ib.second[e_inv[q] - 1];
                    if (dist == -1.) continue;  // skip non-wall
                    D_Phy_DDF df = w[q] * 3.0 * (ex[q]*dM.x + ey[q]*dM.y + ez[q]*dM.z);
                    user_[C_max_level].df.f[q].at(lat_code) += df;
                    df0_adj -= df;
                }
                user_[C_max_level].df.f[0].at(lat_code) += df0_adj;  // ρ conservation
            }
        }
    }

    // Print front/back hemisphere diagnostic every 10 treat() calls
    {
        static D_int call_counter = 0;
        ++call_counter;
        if (call_counter % 10 == 0) {
            D_Phy_real favg_f = front_links > 0 ? front_fcol_sum / front_links : 0.;
            D_Phy_real favg_b = back_links  > 0 ? back_fcol_sum  / back_links  : 0.;
            D_Phy_real ravg_f = front_cells > 0 ? front_rho_sum  / front_cells  : 0.;
            D_Phy_real ravg_b = back_cells  > 0 ? back_rho_sum   / back_cells   : 0.;
            D_Phy_real uavg_f = front_links > 0 ? front_u_sum    / front_links  : 0.;
            D_Phy_real uavg_b = back_links  > 0 ? back_u_sum     / back_links   : 0.;
            std::cout << "\n[IB front/back] call#" << call_counter
                      << " FRONT: cells=" << front_cells << " links=" << front_links
                      << " fcol=" << favg_f << " rho=" << std::setprecision(8) << ravg_f << std::setprecision(6)
                      << " |u|=" << uavg_f
                      << " | BACK: cells=" << back_cells << " links=" << back_links
                      << " fcol=" << favg_b << " rho=" << std::setprecision(8) << ravg_b << std::setprecision(6)
                      << " |u|=" << uavg_b
                      << " | Δfcol=" << (favg_f - favg_b)
                      << " (rel=" << (favg_f > 0 ? 100.*(favg_f-favg_b)/favg_f : 0.) << "%)"
                      << " Δρ=" << std::setprecision(8) << (ravg_f - ravg_b) << std::setprecision(6) << std::endl;
            std::cout << "[IB post-treat] <|u|>=" << std::scientific << std::setprecision(3)
                      << (postib_count > 0 ? postib_usum / postib_count : 0.)
                      << " max|u|=" << postib_umax
                      << " <ρ>=" << std::fixed << std::setprecision(6)
                      << (postib_count > 0 ? postib_rhosum / postib_count : 0.)
                      << " n=" << postib_count << std::endl;
        }
    }
}
