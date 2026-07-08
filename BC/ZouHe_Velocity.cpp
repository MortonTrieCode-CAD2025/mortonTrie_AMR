/**
// ZouHe Velocity BC: Hecht & Harting (2010), J. Stat. Mech.

#include "BC/ZouHe_Velocity.hpp"
#if (C_Q == 19)
#include "Lat_Manager.hpp"

void nonEquilibrumBounceBack_Velocity_West::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0)); 
    for (auto i_bc_lat : *westBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        // std::cout << "bc_velocity " << bc_velocity.x << "," <<
        // bc_velocity.y << "," << bc_velocity.z << std::endl;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_West::applyBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
    for (auto i_bc_lat : *westBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        // D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        // D_Phy_DDF f2 = f_ptr[2].at(lat_code);
        // D_Phy_DDF f3 = f_ptr[3].at(lat_code);
        // D_Phy_DDF f4 = f_ptr[4].at(lat_code);
        // D_Phy_DDF f5 = f_ptr[5].at(lat_code);
        // D_Phy_DDF f6 = f_ptr[6].at(lat_code);
        // D_Phy_DDF f8 = f_ptr[8].at(lat_code);
        // D_Phy_DDF f10= f_ptr[10].at(lat_code);
        // D_Phy_DDF f12= f_ptr[12].at(lat_code);
        // D_Phy_DDF f14= f_ptr[14].at(lat_code);
        // D_Phy_DDF f15= f_ptr[15].at(lat_code);
        // D_Phy_DDF f16= f_ptr[16].at(lat_code);
        // D_Phy_DDF f17= f_ptr[17].at(lat_code);
        // D_Phy_DDF f18= f_ptr[18].at(lat_code);
        // D_Phy_Rho rho = 1. / (1. + bc_velocity.x) * (f0+f3+f4+f5+f6+f15+f17+f18+f16 + 2.*(f2+f10+f8+f14+f12));
        // D_Phy_DDF ny_x = 1./2. * (f3+f15+f17-f4-f18-f16) - 1./3. * rho * bc_velocity.y;
        // D_Phy_DDF nz_x = 1./2. * (f5+f15+f18-f6-f17-f16) - 1./3. * rho * bc_velocity.z;
        // density_ptr->at(lat_code) = rho;
        // f_ptr[1].at(lat_code) = f2 + 1./3. * rho * bc_velocity.x;
        // f_ptr[7].at(lat_code) = f8 + rho / 6. * (bc_velocity.x + bc_velocity.y) - ny_x;
        // f_ptr[9].at(lat_code) = f10 + rho / 6. * (bc_velocity.x - bc_velocity.y) + ny_x;
        // f_ptr[11].at(lat_code) = f12 + rho / 6. * (bc_velocity.x + bc_velocity.z) - nz_x;
        // f_ptr[13].at(lat_code) = f14 + rho / 6. * (bc_velocity.x - bc_velocity.z) + nz_x;

        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }

}

void nonEquilibrumBounceBack_Velocity_East::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
    for (auto i_bc_lat : *eastBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_East::applyBCStrategy()
{
    D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
    for (auto i_bc_lat : *eastBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        D_Phy_DDF f1 = f_ptr[1].at(lat_code);
        D_Phy_DDF f3 = f_ptr[3].at(lat_code);
        D_Phy_DDF f4 = f_ptr[4].at(lat_code);
        D_Phy_DDF f5 = f_ptr[5].at(lat_code);
        D_Phy_DDF f6 = f_ptr[6].at(lat_code);
        D_Phy_DDF f7 = f_ptr[7].at(lat_code);
        D_Phy_DDF f9 = f_ptr[9].at(lat_code);
        D_Phy_DDF f11= f_ptr[11].at(lat_code);
        D_Phy_DDF f13= f_ptr[13].at(lat_code);
        D_Phy_DDF f15= f_ptr[15].at(lat_code);
        D_Phy_DDF f16= f_ptr[16].at(lat_code);
        D_Phy_DDF f17= f_ptr[17].at(lat_code);
        D_Phy_DDF f18= f_ptr[18].at(lat_code);
        D_Phy_Rho rho = 1. / (1. - bc_velocity.x) * (f0+f3+f4+f5+f6+f15+f17+f18+f16 + 2.*(f1+f7+f9+f11+f13));
        D_Phy_DDF ny_x = 1./2. * (f3+f15+f17-f4-f18-f16) - 1./3. * rho * bc_velocity.y;
        D_Phy_DDF nz_x = 1./2. * (f5+f15+f18-f6-f17-f16) - 1./3. * rho * bc_velocity.z;
        density_ptr->at(lat_code) = rho;
        f_ptr[2].at(lat_code) = f1 - 1./3. * rho * bc_velocity.x;
        f_ptr[8].at(lat_code) = f7 - rho / 6. * (bc_velocity.x + bc_velocity.y) + ny_x;
        f_ptr[10].at(lat_code) = f9 + rho / 6. * (bc_velocity.y - bc_velocity.x) - ny_x;
        f_ptr[12].at(lat_code) = f11 - rho / 6. * (bc_velocity.x + bc_velocity.z) + nz_x;
        f_ptr[14].at(lat_code) = f13 + rho / 6. * (bc_velocity.z - bc_velocity.x) - nz_x;
    }
}

void nonEquilibrumBounceBack_Velocity_South::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
    for (auto i_bc_lat : *southBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_South::applyBCStrategy()
{
    D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
    for (auto i_bc_lat : *southBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        D_Phy_DDF f1 = f_ptr[1].at(lat_code);
        D_Phy_DDF f2 = f_ptr[2].at(lat_code);
        D_Phy_DDF f4 = f_ptr[4].at(lat_code);
        D_Phy_DDF f5 = f_ptr[5].at(lat_code);
        D_Phy_DDF f6 = f_ptr[6].at(lat_code);
        D_Phy_DDF f8 = f_ptr[8].at(lat_code);
        D_Phy_DDF f9 = f_ptr[9].at(lat_code);
        D_Phy_DDF f11= f_ptr[11].at(lat_code);
        D_Phy_DDF f12= f_ptr[12].at(lat_code);
        D_Phy_DDF f13= f_ptr[13].at(lat_code);
        D_Phy_DDF f14= f_ptr[14].at(lat_code);
        D_Phy_DDF f16= f_ptr[16].at(lat_code);
        D_Phy_DDF f18= f_ptr[18].at(lat_code);
        D_Phy_Rho rho = 1. / (1. - bc_velocity.y) * (f0+f1+f2+f5+f6+f11+f13+f14+f12 + 2.*(f4+f9+f8+f18+f16));
        D_Phy_DDF nx_y = 1./2. * (f1+f11+f13-f2-f14-f12) - 1./3. * rho * bc_velocity.x;
        D_Phy_DDF nz_y = 1./2. * (f5+f11+f14-f6-f13-f12) - 1./3. * rho * bc_velocity.z;
        density_ptr->at(lat_code) = rho;
        f_ptr[3].at(lat_code) = f4 + 1./3. * rho * bc_velocity.y;
        f_ptr[7].at(lat_code) = f8 + rho / 6. * (bc_velocity.x + bc_velocity.y) - nx_y;
        f_ptr[10].at(lat_code) = f9 + rho / 6. * (bc_velocity.y - bc_velocity.x) + nx_y;
        f_ptr[15].at(lat_code) = f16 + rho / 6. * (bc_velocity.y + bc_velocity.z) - nz_y;
        f_ptr[17].at(lat_code) = f18 + rho / 6. * (bc_velocity.y - bc_velocity.z) + nz_y;
    }
}

void nonEquilibrumBounceBack_Velocity_North::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
    for (auto i_bc_lat : *northBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_North::applyBCStrategy()
{
    D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
    for (auto i_bc_lat : *northBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        D_Phy_DDF f1 = f_ptr[1].at(lat_code);
        D_Phy_DDF f2 = f_ptr[2].at(lat_code);
        D_Phy_DDF f3 = f_ptr[3].at(lat_code);
        D_Phy_DDF f5 = f_ptr[5].at(lat_code);
        D_Phy_DDF f6 = f_ptr[6].at(lat_code);
        D_Phy_DDF f7 = f_ptr[7].at(lat_code);
        D_Phy_DDF f10= f_ptr[10].at(lat_code);
        D_Phy_DDF f11= f_ptr[11].at(lat_code);
        D_Phy_DDF f12= f_ptr[12].at(lat_code);
        D_Phy_DDF f13= f_ptr[13].at(lat_code);
        D_Phy_DDF f14= f_ptr[14].at(lat_code);
        D_Phy_DDF f15= f_ptr[15].at(lat_code);
        D_Phy_DDF f17= f_ptr[17].at(lat_code);
        D_Phy_Rho rho = 1. / (1. + bc_velocity.y) * (f0+f1+f2+f5+f6+f11+f13+f14+f12 + 2.*(f3+f7+f10+f15+f17));
        D_Phy_DDF nx_y = 1./2. * (f1+f11+f13-f2-f14-f12) - 1./3. * rho * bc_velocity.x;
        D_Phy_DDF nz_y = 1./2. * (f5+f11+f14-f6-f13-f12) - 1./3. * rho * bc_velocity.z;
        density_ptr->at(lat_code) = rho;
        f_ptr[4].at(lat_code) = f3 - 1./3. * rho * bc_velocity.y;
        f_ptr[8].at(lat_code) = f7 - rho / 6. * (bc_velocity.x + bc_velocity.y) + nx_y;
        f_ptr[9].at(lat_code) = f10 + rho / 6. * (bc_velocity.x - bc_velocity.y) - nx_y;
        f_ptr[16].at(lat_code) = f15 - rho / 6. * (bc_velocity.z + bc_velocity.y) + nz_y;
        f_ptr[18].at(lat_code) = f17 + rho / 6. * (bc_velocity.z - bc_velocity.y) - nz_y;
    }
}

void nonEquilibrumBounceBack_Velocity_Top::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
    for (auto i_bc_lat : *topBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_Top::applyBCStrategy()
{
    D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
    for (auto i_bc_lat : *topBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        D_Phy_DDF f1 = f_ptr[1].at(lat_code);
        D_Phy_DDF f2 = f_ptr[2].at(lat_code);
        D_Phy_DDF f3 = f_ptr[3].at(lat_code);
        D_Phy_DDF f4 = f_ptr[4].at(lat_code);
        D_Phy_DDF f5 = f_ptr[5].at(lat_code);
        D_Phy_DDF f7 = f_ptr[7].at(lat_code);
        D_Phy_DDF f8 = f_ptr[8].at(lat_code);
        D_Phy_DDF f9 = f_ptr[9].at(lat_code);
        D_Phy_DDF f10= f_ptr[10].at(lat_code);
        D_Phy_DDF f11= f_ptr[11].at(lat_code);
        D_Phy_DDF f14= f_ptr[14].at(lat_code);
        D_Phy_DDF f15= f_ptr[15].at(lat_code);
        D_Phy_DDF f18= f_ptr[18].at(lat_code);
        D_Phy_Rho rho = 1. / (1. + bc_velocity.z) * (f0+f1+f2+f3+f4+f7+f9+f10+f8 + 2.*(f5+f11+f14+f15+f18));
        D_Phy_DDF nx_z = 1./2. * (f1+f7+f9-f2-f10-f8) - 1./3. * rho * bc_velocity.x;
        D_Phy_DDF ny_z = 1./2. * (f3+f7+f10-f4-f9-f8) - 1./3. * rho * bc_velocity.y;
        density_ptr->at(lat_code) = rho;
        f_ptr[6].at(lat_code) = f5 - 1./3. * rho * bc_velocity.z;
        f_ptr[12].at(lat_code) = f11 - rho / 6. * (bc_velocity.z + bc_velocity.x) + nx_z;
        f_ptr[13].at(lat_code) = f14 + rho / 6. * (bc_velocity.x - bc_velocity.z) - nx_z;
        f_ptr[16].at(lat_code) = f15 - rho / 6. * (bc_velocity.z + bc_velocity.y) + ny_z;
        f_ptr[17].at(lat_code) = f18 + rho / 6. * (bc_velocity.y - bc_velocity.z) - ny_z;
    }
}

void nonEquilibrumBounceBack_Velocity_Bot::initialBCStrategy()
{
    D_Phy_DDF feq_temp[C_Q];
    D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
    for (auto i_bc_lat : *botBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        velocity_ptr->at(lat_code) = bc_velocity;
        calculateFeq(density_ptr->at(lat_code), velocity_ptr->at(lat_code), feq_temp);
        df_ptr->fcol[0].at(lat_code)  = feq_temp[0];
        df_ptr->fcol[1].at(lat_code)  = feq_temp[1];
        df_ptr->fcol[2].at(lat_code)  = feq_temp[2];
        df_ptr->fcol[3].at(lat_code)  = feq_temp[3];
        df_ptr->fcol[4].at(lat_code)  = feq_temp[4];
        df_ptr->fcol[5].at(lat_code)  = feq_temp[5];
        df_ptr->fcol[6].at(lat_code)  = feq_temp[6];
        df_ptr->fcol[7].at(lat_code)  = feq_temp[7];
        df_ptr->fcol[8].at(lat_code)  = feq_temp[8];
        df_ptr->fcol[9].at(lat_code)  = feq_temp[9];
        df_ptr->fcol[10].at(lat_code) = feq_temp[10];
        df_ptr->fcol[11].at(lat_code) = feq_temp[11];
        df_ptr->fcol[12].at(lat_code) = feq_temp[12];
        df_ptr->fcol[13].at(lat_code) = feq_temp[13];
        df_ptr->fcol[14].at(lat_code) = feq_temp[14];
        df_ptr->fcol[15].at(lat_code) = feq_temp[15];
        df_ptr->fcol[16].at(lat_code) = feq_temp[16];
        df_ptr->fcol[17].at(lat_code) = feq_temp[17];
        df_ptr->fcol[18].at(lat_code) = feq_temp[18];

        df_ptr->f[0].at(lat_code)  = feq_temp[0];
        df_ptr->f[1].at(lat_code)  = feq_temp[1];
        df_ptr->f[2].at(lat_code)  = feq_temp[2];
        df_ptr->f[3].at(lat_code)  = feq_temp[3];
        df_ptr->f[4].at(lat_code)  = feq_temp[4];
        df_ptr->f[5].at(lat_code)  = feq_temp[5];
        df_ptr->f[6].at(lat_code)  = feq_temp[6];
        df_ptr->f[7].at(lat_code)  = feq_temp[7];
        df_ptr->f[8].at(lat_code)  = feq_temp[8];
        df_ptr->f[9].at(lat_code)  = feq_temp[9];
        df_ptr->f[10].at(lat_code) = feq_temp[10];
        df_ptr->f[11].at(lat_code) = feq_temp[11];
        df_ptr->f[12].at(lat_code) = feq_temp[12];
        df_ptr->f[13].at(lat_code) = feq_temp[13];
        df_ptr->f[14].at(lat_code) = feq_temp[14];
        df_ptr->f[15].at(lat_code) = feq_temp[15];
        df_ptr->f[16].at(lat_code) = feq_temp[16];
        df_ptr->f[17].at(lat_code) = feq_temp[17];
        df_ptr->f[18].at(lat_code) = feq_temp[18];
    }
}

void nonEquilibrumBounceBack_Velocity_Bot::applyBCStrategy()
{
    D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
    for (auto i_bc_lat : *botBc_ptr)
    {
        D_morton lat_code = i_bc_lat.first;
        D_Phy_DDF f0 = f_ptr[0].at(lat_code);
        D_Phy_DDF f1 = f_ptr[1].at(lat_code);
        D_Phy_DDF f2 = f_ptr[2].at(lat_code);
        D_Phy_DDF f3 = f_ptr[3].at(lat_code);
        D_Phy_DDF f4 = f_ptr[4].at(lat_code);
        D_Phy_DDF f6 = f_ptr[6].at(lat_code);
        D_Phy_DDF f7 = f_ptr[7].at(lat_code);
        D_Phy_DDF f8 = f_ptr[8].at(lat_code);
        D_Phy_DDF f9 = f_ptr[9].at(lat_code);
        D_Phy_DDF f10= f_ptr[10].at(lat_code);
        D_Phy_DDF f12= f_ptr[12].at(lat_code);
        D_Phy_DDF f13= f_ptr[13].at(lat_code);
        D_Phy_DDF f16= f_ptr[16].at(lat_code);
        D_Phy_DDF f17= f_ptr[17].at(lat_code);
        D_Phy_Rho rho = 1. / (1. - bc_velocity.z) * (f0+f1+f2+f3+f4+f7+f9+f10+f8 + 2.*(f6+f13+f12+f17+f16));
        D_Phy_DDF nx_z = 1./2. * (f1+f7+f9-f2-f10-f8) - 1./3. * rho * bc_velocity.x;
        D_Phy_DDF ny_z = 1./2. * (f3+f7+f10-f4-f9-f8) - 1./3. * rho * bc_velocity.y;
        density_ptr->at(lat_code) = rho;
        f_ptr[5].at(lat_code) = f6 + 1./3. * rho * bc_velocity.z;
        f_ptr[11].at(lat_code) = f12 + rho / 6. * (bc_velocity.z + bc_velocity.x) - nx_z;
        f_ptr[14].at(lat_code) = f13 + rho / 6. * (bc_velocity.z - bc_velocity.x) + nx_z;
        f_ptr[15].at(lat_code) = f16 + rho / 6. * (bc_velocity.z + bc_velocity.y) - ny_z;
        f_ptr[18].at(lat_code) = f17 + rho / 6. * (bc_velocity.z - bc_velocity.y) + ny_z;
    }
}
#endif