/**
 * @file Periodic.cpp
 * @brief Periodic Boundary Condition Method, copy Lattice_BC to Ghost_BC
 * @date 2023-09-22
 */
#include "BC/Periodic.hpp"

periodic_West::periodic_West()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;

    // Copy east wall to ghost bc layer in west
    for (auto i_east_bc : Lat_Manager::pointer_me->lat_bc_x.at(1))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_east_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, yint, zint)) << Morton_Assist::bit_otherlevel;
        ghost_bc_x0.insert(make_pair(ghost_code, i_east_bc.first));
    }

    for (D_uint iz = 1; iz < gNz-2; ++iz)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, 0, iz)) << Morton_Assist::bit_otherlevel; // West-South
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, 1, iz)) << Morton_Assist::bit_otherlevel; // 
        ghost_bc_x0.insert(make_pair(ghost_code, bc_code));  // copy bc_edge to West-South
        
        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, gNy-2, iz)) << Morton_Assist::bit_otherlevel; // West-North
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, gNy-3, iz)) << Morton_Assist::bit_otherlevel; // lat_bc_x x1y0
        ghost_bc_x0.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint iy = 1; iy < gNy-2; ++iy)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, iy, gNz-2)) << Morton_Assist::bit_otherlevel; // West-Top
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, iy, 1)) << Morton_Assist::bit_otherlevel; // 
        ghost_bc_x0.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, iy, 0)) << Morton_Assist::bit_otherlevel; // West-Bottom
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, iy, gNz-3)) << Morton_Assist::bit_otherlevel;
        ghost_bc_x0.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_x0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,0)        ) <<Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,gNz-3)) <<Morton_Assist::bit_otherlevel ));
    ghost_bc_x0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,gNz-2)    ) <<Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,1)    ) <<Morton_Assist::bit_otherlevel ));
    ghost_bc_x0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,gNy-2,0)    ) <<Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,1,gNz-3)    ) <<Morton_Assist::bit_otherlevel ));
    ghost_bc_x0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,gNy-2,gNz-2)) <<Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,1,1)        ) <<Morton_Assist::bit_otherlevel ));

/*
    // for (auto i_lat_bc : east_lat_bc) 
    // {
    //     if ( east_lat_bc.find( Morton_Assist::find_y1(i_lat_bc.first,0) ) == east_lat_bc.end() )
    //         ghost_bc_edge_EN.insert( Morton_Assist::find_y1( Morton_Assist::find_x1(i_lat_bc.first,0) ,0) );

    //     else if ( east_lat_bc.find( Morton_Assist::find_y0(i_lat_bc.first,0) ) == east_lat_bc.end() )
    //         ghost_bc_edge_ES.insert( Morton_Assist::find_y0( Morton_Assist::find_x1(i_lat_bc.first,0) ,0) );

    //     else if ( east_lat_bc.find( Morton_Assist::find_z1(i_lat_bc.first,0) ) == east_lat_bc.end() )
    //         ghost_bc_edge_ET.insert( Morton_Assist::find_z1( Morton_Assist::find_x1(i_lat_bc.first,0) ,0) );

    //     else if ( east_lat_bc.find( Morton_Assist::find_z0(i_lat_bc.first,0) ) == east_lat_bc.end() )
    //         ghost_bc_edge_EB.insert( Morton_Assist::find_z0( Morton_Assist::find_x1(i_lat_bc.first,0) ,0) );

    // }

    // for (auto i_edge : ghost_bc_edge_EN) {
    //     if (ghost_bc_edge_EN.find(Morton_Assist::find_z1(i_edge, 0)) == ghost_bc_edge_EN.end() )
    //         ghost_bc_corner_ETN = Morton_Assist::find_z1(i_edge, 0);

    //     else if (ghost_bc_edge_EN.find(Morton_Assist::find_z0(i_edge, 0)) == ghost_bc_edge_EN.end() )
    //         ghost_bc_corner_EBN = Morton_Assist::find_z0(i_edge, 0);
    // }

    // for (auto i_edge : ghost_bc_edge_ES) {
    //     if (ghost_bc_edge_ES.find(Morton_Assist::find_z1(i_edge, 0)) == ghost_bc_edge_ES.end() )
    //         ghost_bc_corner_ETS = Morton_Assist::find_z1(i_edge, 0);

    //     else if (ghost_bc_edge_ES.find(Morton_Assist::find_z0(i_edge, 0)) == ghost_bc_edge_ES.end())
    //         ghost_bc_corner_EBS = Morton_Assist::find_z0(i_edge, 0);
    // }
*/
}

void periodic_West::initialBCStrategy()
{
    D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
    initial_Feq(westBc_ptr);
}

void periodic_West::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_x0) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}

periodic_East::periodic_East()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;
    
    // Copy west wall to ghost bc layer in east
    for (auto i_west_bc : Lat_Manager::pointer_me->lat_bc_x.at(0))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_west_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        // std::cout << "lat_xyz " << lat_xyz.x << " , " << lat_xyz.y << " , " << lat_xyz.z << std::endl;
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, yint, zint))<< Morton_Assist::bit_otherlevel;

        // D_vec temp_xyz;
        // Morton_Assist::compute_coordinate(ghost_code, 0, temp_xyz.x, temp_xyz.y, temp_xyz.z);
        // std::cout << "temp_xyz " << temp_xyz.x << " , " << temp_xyz.y << " , " << temp_xyz.z << std::endl;
        ghost_bc_x1.insert(make_pair(ghost_code, i_west_bc.first));
    }

    for (D_uint iz = 1; iz < gNz-2; ++iz)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, 0, iz))<< Morton_Assist::bit_otherlevel; // East-South
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1, gNy-3 , iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_x1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, gNy-2, iz))<< Morton_Assist::bit_otherlevel; // East-North
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1, 1, iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_x1.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint iy = 1; iy < gNy-2; ++iy)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, iy, gNz-2))<< Morton_Assist::bit_otherlevel; // East-Top
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1, iy, 1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_x1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, iy, 0))<< Morton_Assist::bit_otherlevel; // East-Bottom
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,iy,gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_x1.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_x1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,0)        )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,gNz-3))<< Morton_Assist::bit_otherlevel));
    ghost_bc_x1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,gNz-2)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,1)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_x1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,gNy-2,0)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,gNz-3)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_x1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,gNy-2,gNz-2))<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,1)        )<< Morton_Assist::bit_otherlevel));

    for (auto i : ghost_bc_x1) {
        if (Lat_Manager::pointer_me->ghost_bc.find(i.first) == Lat_Manager::pointer_me->ghost_bc.end()) {
            D_vec xyz_temp;
            Morton_Assist::compute_coordinate(i.first, 0, xyz_temp.x, xyz_temp.y, xyz_temp.z);
            std::cout << "i.first " << i.first << " Coord: (" << xyz_temp.x << " , " << xyz_temp.y << " , " << xyz_temp.z << ")" << std::endl;
        }
        if (Lat_Manager::pointer_me->lat_f.at(0).find(i.second) == Lat_Manager::pointer_me->lat_f.at(0).end()) {
            D_vec xyz_temp;
            Morton_Assist::compute_coordinate(i.second, 0, xyz_temp.x, xyz_temp.y, xyz_temp.z);
            std::cout << "i.second " << i.second << " Coord: (" << xyz_temp.x << " , " << xyz_temp.y << " , " << xyz_temp.z << ")" << std::endl;
        }
    }
}

void periodic_East::initialBCStrategy()
{
    D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
    initial_Feq(eastBc_ptr);
}

void periodic_East::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_x1) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}

periodic_South::periodic_South()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;

    // Copy north wall to ghost bc layer in south
    for (auto i_north_bc : Lat_Manager::pointer_me->lat_bc_y.at(1))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_north_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(xint, 0, zint))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y0.insert(make_pair(ghost_code, i_north_bc.first));
    }

    for (D_uint iz = 1; iz < gNz-2; ++iz)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,iz))<< Morton_Assist::bit_otherlevel; // South-East
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y0.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,iz))<< Morton_Assist::bit_otherlevel; // South-West
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y0.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint ix = 1; ix < gNx-2; ++ix)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,0,gNz-2))<< Morton_Assist::bit_otherlevel; // South-Top
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,gNy-3,1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y0.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,0,0))<< Morton_Assist::bit_otherlevel; // South-Bottom
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,gNy-3,gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y0.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_y0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,0)        )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,gNz-3))<< Morton_Assist::bit_otherlevel));
    ghost_bc_y0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,0)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,gNz-3)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_y0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,gNz-2)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,1)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_y0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,gNz-2))<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,1)        )<< Morton_Assist::bit_otherlevel));

}

void periodic_South::initialBCStrategy()
{
    D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
    initial_Feq(southBc_ptr);
}

void periodic_South::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_y0) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}

periodic_North::periodic_North()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;

    // Copy south wall to ghost bc layer in north
    for (auto i_south_bc : Lat_Manager::pointer_me->lat_bc_y.at(0))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_south_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(xint, gNy-2, zint))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y1.insert(make_pair(ghost_code, i_south_bc.first));
    }

    for (D_uint iz = 1; iz < gNz-2; ++iz)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, gNy-2, iz))<< Morton_Assist::bit_otherlevel; // North-East
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1, iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, gNy-2, iz))<< Morton_Assist::bit_otherlevel; // North-West
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, 1, iz))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y1.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint ix = 1; ix < gNx-2; ++ix)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-2, gNz-2))<< Morton_Assist::bit_otherlevel; // North-Top
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,1, 1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-2, 0))<< Morton_Assist::bit_otherlevel; // North-Bottom
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix,1,gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_y1.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_y1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, gNy-2, 0)        )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3, 1, gNz-3))<< Morton_Assist::bit_otherlevel));
    ghost_bc_y1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0, gNy-2, gNz-2)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3 ,1,1)     )<< Morton_Assist::bit_otherlevel));
    ghost_bc_y1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, gNy-2, 0)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,gNz-3)      )<< Morton_Assist::bit_otherlevel));
    ghost_bc_y1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2, gNy-2, gNz-2))<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,1)          )<< Morton_Assist::bit_otherlevel));

}

void periodic_North::initialBCStrategy()
{
    D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
    initial_Feq(northBc_ptr);
}

void periodic_North::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_y1) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}

periodic_Top::periodic_Top()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;
    
    // Copy bottom wall to ghost bc layer in top
    for (auto i_bot_bc : Lat_Manager::pointer_me->lat_bc_z.at(0))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_bot_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint, gNz-2))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z1.insert(make_pair(ghost_code, i_bot_bc.first));
    }

    for (D_uint ix = 1; ix < gNx-2; ++ix)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-2 , gNz-2))<< Morton_Assist::bit_otherlevel; // Top-North
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, 1, 1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, 0, gNz-2))<< Morton_Assist::bit_otherlevel; // Top-South
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-3, 1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z1.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint iy = 1; iy < gNy-2; ++iy)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,iy,gNz-2))<< Morton_Assist::bit_otherlevel; // Top-West
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,iy,1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z1.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,iy,gNz-2))<< Morton_Assist::bit_otherlevel; // Top-East
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,iy,1))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z1.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_z1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,gNz-2)        )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,1))<< Morton_Assist::bit_otherlevel));
    ghost_bc_z1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,gNz-2)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,1)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_z1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,gNy-2,gNz-2)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,1,1)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_z1.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,gNy-2,gNz-2))<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,1)        )<< Morton_Assist::bit_otherlevel));


}

void periodic_Top::initialBCStrategy()
{
    D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
    initial_Feq(topBc_ptr);
}

void periodic_Top::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_z1) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}

periodic_Bot::periodic_Bot()
{
    D_uint gNx = Grid_Manager::pointer_me->nx;
    D_uint gNy = Grid_Manager::pointer_me->ny;
    D_uint gNz = Grid_Manager::pointer_me->nz;

    // Copy top wall to ghost bc layer in bottom
    for (auto i_bot_bc : Lat_Manager::pointer_me->lat_bc_z.at(1))
    {
        D_vec lat_xyz;
        Morton_Assist::compute_coordinate(i_bot_bc.first, 0, lat_xyz.x, lat_xyz.y, lat_xyz.z);
        D_vec lat_idx = lat_xyz / Lat_Manager::pointer_me->dx[0];
        D_uint xint = static_cast<D_uint>(lat_idx.x);
        D_uint yint = static_cast<D_uint>(lat_idx.y);
        D_uint zint = static_cast<D_uint>(lat_idx.z);

        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(xint, yint, 0))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z0.insert(make_pair(ghost_code, i_bot_bc.first));
    }

    for (D_uint ix = 1; ix < gNx-2; ++ix)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-2 , 0))<< Morton_Assist::bit_otherlevel; // Bot-North
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, 1, gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z0.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, 0, 0))<< Morton_Assist::bit_otherlevel; // Bot-South
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(ix, gNy-3, gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z0.insert(make_pair(ghost_code, bc_code));
    }

    for (D_uint iy = 1; iy < gNy-2; ++iy)
    {
        D_morton ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,iy,0))<< Morton_Assist::bit_otherlevel; // Bot-West
        D_morton bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,iy,gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z0.insert(make_pair(ghost_code, bc_code));

        ghost_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,iy,0))<< Morton_Assist::bit_otherlevel; // Bot-East
        bc_code = static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,iy,gNz-3))<< Morton_Assist::bit_otherlevel;
        ghost_bc_z0.insert(make_pair(ghost_code, bc_code));
    }

    ghost_bc_z0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,0,0)        )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,gNy-3,gNz-3))<< Morton_Assist::bit_otherlevel));
    ghost_bc_z0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(0,gNy-2,0)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-3,1,gNz-3)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_z0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,0,0)    )<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,gNy-3,gNz-3)    )<< Morton_Assist::bit_otherlevel));
    ghost_bc_z0.insert(make_pair(static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(gNx-2,gNy-2,0))<< Morton_Assist::bit_otherlevel, static_cast<D_morton>(Morton_Assist::pointer_me->morton_encode(1,1,gNz-3)        )<< Morton_Assist::bit_otherlevel));
}

void periodic_Bot::initialBCStrategy()
{
    D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
    initial_Feq(botBc_ptr);
}

void periodic_Bot::applyBCStrategy()
{
    for (auto i_bc : ghost_bc_z0) {
        for (int i_q = 0; i_q < C_Q; ++i_q) {
            df_ptr->f[i_q].at(i_bc.first) = df_ptr->f[i_q].at(i_bc.second);
            df_ptr->fcol[i_q].at(i_bc.first) = df_ptr->fcol[i_q].at(i_bc.second);
        }
    }
}