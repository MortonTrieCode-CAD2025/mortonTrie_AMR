/**
 * @file outflow.h
 * @brief simple exterpolation
 * @date 2023-07-24
 * @note Hard-coded D3Q19 inward-pointing direction set; gated to D3Q19 builds only.
 */
#include "BC/Outflow.hpp"
#if (C_Q == 19)
#include "Lat_Manager.hpp"
#include "Morton_Assist.h"

void outflow_West::initialBCStrategy()
{
	D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
	initial_Feq(westBc_ptr);
}

void outflow_West::applyBCStrategy()
{
	D_mapLat* westBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(0));
	for (auto i_bc_lat : *westBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_x1(loc_code, 0);
		f_ptr[1].at(loc_code) = f_ptr[1].at(nbr_code);
		f_ptr[7].at(loc_code) = f_ptr[7].at(nbr_code);
		f_ptr[9].at(loc_code) = f_ptr[9].at(nbr_code);
		f_ptr[11].at(loc_code) = f_ptr[11].at(nbr_code);
		f_ptr[13].at(loc_code) = f_ptr[13].at(nbr_code);
	}
}

void outflow_East::initialBCStrategy()
{
	D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
	initial_Feq(eastBc_ptr);
}

void outflow_East::applyBCStrategy()
{
	D_mapLat* eastBc_ptr = &(Lat_Manager::pointer_me->lat_bc_x.at(1));
	for (auto i_bc_lat : *eastBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_x0(loc_code, 0);
		f_ptr[2].at(loc_code) = f_ptr[2].at(nbr_code);
		f_ptr[8].at(loc_code) = f_ptr[8].at(nbr_code);
		f_ptr[10].at(loc_code) = f_ptr[10].at(nbr_code);
		f_ptr[12].at(loc_code) = f_ptr[12].at(nbr_code);
		f_ptr[14].at(loc_code) = f_ptr[14].at(nbr_code);
	}
}

void outflow_South::initialBCStrategy()
{
	D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
	initial_Feq(southBc_ptr);
}

void outflow_South::applyBCStrategy()
{
	D_mapLat* southBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(0));
	for (auto i_bc_lat : *southBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_y1(loc_code, 0);
		f_ptr[4].at(loc_code) = f_ptr[4].at(nbr_code);
		f_ptr[8].at(loc_code) = f_ptr[8].at(nbr_code);
		f_ptr[9].at(loc_code) = f_ptr[9].at(nbr_code);
		f_ptr[16].at(loc_code) = f_ptr[16].at(nbr_code);
		f_ptr[18].at(loc_code) = f_ptr[18].at(nbr_code);
	}
}

void outflow_North::initialBCStrategy()
{
	D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
	initial_Feq(northBc_ptr);
}

void outflow_North::applyBCStrategy()
{
	D_mapLat* northBc_ptr = &(Lat_Manager::pointer_me->lat_bc_y.at(1));
	for (auto i_bc_lat : *northBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_y0(loc_code, 0);
		f_ptr[3].at(loc_code) = f_ptr[3].at(nbr_code);
		f_ptr[7].at(loc_code) = f_ptr[7].at(nbr_code);
		f_ptr[10].at(loc_code) = f_ptr[10].at(nbr_code);
		f_ptr[15].at(loc_code) = f_ptr[15].at(nbr_code);
		f_ptr[17].at(loc_code) = f_ptr[17].at(nbr_code);
	}
}

void outflow_Bot::initialBCStrategy()
{
	D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
	initial_Feq(botBc_ptr);
}

void outflow_Bot::applyBCStrategy()
{
	D_mapLat* botBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(0));
	for (auto i_bc_lat : *botBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_z1(loc_code, 0);
		f_ptr[5].at(loc_code) = f_ptr[5].at(nbr_code);
		f_ptr[11].at(loc_code) = f_ptr[11].at(nbr_code);
		f_ptr[14].at(loc_code) = f_ptr[14].at(nbr_code);
		f_ptr[15].at(loc_code) = f_ptr[15].at(nbr_code);
		f_ptr[18].at(loc_code) = f_ptr[18].at(nbr_code);
	}
}

void outflow_Top::initialBCStrategy()
{
	D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
	initial_Feq(topBc_ptr);
}

void outflow_Top::applyBCStrategy()
{
	D_mapLat* topBc_ptr = &(Lat_Manager::pointer_me->lat_bc_z.at(1));
	for (auto i_bc_lat : *topBc_ptr)
	{
		D_morton loc_code = i_bc_lat.first;
		D_morton nbr_code = Morton_Assist::find_z0(loc_code, 0);
		f_ptr[6].at(loc_code) = f_ptr[6].at(nbr_code);
		f_ptr[12].at(loc_code) = f_ptr[12].at(nbr_code);
		f_ptr[13].at(loc_code) = f_ptr[13].at(nbr_code);
		f_ptr[16].at(loc_code) = f_ptr[16].at(nbr_code);
		f_ptr[17].at(loc_code) = f_ptr[17].at(nbr_code);
	}
}
#endif