/**
* @file
* @brief Define objective manager class.
* @note .
*/

#ifndef OBJ_MANAGER_H
#define OBJ_MANAGER_H
#include "General.h"
#include "Grid_Manager.h"
#include "Lat_Manager.hpp"
#include "IO_Manager.h"
#include "Solid_Manager.h"

class Obj_Manager
{
private:
	std::vector<unsigned int> run_order;
	
public:
	static Obj_Manager* pointer_me;         ///< pointer points to the class itself
	void initial();
	void time_marching_management();
	void output();
private:
	void time_marching(D_real sum_t, std::array<D_mapint, C_max_level + 1>  &map_add_nodes, std::array<D_mapint, C_max_level + 1>  &map_remove_nodes);
};

#endif
