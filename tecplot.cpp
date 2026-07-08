/**
 * @file tecplot.cpp
 * @date 2023-09-21
 */

#include "tecplot.h"
#include "Lat_Manager.hpp"
#include "LBM_Manager.hpp"
#include "Morton_Assist.h"
#include "Grid_Manager.h"
#include "Solid_Manager.h"
#include <iostream>
#include <sstream>
#include <map>

void IO_TECPLOTH::write(const D_int iter)
{
    char fileName[1024];
    sprintf(fileName, "%s/sphereFlow_%f.dat", C_OUTPUT_DIR, iter * LBM_Manager::pointer_me->CF_T[0]);
    
    D_Phy_real Convert_U = LBM_Manager::pointer_me->CF_U;

    // calculate the node numbers
    int elemNums = 0;
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        elemNums += Lat_Manager::pointer_me->lat_f.at(i_level).size() + Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level).size();
    }
    
    D_map_define<std::array<int,2>> node_idx;
    std::vector<D_morton> mortonCode_list;
    int node_idx_count = 0;
    for (auto i_level = 0; i_level <= C_max_level; ++i_level)
    {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
        {
            std::array<D_morton, 8> vertex ={
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),
                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level),i_level)
            };
            for (auto i_vet : vertex) {
                if (node_idx.find(i_vet) == node_idx.end()) {
                    ++node_idx_count;
                    node_idx.insert(make_pair(i_vet, std::array<int,2>{node_idx_count, i_level}));
                    mortonCode_list.push_back(i_vet);
                }
            }
        }
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
        {
            std::array<D_morton, 8> vertex = {
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),
                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level)
            };
            for (auto i_vet : vertex) {
                if (node_idx.find(i_vet) == node_idx.end()) {
                    ++node_idx_count;
                    node_idx.insert(make_pair(i_vet, std::array<int,2>{node_idx_count, i_level}));
                    mortonCode_list.push_back(i_vet);
                }
            }
        }
    }
    int nodeNums = node_idx.size();
    std::ofstream plt_file(fileName);
    if (!plt_file.is_open()) {
        std::cerr << "Failed to open .plt file for writing." << std::endl;
        exit(0);
    }

    plt_file << "TITLE = \"Flow Field Example\"" << std::endl;
    plt_file << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\"" << std::endl;
    plt_file << "ZONE, NODES=" << nodeNums << ", ELEMENTS=" << elemNums << ", DATAPACKING=BLOCK, VARLOCATION=([4,5,6]=CELLCENTERED), ZONETYPE=FEBRICK" << std::endl;

    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << x << std::endl;
    }
    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << y << std::endl;
    }
    for (auto i_code : mortonCode_list)
    {
        double x,y,z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        plt_file << z << std::endl;
    }

    // u
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).x * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).x * Convert_U << std::endl;
    }

    // v
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).y * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).y * Convert_U << std::endl;
    }

    // w
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).z * Convert_U << std::endl;
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level))
            plt_file << LBM_Manager::pointer_me->user[i_level].velocity.at(i_lat.first).z * Convert_U << std::endl;
    }

    // cell information
    for (int i_level = 0; i_level <= C_max_level; ++i_level) {
        for (auto i_lat : Lat_Manager::pointer_me->lat_f.at(i_level)) 
        {
            std::array<D_morton, 8> vertex = {
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),

                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level)
                
            };
            for (auto i_vex : vertex) {
                plt_file << node_idx.at(i_vex)[0] << " ";
            }
            plt_file << std::endl;
        }
        for (auto i_lat : Lat_Manager::pointer_me->lat_overlap_F2C.at(i_level)) 
        {
            std::array<D_morton, 8> vertex ={
                i_lat.first,
                Morton_Assist::find_x1(i_lat.first, i_level),
                Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_y1(i_lat.first, i_level),

                Morton_Assist::find_z1(i_lat.first, i_level),
                Morton_Assist::find_z1(Morton_Assist::find_x1(i_lat.first,i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(Morton_Assist::find_x1(i_lat.first,i_level),i_level), i_level),
                Morton_Assist::find_z1(Morton_Assist::find_y1(i_lat.first,i_level), i_level)
                
            };
            for (auto i_vex : vertex) {
                plt_file << node_idx.at(i_vex)[0] << " ";
            }
            plt_file << std::endl;
        }
    }

    plt_file.close();
}

/**
 * @brief Main function to manage mesh output in Tecplot format
 * @param[in] outfile name for output file
 */
int IO_TECPLOTH::write_mesh_manager(std::string outfile)
{
    write_mesh(outfile);
    return 0;
}

/**
 * @brief Function to write mesh data in Tecplot format
 * @param[in] outfile name for output file
 */
void IO_TECPLOTH::write_mesh(std::string outfile)
{
    char fileName[1024];
    sprintf(fileName, "%s/%s_%s.dat", C_OUTPUT_DIR, meshFileName_prefix, outfile.c_str());

    std::ofstream plt_file(fileName);
    if (!plt_file.is_open()) {
        std::cerr << "Failed to open mesh .plt file for writing: " << fileName << std::endl;
        return;
    }

    // Clear previous data
    grid_npoints.clear();
    grid_nelements.clear();

    plt_file << "TITLE = \"AMR Mesh Structure\"" << std::endl;
    plt_file << "VARIABLES = \"X\", \"Y\", \"Z\", \"Level\", \"Flag\"" << std::endl;

    int zone_count = 0;

    // Write mesh grids for each level
    for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
    {
        if (grid_vlevel.at(level_index) == C_max_level)
        {
            write_mesh_grid(Grid_Manager::pointer_me->gr_inner, level_index, plt_file, zone_count);
        }
        else
        {
            write_mesh_grid(Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)), level_index, plt_file, zone_count);
        }
    }

    // Write boundary information
    write_mesh_boundaries(plt_file, zone_count);

    // Write solid boundaries if they exist
    if (Solid_Manager::pointer_me->numb_solids > 0)
    {
        write_mesh_solid_boundaries(plt_file, zone_count);
    }

    plt_file.close();
    std::cout << "Mesh output written to: " << fileName << std::endl;
}

/**
 * @brief Function to write grid information for a specific refinement level
 * @param[in] grid_ptr pointer to the grid data structure
 * @param[in] level_index index of the refinement level
 * @param[in] plt_file output file stream
 * @param[in,out] zone_count zone counter for Tecplot zones
 */
template <class T_grptr>
void IO_TECPLOTH::write_mesh_grid(T_grptr &grid_ptr, int level_index, std::ofstream &plt_file, int &zone_count)
{
    unsigned int ilevel = grid_vlevel.at(level_index);

    // Calculate node and element counts
    unsigned int npoints = static_cast<unsigned int>(grid_ptr.grid.size());
    unsigned int nelements = 0;

    // Build node index mapping and count valid elements
    D_map_define<std::array<int,2>> node_idx;
    std::vector<D_morton> mortonCode_list;
    int node_idx_count = 0;

    // First pass: build node index and count elements
    for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
    {
        D_morton key_current = iter->first;

        // Add node to index if not already present
        if (node_idx.find(key_current) == node_idx.end()) {
            ++node_idx_count;
            node_idx.insert(std::make_pair(key_current, std::array<int,2>{node_idx_count, static_cast<int>(ilevel)}));
            mortonCode_list.push_back(key_current);
        }

        // Check if this node can form a complete element
#if (C_DIMS==2)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);

        if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
        {
            nelements++;
        }
#endif
#if (C_DIMS==3)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

        bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

        if (bool_xyz)
        {
            nelements++;
        }
#endif
    }

    grid_npoints.push_back(npoints);
    grid_nelements.push_back(nelements);

    // Write zone header
    zone_count++;
#if (C_DIMS==2)
    plt_file << "ZONE T=\"Level_" << ilevel << "\", NODES=" << npoints << ", ELEMENTS=" << nelements
             << ", DATAPACKING=BLOCK, VARLOCATION=([4,5]=CELLCENTERED), ZONETYPE=FEQUADRILATERAL" << std::endl;
#endif
#if (C_DIMS==3)
    plt_file << "ZONE T=\"Level_" << ilevel << "\", NODES=" << npoints << ", ELEMENTS=" << nelements
             << ", DATAPACKING=BLOCK, VARLOCATION=([4,5]=CELLCENTERED), ZONETYPE=FEBRICK" << std::endl;
#endif

    // Write coordinates (X, Y, Z)
    for (auto i_code : mortonCode_list)
    {
        double x, y, z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        // Apply unified domain offset for arbitrary solid positioning
        x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
        y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
        z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
        plt_file << x << std::endl;
    }

    for (auto i_code : mortonCode_list)
    {
        double x, y, z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        // Apply unified domain offset for arbitrary solid positioning
        x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
        y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
        z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
        plt_file << y << std::endl;
    }

    for (auto i_code : mortonCode_list)
    {
        double x, y, z;
        Morton_Assist::compute_coordinate(i_code, node_idx.at(i_code)[1], x, y, z);
        // Apply unified domain offset for arbitrary solid positioning
        x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
        y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
        z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
        plt_file << z << std::endl;
    }

    // Write level information (cell-centered)
    for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
    {
        D_morton key_current = iter->first;

        // Check if this node forms a valid element
#if (C_DIMS==2)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);

        if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
        {
            plt_file << ilevel << std::endl;
        }
#endif
#if (C_DIMS==3)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

        bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

        if (bool_xyz)
        {
            plt_file << ilevel << std::endl;
        }
#endif
    }

    // Write flag information (cell-centered)
    for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
    {
        D_morton key_current = iter->first;

        // Check if this node forms a valid element
#if (C_DIMS==2)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);

        if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
        {
            plt_file << iter->second.flag << std::endl;
        }
#endif
#if (C_DIMS==3)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

        bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

        if (bool_xyz)
        {
            plt_file << iter->second.flag << std::endl;
        }
#endif
    }

    // Write connectivity information
    for (auto iter = grid_ptr.grid.begin(); iter != grid_ptr.grid.end(); ++iter)
    {
        D_morton key_current = iter->first;

#if (C_DIMS==2)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);

        if ((grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
            (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()))
        {
            plt_file << node_idx.at(key_current)[0] << " "
                     << node_idx.at(key_x1)[0] << " "
                     << node_idx.at(key_xy1)[0] << " "
                     << node_idx.at(key_y1)[0] << std::endl;
        }
#endif
#if (C_DIMS==3)
        D_morton key_x1 = Morton_Assist::find_x1(key_current, ilevel);
        D_morton key_y1 = Morton_Assist::find_y1(key_current, ilevel);
        D_morton key_z1 = Morton_Assist::find_z1(key_current, ilevel);
        D_morton key_xy1 = Morton_Assist::find_y1(key_x1, ilevel);
        D_morton key_xz1 = Morton_Assist::find_z1(key_x1, ilevel);
        D_morton key_yz1 = Morton_Assist::find_z1(key_y1, ilevel);
        D_morton key_xyz1 = Morton_Assist::find_z1(key_xy1, ilevel);

        bool bool_xyz = (grid_ptr.grid.find(key_x1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_y1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_z1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xy1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_yz1) != grid_ptr.grid.end()) &&
                       (grid_ptr.grid.find(key_xyz1) != grid_ptr.grid.end());

        if (bool_xyz)
        {
            plt_file << node_idx.at(key_current)[0] << " "
                     << node_idx.at(key_x1)[0] << " "
                     << node_idx.at(key_xy1)[0] << " "
                     << node_idx.at(key_y1)[0] << " "
                     << node_idx.at(key_z1)[0] << " "
                     << node_idx.at(key_xz1)[0] << " "
                     << node_idx.at(key_xyz1)[0] << " "
                     << node_idx.at(key_yz1)[0] << std::endl;
        }
#endif
    }
}

/**
 * @brief Function to write numerical boundary information
 * @param[in] plt_file output file stream
 * @param[in,out] zone_count zone counter for Tecplot zones
 */
void IO_TECPLOTH::write_mesh_boundaries(std::ofstream &plt_file, int &zone_count)
{
    // Write numerical boundaries (interfaces between refinement levels)
    for (int level_index = 0; level_index < grid_vlevel.size(); ++level_index)
    {
        Grid * grid_ptr;
        if (grid_vlevel.at(level_index) == C_max_level)
        {
            grid_ptr = &(Grid_Manager::pointer_me->gr_inner);
        }
        else
        {
            grid_ptr = &(Grid_Manager::pointer_me->gr_NoIB.at(grid_vlevel.at(level_index)));
        }

        unsigned int ilevel = grid_vlevel.at(level_index);

        // Write fine2coarse boundaries
        for (int iboundary = 1; iboundary < 3; ++iboundary)
        {
            if (grid_ptr->fine2coarse.at(iboundary).size() > 0)
            {
                zone_count++;
                unsigned int npoints = static_cast<unsigned int>(grid_ptr->fine2coarse.at(iboundary).size());

                plt_file << "ZONE T=\"Boundary_L" << ilevel << "_F2C_" << iboundary << "\", NODES=" << npoints
                         << ", ELEMENTS=0, DATAPACKING=BLOCK, ZONETYPE=ORDERED" << std::endl;

                // Write coordinates
                for (auto iter = grid_ptr->fine2coarse.at(iboundary).begin();
                     iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
                {
                    double x, y, z;
                    Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
                    // Apply unified domain offset for arbitrary solid positioning
                    x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
                    y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
                    z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
                    plt_file << x << std::endl;
                }

                for (auto iter = grid_ptr->fine2coarse.at(iboundary).begin();
                     iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
                {
                    double x, y, z;
                    Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
                    // Apply unified domain offset for arbitrary solid positioning
                    x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
                    y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
                    z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
                    plt_file << y << std::endl;
                }

                for (auto iter = grid_ptr->fine2coarse.at(iboundary).begin();
                     iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
                {
                    double x, y, z;
                    Morton_Assist::compute_coordinate(iter->first, ilevel, x, y, z);
                    // Apply unified domain offset for arbitrary solid positioning
                    x -= Solid_Manager::pointer_me->shape_offset_x0_grid;
                    y -= Solid_Manager::pointer_me->shape_offset_y0_grid;
#if (C_DIMS==3)
                    z -= Solid_Manager::pointer_me->shape_offset_z0_grid;
#endif
                    plt_file << z << std::endl;
                }

                // Write level and flag information
                for (auto iter = grid_ptr->fine2coarse.at(iboundary).begin();
                     iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
                {
                    plt_file << ilevel << std::endl;
                }

                for (auto iter = grid_ptr->fine2coarse.at(iboundary).begin();
                     iter != grid_ptr->fine2coarse.at(iboundary).end(); ++iter)
                {
                    plt_file << iter->second << std::endl;  // boundary flag
                }
            }
        }
    }
}

/**
 * @brief Function to write solid boundary information
 * @param[in] plt_file output file stream
 * @param[in,out] zone_count zone counter for Tecplot zones
 */
void IO_TECPLOTH::write_mesh_solid_boundaries(std::ofstream &plt_file, int &zone_count)
{
    // Write solid boundary information
    for (unsigned int ishape = 0; ishape < Solid_Manager::pointer_me->numb_solids; ++ishape)
    {
        zone_count++;

        // Get solid boundary points
        const auto& solid_shape = Solid_Manager::pointer_me->shape_solids.at(ishape);

        // For now, write a simple marker zone for solid boundaries
        // This can be expanded based on the specific solid boundary data structure
        plt_file << "ZONE T=\"Solid_Boundary_" << ishape << "\", NODES=1, ELEMENTS=0, DATAPACKING=BLOCK, ZONETYPE=ORDERED" << std::endl;

        // Write a single point as a marker (this should be expanded based on actual solid data)
        plt_file << "0.0" << std::endl;  // X
        plt_file << "0.0" << std::endl;  // Y
        plt_file << "0.0" << std::endl;  // Z
        plt_file << "-1" << std::endl;   // Level (marker for solid)
        plt_file << "999" << std::endl;  // Flag (solid boundary marker)
    }
}