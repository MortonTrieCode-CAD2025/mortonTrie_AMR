/**
 * @file Lat_Manager.h
 * @date 2023-06-11
 */
#ifndef LAT_MANAGER_H
#define LAT_MANAGER_H

#include "Constants.h"
#include "Grid_Class.h"
#include "Grid_Manager.h"
class Solid_Face;
// #include <set>
#include "user.h"
// #define SPHERE_TEST
// #define OLD_VERSION_LAT
#define BC_NO_GHOST

// enum class Cell_Flag
enum Cell_Flag : int
{

/*
    FLUID (Pure fluid)

        7________8
       /:       /|
      5_:______6 |
      | :  *   | |
      | 3------|-4
      | /      | /
      1/_______2/   

    IB (stl Triangle face cutting a corner of lattice)
           ______________
        7__\            /
       /:   \          /
      5_:____\        /
      | :  *  \      /
      | 3------\    /
      | /      |\  /
      1/_______2/\/  

    SURFACE (stl Triangle face immersing over the center of lattice)
      ______________
      \            /
      /\          /
      | \        /
      |  \      /
      | / \    /
      |/___\  /
            \/

    SOLID (inside the shape body)    
      ______________  
      \            /    7________8
       \          /    /:       /|
        \        /    5_:______6 |
         \      /     | :  *   | |
          \    /      | 3------|-4
           \  /       | /      | /
            \/        1/_______2/   

*/


    FLUID = 2, // flag_refine
    IB = 16,   
    SOLID = 4, // flag_ghost (inside shape)
    SURFACE = 8, 
    OVERLAP_C2F = -1, 
    OVERLAP_F2C = -3,

// Doundary set
#if C_DIMS == 3
    BDRY_FACE_WEST  = 10,
    BDRY_FACE_EAST  = 11,
    BDRY_FACE_SOUTH = 12,
    BDRY_FACE_NORTH = 13,
    BDRY_FACE_BOT   = 14,
    BDRY_FACE_TOP   = 15,

    BDRY_GHOST = -4,
    OVERLAP_GHOST = -5,

    // BDRY_GHOST_WEST  = -4,
    // BDRY_GHOST_EAST  = -5,
    // BDRY_GHOST_SOUTH = -6,
    // BDRY_GHOST_NORTH = -7,
    // BDRY_GHOST_BOT   = -8,
    // BDRY_GHOST_TOP   = -9,

    // BDRY_GHOST_NE    = -10,
    // BDRY_GHOST_NW    = -11,
    // BDRY_GHOST_SW    = -12,
    // BDRY_GHOST_SE    = -13,

    // BDRY_GHOST_TN    = -14,
    // BDRY_GHOST_TS    = -15,
    // BDRY_GHOST_BN    = -16,
    // BDRY_GHOST_BS    = -17,

    // BDRY_GHOST_TE    = -18,
    // BDRY_GHOST_TW    = -19,
    // BDRY_GHOST_BE    = -20,
    // BDRY_GHOST_BW    = -21,

    // BDRY_GHOST_NET   = -22,
    // BDRY_GHOST_NWT   = -23,
    // BDRY_GHOST_SWT   = -24,
    // BDRY_GHOST_SET   = -25,
    // BDRY_GHOST_NEB   = -26,
    // BDRY_GHOST_NWB   = -27,
    // BDRY_GHOST_SWB   = -28,
    // BDRY_GHOST_SEB   = -29,

        // Considering Edge & Corner
    #ifdef BC_EDGE
        BDRY_EDGE_NE = 16, // North-East   ++0
        BDRY_EDGE_NW = 17, // North-West   +-0
        BDRY_EDGE_SW = 18, // South-West   --0
        BDRY_EDGE_SE = 19, // South-East   -+0
        BDRY_EDGE_TN = 20, // Top-North    0++
        BDRY_EDGE_TS = 21, // Top-South    0-+
        BDRY_EDGE_BN = 22, // Bottom-North 0+-
        BDRY_EDGE_BS = 23, // Bottom-South 0--
        BDRY_EDGE_TE = 24, // Top-East     +0+
        BDRY_EDGE_TW = 25, // Top-West     -0+
        BDRY_EDGE_BE = 26, // Bottom-East  +0-
        BDRY_EDGE_BW = 27, // Bottom-West  -0-

        // Convex corner
        BDRY_CORNER_NET = 28, // North-East-Top    +++
        BDRY_CORNER_NWT = 29, // North-West-Top    -++
        BDRY_CORNER_SWT = 30, // South-West-Top    --+
        BDRY_CORNER_SET = 31, // South-East-Top    +-+
        BDRY_CORNER_NEB = 32, // North-East-Bottom ++-
        BDRY_CORNER_NWB = 33, // North-West-Bottom -+-
        BDRY_CORNER_SWB = 34, // South-West-Bottom ---
        BDRY_CORNER_SEB = 35, // South-East-Bottom +--

        // Concave corner
        /**  @todo 
        */
    #endif

#endif

#if C_DIMS == 2
    BDRY_EDGE_LEFT  = 10, // Left
    BDRY_EDGE_RIGHT = 11, // Right
    BDRY_EDGE_TOP   = 12, // Top
    BDRY_EDGE_BOT   = 13, // Bottom

    BDRY_CORNER_RT  = 14, // Right-Top    ++
    BDRY_CORNER_RB  = 15, // Right-Bottom +-
    BDRY_CORNER_LT  = 16, // Left-Top     -+
    BDRY_CORNER_LB  = 17, // Left-Bottom  --
#endif

    TBD = 9999, // to be determined
};

typedef struct Cell
{
public:
    int flag = 0;
    // D_uint index;
    /**
               7________8
              /:       /|
             5_:______6 |
             | :      | |
             | 3------|-4
             | /      | /
base_point-> 1/_______2/   
     */

    Cell(const int& flag = 0/*, const D_uint& index*/)
        : flag(flag)/*, index(index) */
    {}
} Cell;

using D_mapLat = D_map_define<Cell>;

class Lat_Manager
{
private:
    D_map_define<std::vector<std::array<D_int,2>>> triFaceID_in_srf; // record the triFace id inside the cell [MortonCode][all_triFace][solid_id, triFace]
#ifdef OLD_VERSION_LAT
    void initial_IBdistance(const D_map_define<std::vector<std::set<unsigned int>>> &);
#endif
#ifdef SPHERE_TEST
    bool get_intersectPoint_withSphere(const D_vec& IBpoint, const D_vec& solidpoint, std::vector<D_vec>& intersectPoint, const D_vec& sphere_center = D_vec(4.,4.,4.), double r=0.5);
#endif
    // void initial_triface_in_Lat(D_map_define<std::vector<std::set<unsigned int>>> &);

    void update_cutting_lattices();
    void propagate_fluid_lattices();
    void update_solid_lattices();
    void update_innerLevel_lattices();
    void update_outerLevel_lattices();
    void update_boundary_lattices();
    void update_overlap_lattices();

	void compute_distance(Solid_Face triFace, D_morton local_code, D_vec lat_center);

    D_morton dequeue(D_int &front, const D_int &rear, const std::vector<D_morton> &queue);
    void enqueue(D_morton cur_code, D_int &front, D_int &rear, std::vector<D_morton> &queue, D_map_define<bool> &visited);

private:
    enum Bdry_Type : uint {
        BDRY_FACE_WEST  = 0,
        BDRY_FACE_EAST  = 1,
        BDRY_FACE_SOUTH = 2,
        BDRY_FACE_NORTH = 3,
        BDRY_FACE_BOT   = 4,
        BDRY_FACE_TOP   = 5,

        BDRY_EDGE_EN = 6,
        BDRY_EDGE_WN = 7,
        BDRY_EDGE_WS = 8,
        BDRY_EDGE_ES = 9,
        BDRY_EDGE_NT = 10,
        BDRY_EDGE_ST = 11,
        BDRY_EDGE_NB = 12,
        BDRY_EDGE_SB = 13,
        BDRY_EDGE_ET = 14,
        BDRY_EDGE_WT = 15,
        BDRY_EDGE_EB = 16,
        BDRY_EDGE_WB = 17,

        BDRY_CORNER_ENT = 18,
        BDRY_CORNER_WNT = 19,
        BDRY_CORNER_WST = 20,
        BDRY_CORNER_EST = 21,
        BDRY_CORNER_ENB = 22,
        BDRY_CORNER_WNB = 23,
        BDRY_CORNER_WSB = 24,
        BDRY_CORNER_ESB = 25 };

public:
    static Lat_Manager* pointer_me;

    // lat_f, lat_overlap & lat_sf play the role as the lattice index lists
    std::array<D_mapLat, C_max_level+1> lat_f;  ///> The Hash table storing the fluid (FLUID & IB)

    /**
     *    Real Lattice domain
     *       _________ _________ _________ ____ ____ ____ ____ __ __ __ __ 
     *      |         |         |         | F  | F  | F  | F  |__|__|__|__|
     *      |    C    |    C    |    C    |____|____|____|____|__|__|__|__|
     *      |         |         |         | F  | F  | F  | F  |__|__|__|__|
     *      |_________|_________|_________|____|____|____|____|__|__|__|__|
     *
     * 
     *    Overlap for interface interpolation from Coarse to Fine
     *     Ghost_C2F layer to be a capacity of providing the DDF to overlap_C2F during stream
     *     ghost_C2F also need explosion from a inner Real Coarse Lattice, but donot execute collision & stream
     *      ............... ____ ____ ____ ............ __ __ __ ............
     *      :         :    |ghos|C2F |C2F |    :    :..|__|_C|_C|..:..:..:..:
     *      :         :    |____|____|____|....:....:..|__|_2|_2|..:..:..:..:
     *      :         :    |tC2F|C2F |C2F |    :    :..|__|_F|_F|..:..:..:..:
     *      :.........:....|____|____|____|....:....:..|__|__|__|..:..:..:..:
     *                                                  ^
     *                                               ghost_C2F
     * 
     * 
     *    Overlap to find C2F belonging to which Real Coarse Lattice during interface interpolation
     *      overlap_F2C only collision, stream by coalescence implicitly, this makes overlap_F2C be removed from lat_f[coarse]
     *        overlap_F2C is assigned to find overlap_C2F
     *        ghost_F2C is assigned to find overlap ghost_C2F
     *      .................... _________ .............. ____ ............
     *      :         :         |         |    :    :    |F2C |..:..:..:..:
     *      :         :         |   F2C   |....:....:....|____|..:..:..:..:
     *      :         :         |         |    :    :    |F2C |..:..:..:..:
     *      :.........:.........|_________|....:....:....|____|..:..:..:..:
     * 
     *      .......... _________ ................... ____ .................
     *      :         |         |         :    :    |    |    :..:..:..:..:
     *      :         |ghost_F2C|         :....:....|____|....:..:..:..:..:
     *      :         |         |         :    :    |    |    :..:..:..:..:
     *      :.........|_________|.........:....:....|____|....:..:..:..:..:
     *                                                ^
     *                                             ghost_F2C
     * 
     */

    std::array<D_mapLat, C_max_level+1> lat_overlap_F2C; ///> The Hash table storing the overlapping layer (OVERLAP_F2C)
    std::array<D_mapLat, C_max_level+1> lat_overlap_C2F; ///> The Hash table storing the overlapping layer (OVERLAP_C2F)
    std::array<D_mapLat, C_max_level+1> ghost_overlap_C2F;
    std::array<D_mapLat, C_max_level+1> ghost_overlap_F2C;
    D_mapLat lat_sf; ///> The Hash table storing the SURFACE & SOLID, they don't work in the collision&stream, but alloced in the user[].df memory.

#ifndef BC_NO_GHOST
    std::array<D_mapLat, 2> lat_bc_x;  // lat_bc_xFace[0]:x=0(West)  lat_bc_xFace[1]:x=C_xb(East)
    std::array<D_mapLat, 2> lat_bc_y;  // lat_bc_yFace[0]:y=0(South) lat_bc_yFace[1]:y=C_yb(North)
    std::array<D_mapLat, 2> lat_bc_z;  // lat_bc_zFace[0]:z=0(Bot)   lat_bc_zFace[1]:z=C_zb(Top)
#else
    // std::tuple<std::array<D_setLat, 6>, std::array<D_setLat, 12>, std::array<D_setLat, 8>> lat_bc;
    std::array<D_setLat, 26> lat_bc;
#endif

#ifndef BC_NO_GHOST
    D_mapLat ghost_bc;
#endif

#ifdef OLD_VERSION_LAT
    DEPRECATED void initial();
#endif
    void voxelize();
    // void output_hdf5();

    D_map_define<std::array<D_real, C_Q-1>> dis;

    D_int nNodeX, nNodeY, nNodeZ;
    std::array<D_real, C_max_level+1> dx;
};

#endif