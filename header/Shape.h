/**
* @file
* @brief Define objective manager class.
*/

#pragma once

#include "General.h"
#include "user.h"
#include <math.h>
#include <unordered_set>

struct Solid_Node
// class Solid_Node
{
public:
	D_real x, y;
#if (C_FSI_INTERFACE == 1)
	D_real xref, yref; ///< reference coordinate at previous time step, to check if the solid point move to adjacent cells at C_max_level
	D_real xref2, yref2; ///< reference coordinate at previous time step, to check if the solid point move to adjacent cells at (C_max_level - 1)
	D_real area;     ///< area near the node, used to calculate IB force
	D_real u, v;    /// velocity of the solid point
	#endif
#if (C_DIMS==3)
	D_real z;
	// std::array<D_real,3> dir;droplet
#if (C_FSI_INTERFACE == 1)
	D_real zref, zref2, w;
#endif
#endif
// #if (C_DIMS==2)
// 	std::array<D_real, 2> dir;
// #endif
};

enum faceDir_enum
{
	non_orth = 0,
	x_plane = 1, // (+/-1, 0, 0)
	y_plane = 2, // (0, +/-1, 0)
	z_plane = 3, // (0, 0, +/-1)
};

class Box
{
public:
	Box(Solid_Node p1, Solid_Node p2)
	 : leftSouthBot({min_of_two(p1.x, p2.x), min_of_two(p1.y, p2.y), min_of_two(p1.z, p2.z)}),
	   rightNorthTop({max_of_two(p1.x, p2.x), max_of_two(p1.y, p2.y), max_of_two(p1.z, p2.z)})
	{
		extent = rightNorthTop - leftSouthBot;
		center.x = (rightNorthTop.x + leftSouthBot.x) * 0.5f;
		center.y = (rightNorthTop.y + leftSouthBot.y) * 0.5f;
		center.z = (rightNorthTop.z + leftSouthBot.z) * 0.5f;
	};

	Box(D_vec p1, D_vec p2)
	 : leftSouthBot({min_of_two(p1.x, p2.x), min_of_two(p1.y, p2.y), min_of_two(p1.z, p2.z)}),
	   rightNorthTop({max_of_two(p1.x, p2.x), max_of_two(p1.y, p2.y), max_of_two(p1.z, p2.z)})
	{
		extent = rightNorthTop - leftSouthBot;
		center.x = (rightNorthTop.x + leftSouthBot.x) * 0.5f;
		center.y = (rightNorthTop.y + leftSouthBot.y) * 0.5f;
		center.z = (rightNorthTop.z + leftSouthBot.z) * 0.5f;
	};


protected:
	D_vec leftSouthBot, rightNorthTop;
	D_vec extent;
	D_vec center;
};

class LineSegement
{
private:
	D_vec p0, p1;
	D_vec dir;
	D_real length;

public:
	LineSegement(D_vec p0, D_vec p1)
	 : dir((p1-p0).norm()), length((D_real)sqrtf64((p1.x-p0.x)*(p1.x-p0.x) + (p1.y-p0.y)*(p1.y-p0.y) + (p1.z-p0.z)*(p1.z-p0.z))) {};
};

class VertexBox : public Box
{
public:
	VertexBox() = default;

	VertexBox(D_vec p1, D_vec p2)
	 : Box(p1, p2) {};

	VertexBox(D_vec p1, D_vec p2, std::array<std::array<D_int, 3>,2> box_vertex, std::array<D_morton,2> box_vertex_code)
	 : Box(p1, p2), box_p0_x(box_vertex[0][0]), box_p0_y(box_vertex[0][1]), box_p0_z(box_vertex[0][2]), 
	   box_p1_x(box_vertex[1][0]), box_p1_y(box_vertex[1][1]), box_p1_z(box_vertex[1][2]),
	   leftSouthBot_code(box_vertex_code[0]), rightNorthTop_code(box_vertex_code[1]) {};

	D_setLat lat_in_aabb;

	D_int get_p0_x() {return box_p0_x;};
	D_int get_p0_y() {return box_p0_y;};
	D_int get_p0_z() {return box_p0_z;};

	D_int get_p1_x() {return box_p1_x;};
	D_int get_p1_y() {return box_p1_y;};
	D_int get_p1_z() {return box_p1_z;};

	D_vec get_center() {return center;};
	D_vec get_extent() {return extent;};

private:
	D_morton leftSouthBot_code, rightNorthTop_code;
	D_int box_p0_x, box_p0_y, box_p0_z, box_p1_x, box_p1_y, box_p1_z;
};

class Solid_Face
{
public:
	Solid_Node vertex1, vertex2, vertex3; ///< Three verteices of the triangle of stl.
	std::array<D_real,3> faceNorm;  ///< The norm direction of the triangle, pointing to the inside of solid.
	faceDir_enum faceDir_type;

public:
	Solid_Face triFace_offset(D_vec dir, D_real offset);
	VertexBox create_triFace_AABB();
};

/**
* @brief This enumerate will be used to classify different solid shapes.
*/
enum Shape_enum
{
	circle = 1,      // 2D
	line_fillx = 2,  // 2D and 3D
	channel = 3,     // 3D
	geofile = 4,     // read from point cloud data (2D and 3D)
	geofile_stl = 5, // read from stl file
};


/**
* @brief This structure used to store intial inforamtion of solid.
*/
struct Ini_Shape
{
	Shape_enum shape_type;
public:
	bool bool_moving = false;
	bool bool_enclosed = false;    ///< when set as true, it's a enclosed shape
	D_real x0 = 0, y0 = 0;    ///< center of each solid
#if (C_DIMS==3)
	D_real z0 = 0;
#endif
	D_uint numb_nodes = 0; ///< number of solid nodes
	D_uint numb_triFaces = 0; ///< number of solid nodes
	std::vector<D_real> length; ///< used to store character length of the solid. I.e. cycle: length[0] is the radius
};

// typedef Coord<D_real> D_vec;

/**
* @brief This class used to store all inforamtion of solid nodes.
*/
class Shape : public Ini_Shape
{
	friend class Solid_Manager;
public:
	std::vector<Solid_Node> node;
	std::vector<Solid_Face> triFace;
	// std::vector<unsigned int> triFaceIdx_of_node;  // nodes belong to which TriFace
	D_real shape_offset_x0 = 0, shape_offset_y0 = 0;
#if (C_DIMS==3)
	D_real shape_offset_z0 = 0;
#endif
	static bool judge_point_within_Shape(D_vec point, Solid_Face triFace, D_vec &intersectPoint);
	static bool judge_point_within_Shape(D_vec point, Solid_Face triFace);
	static bool intersect_line_with_triangle(D_vec p1, D_vec p2, Solid_Face triFace, Solid_Node& intersectPoint);
	static bool intersect_line_with_triangle(D_vec p1, D_vec p2, Solid_Face triFace);
	static bool intersect_line_with_triangle(D_vec startPoint, D_vec p2, Solid_Face triFace, D_real& dis_to_startPoint);
	static bool triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3]);
	static D_real two_points_length(D_vec, D_vec);
	// Shortest distance from a point to a triangle in 3D (Ericson,
	// Real-Time Collision Detection, ClosestPtPointTriangle). Used by
	// update_cutting_lattices to decide if a triangle cuts a cell — the
	// SAT-based triBoxOverlap returns false negatives for small triangles
	// (sphere.stl edge ~0.035 < dx=0.05), which starves the IB band.
	static D_real point_triangle_distance(D_vec p, Solid_Face triFace);

	void print_pointCloud();

private:
	Shape& operator=(const Ini_Shape &c1);

	void cycle(std::vector<Solid_Node> &node, D_real t);
	void line_fillx(D_real a, D_real b, std::vector<Solid_Node> &node, D_real t);
	void channel(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t);
	void geofile(D_real x0, D_real y0, D_real z0, std::vector<Solid_Node> &node, D_real t);
	void geofile_stl(D_real x0, D_real y0, D_real z0, D_real t);
	void geofile(std::vector<Solid_Node> &node, D_real t);
	// void sphere(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t);
	void sample_solidPoints_onTriFace();
};
