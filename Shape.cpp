/**
* @file
* @brief Main function.
* @note  If shape.enclosed is true, to use flood fill to delte nodes inside the solid boundary, 
*         the line connecting the first solid node and the center of the solid should be inside the solid boundary
*/
#include "Shape.h"
#include "user.h"
#include <algorithm>
#include "Morton_Assist.h"
#ifdef __GNUC__
#define unreachable() (__builtin_unreachable())
#endif
#ifdef __MSC_VER
#define unreachable()(__assume(false))
#endif

Shape& Shape::operator=(const Ini_Shape &c2)
{
	this->shape_type = c2.shape_type;
	this->bool_moving = c2.bool_moving;
	this->x0 = c2.x0;
	this->y0 = c2.y0;
#if (C_DIMS==3)
	this->z0 = c2.z0;
#endif
	this->numb_nodes = c2.numb_nodes;
	this->length = c2.length;
	return *this;
}

#if (C_DIMS == 2)
/**
* @brief      function to generate 2D circle.
* @param[in]  radius     radius of the circle.
* @param[out] node  points on the circle.
* @param[in]  t     time.
*/
void Shape::cycle(std::vector<Solid_Node> &node, D_real t)
{
	bool_enclosed = true;
	x0 = length.at(0);
	y0 = length.at(1);
	D_real radius = length.at(2);

	if ((2 * C_pi*radius / static_cast<D_real>(numb_nodes)) > (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))))
	{
		std::stringstream warning;
		warning << "the distance between two solid points of the cycle is greater than the grid space at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	D_real radius_final = 2*radius - radius * sin(2. * C_pi*t/C_dx/80);

	for (D_uint i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x = radius_final * cos(2. * C_pi * (static_cast<D_real>(i) / static_cast<D_real>(numb_nodes))) + x0 + shape_offset_x0;
		node.at(i).y = radius_final * sin(2. * C_pi * (static_cast<D_real>(i) / static_cast<D_real>(numb_nodes))) + y0 + shape_offset_y0;
#if (C_FSI_INTERFACE == 1)
		node.at(i).u = 0;
		node.at(i).v = 0;
		node.at(i).area = 2 * C_pi*radius_final / static_cast<D_real>(numb_nodes);
#endif
	}

}
#endif

/**
* @brief      function to generate 2D line.
* @param[in]  a     constant for y = ax + b.
* @param[in]  b     constant for y = ax + b.
* @param[out] node  points on the line.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. When in 3D, the z is set as C_zb / 2.
*/
void Shape::line_fillx(D_real a, D_real b, std::vector<Solid_Node> &node, D_real t)
{
	bool_enclosed = false;

	if ((C_xb / static_cast<D_real>(numb_nodes)) > (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))))
	{
		std::stringstream warning;
		warning << "the distance between two solid points of the line is greater than the grid space at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}


    D_real dx = C_xb / static_cast<D_real>(numb_nodes);
	node.at(0).x = 0. + shape_offset_x0;
	node.at(0).y = b + shape_offset_y0;

#if (C_FSI_INTERFACE == 1)
	D_real arc = dx * sqrt(SQ(a) + 1);
	node.at(0).u = 0;
	node.at(0).v = 0;
	node.at(0).area = arc;
#endif

	D_uint end_i = 0;
	for (D_uint i = 1; i < numb_nodes; ++i)
	{
		node.at(i).x = node.at(i - 1).x + dx + shape_offset_x0;
		node.at(i).y = a * node.at(i).x + b + shape_offset_y0;
		if ((node.at(i).y) < 0)
		{
			node.at(i).x = -b/ a + shape_offset_x0;
			node.at(i).y = 0 + shape_offset_y0;
			end_i = i + 1;
			break;
		}
#if (C_FSI_INTERFACE == 1)
		node.at(i).u = 0;
		node.at(i).v = 0;
		node.at(i).area = arc;
#endif
	}

	if (end_i > 0)
	{
		for (D_uint i = end_i; i < numb_nodes; ++i)
		{
			node.at(i).x = node.at(0).x + shape_offset_x0;
			node.at(i).y = node.at(0).y + shape_offset_y0;
#if (C_FSI_INTERFACE == 1)
			node.at(i).u = 0;
			node.at(i).v = 0;
			node.at(i).area = arc;
#endif
		}
	}

#if(C_DIMS == 3)
	for (D_uint i = 0; i < numb_nodes; ++i)
	{
		node.at(i).z = (C_zb + shape_offset_z0) / 2;
	}
#endif

	if (node.at(numb_nodes - 1).y > (C_yb + shape_offset_y0))
	{
		std::stringstream warning;
		warning << "coordinate y = "<< node.at(numb_nodes - 1).y<< " of the last node exceeds the top compuational domain C_yb = "<<C_yb << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	x0 = node.at(numb_nodes / 2).x;
	y0 = node.at(numb_nodes / 2).y;

#if(C_DIMS==3)
	z0 = node.at(numb_nodes / 2).z;
#endif

}

#if (C_DIMS == 3)
/**
* @brief      function to generate 3D channel.
* @param[in]  a     constant for y = a*x + b, centerline of the channel.
* @param[in]  b     constant for y = a*x + b.
* @param[in]  c     constant for z = c*x + d.
* @param[in]  d     constant for z = c*x + d.
* @param[out] node  points on the line.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::channel(D_real a, D_real b, D_real c, D_real d, D_real radius, std::vector<Solid_Node> &node, D_real t)
{

	bool_enclosed = true;

	D_real circ = 2 * C_pi * radius;  // circumference
	D_real dx = (circ + sqrt(SQ(circ) + 4 * static_cast<D_real>(numb_nodes)*C_xb *circ)) / 2 / static_cast<D_real>(numb_nodes);
	if ((static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))) < dx)
	{
		std::stringstream warning;
		warning << "the distance (" << dx << ") between two adjacent solid points of the channel is greater than the grid space(" << (static_cast<D_real>(C_dx) / static_cast<D_real>(two_power_n(C_max_level))) << ") at the finest refinement level" << std::endl;
		log_warning(warning.str(), Log_function::logfile);
	}

	D_uint nx = static_cast<D_uint>(C_xb / dx + C_eps) + 1;
	D_uint ncirc = static_cast<D_uint>(circ / dx + C_eps);
	D_real xtemp0, ytemp0, ytemp1;
	D_real ztemp0, ztemp1;
	D_real theta;
	D_uint icount = 0;
	
	for (unsigned int ix = 0; ix < nx; ++ix)
	{
		xtemp0 = static_cast<D_real>(ix) * dx;
		ytemp0 = a * xtemp0 + b;
		ztemp0 = c * xtemp0 + d;
		for (unsigned int icirc = 0; icirc < ncirc; ++icirc)
		{
			theta = 2 * C_pi*static_cast<D_real>(icirc)/(ncirc + 1);
			ytemp1 = radius * cos(2 * theta);
			ztemp1 = radius * sin(2 * theta);
			node.at(icount).x = xtemp0 + shape_offset_x0;
			node.at(icount).y = ytemp0 + ytemp1 + shape_offset_y0;
			node.at(icount).z = ztemp0 + ztemp1 + shape_offset_z0;
#if (C_FSI_INTERFACE == 1)
			node.at(icount).u = 0;
			node.at(icount).v = 0;
			node.at(icount).area = dx;
#endif
			++icount;
		}
	}
	x0 = C_xb / 2;
	y0 = a * x0 + b;
	z0 = c * x0 + d;

	if (icount < numb_nodes)
	{
		node.resize(icount);
		std::stringstream warning;
		warning << "the number of nodes is not a multiplier of number of segments in x direction, set numb_nodes = "<< numb_nodes <<" as "<<icount << std::endl;
		log_warning(warning.str(), Log_function::logfile);
		numb_nodes = icount;
	}

}

 #endif

/**
* @brief      function to read geometry data (cloud point).
* @param[in]  x0     geometry center, must inside the geometry for fool fill method.
* @param[in]  y0     geometry center.
* @param[in]  z0     geometry center.
* @param[out] node   infomration of solid points.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::geofile(D_real xc, D_real yc, D_real zc, std::vector<Solid_Node> &node, D_real t)
{
	bool_enclosed = true;
	D_real x_offset = shape_offset_x0, y_offset = shape_offset_y0;
#if (C_DIMS == 3)
	D_real z_offset = shape_offset_z0;
#endif

	std::istringstream istr;
	std::string s;
	// std::ifstream file_in("./stl/airplane.txt", std::ios::in);
	std::ifstream file_in(F_model_path, std::ios::in);
	if (!file_in.is_open())
	{
		std::stringstream error;
		char* buffer;
		// if ((buffer = _getcwd(NULL, 0)) == NULL)
		if ((buffer = getcwd(NULL, 0)) == NULL)
		{
			perror("getcwd error");
		}
		else
		{
			printf("%s\n", buffer);
			free(buffer);
		}
		error << "Can't open geometry file " << std::endl;
		log_error(error.str(), Log_function::logfile);
	}

	// check number of vertices and resise vector (nodes)
	numb_nodes = std::count(std::istreambuf_iterator<char>(file_in), std::istreambuf_iterator<char>(), '\n');;
	node.resize(numb_nodes);
	file_in.seekg(0, std::ios::beg);
	//file_in.seekg(0, std::ios::end);
	//std::streampos fp = file_in.tellg();
	//if (int(fp) == 0)
	//{
	//	std::cout << "test" << std::endl;
	//}

	D_real x_min = 0, y_min = 0;
#if (C_DIMS == 3)
	D_real z_min = 0;
#endif
	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		std::getline(file_in, s); istr.str(s);

		istr >> node.at(i).x >> node.at(i).y;
#if (C_DIMS == 3)
		istr >> node.at(i).z;
#endif
		if (node.at(i).x < x_min)
		{
			x_min = node.at(i).x;
		}

		if (node.at(i).y < y_min)
		{
			y_min = node.at(i).y;
		}
#if (C_DIMS == 3)
		if (node.at(i).z < z_min)
		{
			z_min = node.at(i).z;
		}
#endif
		istr.clear();
	}

	x_offset += xc;
	y_offset += yc;
#if (C_DIMS == 3)
	z_offset += zc;
#endif

	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x += x_offset;
		node.at(i).y += y_offset;
#if (C_DIMS == 3)
		node.at(i).z += z_offset;
#endif
	}

	x0 = x_offset;
	y0 = y_offset;
#if (C_DIMS == 3)
	z0 = z_offset;
#endif

	file_in.close();

	//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
		//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
}

/**
* @brief      function to read geometry data (STL format).
* @param[in]  x0     geometry center, must inside the geometry for fool fill method.
* @param[in]  y0     geometry center.
* @param[in]  z0     geometry center.
* @param[in]  t     time.
* @note       Supports arbitrary solid positioning within the computational domain.
*/
void Shape::geofile_stl(D_real xc, D_real yc, D_real zc, D_real t)
{
	std::stringstream error;

	bool_enclosed = false;
	D_real x_offset = shape_offset_x0, y_offset = shape_offset_y0;
#if (C_DIMS == 3)
	D_real z_offset = shape_offset_z0;
#endif

	std::istringstream istr;
	std::string s;
	std::ifstream file_in(F_model_path, std::ios::in);
	if (!file_in.is_open())
	{
		std::stringstream error;
		char* buffer;
		// if ((buffer = _getcwd(NULL, 0)) == NULL)
		if ((buffer = getcwd(NULL, 0)) == NULL)
		{
			perror("getcwd error");
		}
		else
		{
			printf("%s\n", buffer);
			free(buffer);
		}
		error << "Can't open geometry file " << std::endl;
		log_error(error.str(), Log_function::logfile);
	}

	D_real x_min = 0, y_min = 0;
#if (C_DIMS == 3)
	D_real z_min = 0;
#endif
	
	while (std::getline(file_in, s))
	{
		std::string str_split,str_split2;
		istr.str(s); istr >> str_split;
		if ("facet" == str_split)
		{
			Solid_Face triFace_tmp;
			istr >> str_split; // "normal"
			if (str_split != "normal") {  
				error << "Not read 'facet normal' but: " << str_split << std::endl; 
				log_error(error.str(), Log_function::logfile);
			}
			istr >> triFace_tmp.faceNorm.at(0) >> triFace_tmp.faceNorm.at(1);
			#if (C_DIMS == 3)
			istr >> triFace_tmp.faceNorm.at(2);
			#endif

			if (fabs(triFace_tmp.faceNorm.at(0)) < C_eps && fabs(triFace_tmp.faceNorm.at(1)) < C_eps)
				triFace_tmp.faceDir_type = z_plane;
			else if (fabs(triFace_tmp.faceNorm.at(0)) < C_eps && fabs(triFace_tmp.faceNorm.at(2)) < C_eps)
				triFace_tmp.faceDir_type = y_plane;
			else if (fabs(triFace_tmp.faceNorm.at(1)) < C_eps && fabs(triFace_tmp.faceNorm.at(2)) < C_eps)
				triFace_tmp.faceDir_type = x_plane;
			else
				triFace_tmp.faceDir_type = non_orth;
			

			std::getline(file_in, s); istr.str(""); istr.clear(); istr.str(s);
			istr >> str_split >> str_split2;
			if (str_split != "outer" || str_split2 != "loop") {
				error << "Not read 'facet normal' but: " << str_split << std::endl;
				log_error(error.str(), Log_function::logfile);
			}

			// Read three verteies of the face triangle
			{
				std::getline(file_in,s); istr.str(""); istr.clear(); istr.str(s); 
				istr >> str_split;
				if (str_split != "vertex") {
					error << "Not read 'vortex' but: " << str_split << std::endl;
					log_error(error.str(), Log_function::logfile);
				}
				istr >> triFace_tmp.vertex1.x >> triFace_tmp.vertex1.y;
				#if (C_DIMS == 3)
				istr >> triFace_tmp.vertex1.z;
				#endif
			}
			{
				std::getline(file_in,s); istr.str(""); istr.clear(); istr.str(s);
				istr >> str_split;
				if (str_split != "vertex") {
					error << "Not read 'vortex' but: " << str_split << std::endl;
					log_error(error.str(), Log_function::logfile);
				}
				istr >> triFace_tmp.vertex2.x >> triFace_tmp.vertex2.y;
				#if (C_DIMS == 3)
				istr >> triFace_tmp.vertex2.z;
				#endif
			}
			{
				std::getline(file_in,s); istr.str(""); istr.clear(); istr.str(s);
				istr >> str_split;
				if (str_split != "vertex") {
					error << "Not read 'vortex' but: " << str_split << std::endl;
					log_error(error.str(), Log_function::logfile);
				}
				istr >> triFace_tmp.vertex3.x >> triFace_tmp.vertex3.y;
				#if (C_DIMS == 3)
				istr >> triFace_tmp.vertex3.z;
				#endif
			}
			std::getline(file_in,s); istr.str(""); istr.clear(); istr.str(s); 
			istr >> str_split;
			if (str_split != "endloop") {
					error << "Not read 'endloop' but: " << str_split << std::endl;
					log_error(error.str(), Log_function::logfile);
			}
			std::getline(file_in,s); istr.str(""); istr.clear(); istr.str(s); 
			istr >> str_split;
			if (str_split != "endfacet") {
					error << "Not read 'endfacet' but: " << str_split << std::endl;
					log_error(error.str(), Log_function::logfile);
			}
			triFace.push_back(triFace_tmp);
		}
		istr.str(""); istr.clear();
	}

	// Sample points on the triangle surface to add the solid points
	sample_solidPoints_onTriFace();

	numb_nodes = node.size();

	printf("------ Read STL File and Sample Nodes ------\n");
	std::cout << "  Num of Node for construct AMR: " << numb_nodes << std::endl;
	std::cout << "  Num of triFace: " << triFace.size() << std::endl;
	printf("------ ------------------------------ ------\n");


	for (Solid_Node inode : node) {
		if (inode.x < x_min)
			x_min = inode.x;
		if (inode.y < y_min)
			y_min = inode.y;
		if (inode.z < z_min)
			z_min = inode.z;
	}

	x_offset += xc;
	y_offset += yc;
#if (C_DIMS == 3)
	z_offset += zc;
#endif

	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x += x_offset;
		node.at(i).y += y_offset;
#if (C_DIMS == 3)
		node.at(i).z += z_offset;
#endif
	}

	for (std::vector<Solid_Face>::iterator iTri = triFace.begin(); iTri != triFace.end(); ++iTri) {
		iTri->vertex1.x += x_offset; iTri->vertex1.y += y_offset; iTri->vertex1.z += z_offset;
		iTri->vertex2.x += x_offset; iTri->vertex2.y += y_offset; iTri->vertex2.z += z_offset;
		iTri->vertex3.x += x_offset; iTri->vertex3.y += y_offset; iTri->vertex3.z += z_offset;
	}

		// Validate all triangle vertices are within domain bounds
	for (std::vector<Solid_Face>::const_iterator iTri = triFace.cbegin(); iTri != triFace.cend(); ++iTri) {
		const Solid_Node* verts[3] = {&iTri->vertex1, &iTri->vertex2, &iTri->vertex3};
		for (int iv = 0; iv < 3; ++iv) {
			const Solid_Node& v = *verts[iv];
			if (v.x < (C_domain[0] - C_eps) || v.x > (C_domain[3] + C_eps) ||
			    v.y < (C_domain[1] - C_eps) || v.y > (C_domain[4] + C_eps) ||
			    v.z < (C_domain[2] - C_eps) || v.z > (C_domain[5] + C_eps))
			{
				std::stringstream warning;
				warning << "STL vertex (" << v.x << ", " << v.y << ", " << v.z
				        << ") after origin placement is outside C_domain ["
				        << C_domain[0] << "," << C_domain[3] << "] x ["
				        << C_domain[1] << "," << C_domain[4] << "] x ["
				        << C_domain[2] << "," << C_domain[5]
				        << "]. Check C_solid_origin and STL model size." << std::endl;
				log_error(warning.str(), Log_function::logfile);
			}
		}
	}

	x0 = x_offset;
	y0 = y_offset;
#if (C_DIMS == 3)
	z0 = z_offset;
#endif

	std::cout << "(x0, y0, z0) " << x0 << " " << y0 << " " << z0 << std::endl;
	// std::cout << "------------------------------------------------------\n" << std::endl;

	file_in.close();

	//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
		//std::cout <<"center of the geometry: " << xc << ", " << yc << ", " << zc << std::endl;
}


/**
* @brief      function to update information of solid nodes.
* @param[in] node  solid points.
* @param[in]  t     time.
* @note       it is assumed that the range of x is from 0 to C_xb. Both y and z rely on x. Points are distributed with equal distance.
*/
void Shape::geofile(std::vector<Solid_Node> &node, D_real t)
{
	x0+= C_dx / two_power_n(C_max_level) * 0.999;
	for (unsigned int i = 0; i < numb_nodes; ++i)
	{
		node.at(i).x += C_dx / two_power_n(C_max_level) * 0.999;
		node.at(i).y += 0;
#if (C_DIMS == 3)
		node.at(i).z += 0;
#endif
	}
}

/**
 * @brief Sample points on the triangle surface to add the solid points
 * @date 2023/5/31
 */
void Shape::sample_solidPoints_onTriFace()
{
	D_real sample_dx = C_dx / two_power_n(C_max_level) / 10.;

	// Build the cube inscribing triangle
	// for (Solid_Face iTriFace : triFace)
	for (auto iTriFace = triFace.begin(); iTriFace != triFace.end(); ++iTriFace)
	{
		std::vector<Solid_Node> node_onSurface;
		std::vector<unsigned int> node_on_which_Surface;
		unsigned int itriFace_idx = std::distance(triFace.begin(), iTriFace);
		if (itriFace_idx >= triFace.size()) 
			std::cout << "itriFace_idx " << itriFace_idx << std::endl;

		switch (iTriFace->faceDir_type)
		{
		case x_plane:
		{
			D_real pane1_y = min_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);
			D_real pane1_z = min_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);

			D_real pane2_y = max_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);
			D_real pane2_z = max_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);

			Solid_Node node_temp;
			// node_temp.dir = iTriFace->faceNorm;

			for (int j_seg = 1; j_seg < (int)((pane2_y-pane1_y)/sample_dx); ++j_seg)
			{
				D_real y_cut = j_seg * sample_dx + pane1_y;
				for (int k_seg = 1; k_seg < (int)((pane2_z-pane1_z)/sample_dx); ++k_seg)
				{
					D_real z_cut = k_seg * sample_dx + pane1_z;
					if (intersect_line_with_triangle(D_vec(iTriFace->vertex1.x+1, y_cut, z_cut), D_vec(iTriFace->vertex1.x, y_cut, z_cut), *iTriFace, node_temp)) {
						node_onSurface.push_back(node_temp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
			}
		}
			break;

		case y_plane:
		{
			D_real pane1_x = min_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real pane1_z = min_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);

			D_real pane2_x = max_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real pane2_z = max_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);

			Solid_Node node_temp;
			// node_temp.dir = iTriFace->faceNorm;

			for (int i_seg = 1; i_seg < (int)((pane2_x-pane1_x)/sample_dx); ++i_seg)
			{
				D_real x_cut = i_seg * sample_dx + pane1_x;
				for (int k_seg = 1; k_seg < (int)((pane2_z-pane1_z)/sample_dx); ++k_seg)
				{
					D_real z_cut = k_seg * sample_dx + pane1_z;
					if (intersect_line_with_triangle(D_vec(x_cut, iTriFace->vertex1.y+1, z_cut), D_vec(x_cut, iTriFace->vertex1.y, z_cut), *iTriFace, node_temp)) {
						node_onSurface.push_back(node_temp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
			}
		}
			break;

		case z_plane:
		{
			D_real pane1_x = min_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real pane1_y = min_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);

			D_real pane2_x = max_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real pane2_y = max_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);

			Solid_Node node_temp;
			// node_temp.dir = iTriFace->faceNorm;

			for (int i_seg = 1; i_seg < (int)((pane2_x-pane1_x)/sample_dx); ++i_seg)
			{
				D_real x_cut = i_seg * sample_dx + pane1_x;
				for (int j_seg = 1; j_seg < (int)((pane2_y-pane1_y)/sample_dx); ++j_seg)
				{
					D_real y_cut = j_seg * sample_dx + pane1_y;
					if (intersect_line_with_triangle(D_vec(x_cut, y_cut, iTriFace->vertex1.z+1), D_vec(x_cut, y_cut, iTriFace->vertex1.z), *iTriFace, node_temp)) {
						node_onSurface.push_back(node_temp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
			}
		}
			break;

		case non_orth:
		{
			D_real cube1_x = min_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real cube1_y = min_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);
			D_real cube1_z = min_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);
			
			D_real cube2_x = max_of_three(iTriFace->vertex1.x, iTriFace->vertex2.x, iTriFace->vertex3.x);
			D_real cube2_y = max_of_three(iTriFace->vertex1.y, iTriFace->vertex2.y, iTriFace->vertex3.y);
			D_real cube2_z = max_of_three(iTriFace->vertex1.z, iTriFace->vertex2.z, iTriFace->vertex3.z);

			Solid_Node node_tmp;
			// node_tmp.dir = iTriFace->faceNorm;

			for (int i_seg = 1; i_seg < (int)((cube2_x - cube1_x) / sample_dx); ++i_seg)
			{
				D_real x_cut = i_seg * sample_dx + cube1_x;
				for (int j_seg = 1; j_seg < (int)((cube2_y - cube1_y) / sample_dx); ++j_seg)
				{
					D_real y_cut = j_seg * sample_dx + cube1_y;
					if (intersect_line_with_triangle(D_vec(x_cut,y_cut,cube1_z), D_vec(x_cut,y_cut,cube2_z), *iTriFace, node_tmp)) {
						node_onSurface.push_back(node_tmp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
				for (int k_seg = 1; k_seg < (int)((cube2_z - cube1_z) / sample_dx); ++k_seg)
				{
					D_real z_cut = k_seg * sample_dx + cube1_z;
					if (intersect_line_with_triangle(D_vec(x_cut,cube1_y,z_cut), D_vec(x_cut,cube2_y,z_cut), *iTriFace, node_tmp)) {
						node_onSurface.push_back(node_tmp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
			}
			for (int j_seg = 1; j_seg < (int)((cube2_y - cube1_y) / sample_dx); ++j_seg)
			{
				D_real y_cut = j_seg * sample_dx + cube1_y;
				for (int k_seg = 1; k_seg < (int)((cube2_z - cube1_z) / sample_dx); ++k_seg)
				{
					D_real z_cut = k_seg * sample_dx + cube1_z;
					if (intersect_line_with_triangle(D_vec(cube1_x,y_cut,z_cut), D_vec(cube2_x,y_cut,z_cut), *iTriFace, node_tmp)) {
						node_onSurface.push_back(node_tmp);
						node_on_which_Surface.push_back(itriFace_idx);
					}
				}
			}
		}
			break;
		
		default:
			unreachable();
			break;
		}

		node.push_back(iTriFace->vertex1); node.push_back(iTriFace->vertex2); node.push_back(iTriFace->vertex3);
		// triFaceIdx_of_node.push_back(itriFace_idx); triFaceIdx_of_node.push_back(itriFace_idx); triFaceIdx_of_node.push_back(itriFace_idx);
		node.insert(node.end(), node_onSurface.begin(), node_onSurface.end());
		// triFaceIdx_of_node.insert(triFaceIdx_of_node.end(), node_on_which_Surface.begin(), node_on_which_Surface.end());
	}
}

void Shape::print_pointCloud()
{

	std::ofstream solid_nodes_file;
	solid_nodes_file.open("solid_pointCloude.vtk", std::ios::out);
	solid_nodes_file << "# vtk DataFile Version 2.0\n\n";
	solid_nodes_file << "ASCII\nDATASET POLYDATA\n";
	solid_nodes_file << "POINTS " << numb_nodes << " float\n";

	for (unsigned int i = 0; i < numb_nodes; ++i) {
		solid_nodes_file << node.at(i).x << " " << node.at(i).y  << " " << node.at(i).z << "\n";
	}
	solid_nodes_file << "VERTICES " << numb_nodes << " " << 2*numb_nodes << "\n";
	for (unsigned int i = 0; i < numb_nodes; ++i) {
		solid_nodes_file << "1 " << i << "\n";
	}
	solid_nodes_file.close();
}

/**
 * @brief Get the intersection point by Moller-Trumbore algorithm
 * @date 2023-5-31
 * @param[in] p1 
 * @param[in] p2 
 * @param[in] triFace 
 * @param[out] intersectPoint 
 * @return true Ray really intersects with triangle face.
 * @ref https://zhuanlan.zhihu.com/p/451582864
 */
bool Shape::intersect_line_with_triangle(D_vec p1, D_vec p2, Solid_Face triFace, Solid_Node& intersectPoint)
{
	auto ray_vector = [](const Solid_Node& point_a, const Solid_Node& point_b)
	{
		return D_vec(point_b.x-point_a.x, point_b.y-point_a.y, point_b.z-point_a.z);
	};

	// std::cout << "p1 (" << p1.x << " , " << p1.y << " , " << p1.z << ") " <<
	// 			 "p2 (" << p2.x << " , " << p2.y << " , " << p2.z << ") " << std::endl;

	D_vec e1 = ray_vector(triFace.vertex1, triFace.vertex2);
	D_vec e2 = ray_vector(triFace.vertex1, triFace.vertex3);
	D_vec d = (p2 - p1).norm();
	D_vec s = D_vec(p1.x - triFace.vertex1.x, p1.y - triFace.vertex1.y, p1.z - triFace.vertex1.z);
	D_vec s1 = cross_product(d, e2);
	D_vec s2 = cross_product(s, e1);

	D_real s1e1 = dot_product(s1, e1);
	D_real t = dot_product(s2, e2)/s1e1;
	D_real b1 = dot_product(s1, s)/s1e1;
	D_real b2 = dot_product(s2, d)/s1e1;

	if (t >= C_eps && b1 >= C_eps && b2 >= C_eps && (1-b1-b2)>=C_eps)
	{
		intersectPoint.x = (1-b1-b2) * triFace.vertex1.x + b1 * triFace.vertex2.x + b2 * triFace.vertex3.x;
		intersectPoint.y = (1-b1-b2) * triFace.vertex1.y + b1 * triFace.vertex2.y + b2 * triFace.vertex3.y;
		intersectPoint.z = (1-b1-b2) * triFace.vertex1.z + b1 * triFace.vertex2.z + b2 * triFace.vertex3.z;
		return true;
	}
	return false;
}

bool Shape::intersect_line_with_triangle(D_vec p1, D_vec p2, Solid_Face triFace)
{
	auto ray_vector = [](const Solid_Node& point_a, const Solid_Node& point_b)
	{
		return D_vec(point_b.x-point_a.x, point_b.y-point_a.y, point_b.z-point_a.z);
	};

	// std::cout << "p1 (" << p1.x << " , " << p1.y << " , " << p1.z << ") " <<
	// 			 "p2 (" << p2.x << " , " << p2.y << " , " << p2.z << ") " << std::endl;

	D_vec e1 = ray_vector(triFace.vertex1, triFace.vertex2);
	D_vec e2 = ray_vector(triFace.vertex1, triFace.vertex3);
	D_vec d = (p2 - p1).norm();
	D_vec s = D_vec(p1.x - triFace.vertex1.x, p1.y - triFace.vertex1.y, p1.z - triFace.vertex1.z);
	D_vec s1 = cross_product(d, e2);
	D_vec s2 = cross_product(s, e1);

	D_real s1e1 = dot_product(s1, e1);
	D_real t = dot_product(s2, e2)/s1e1;
	D_real b1 = dot_product(s1, s)/s1e1;
	D_real b2 = dot_product(s2, d)/s1e1;

	if (t >= C_eps && b1 >= C_eps && b2 >= C_eps && (1-b1-b2)>=C_eps)
		return true;

	return false;
}

bool Shape::intersect_line_with_triangle(D_vec startPoint, D_vec p2, Solid_Face triFace, D_real& dis_to_startPoint)
{
	auto ray_vector = [](const Solid_Node& point_a, const Solid_Node& point_b)
	{
		return D_vec(point_b.x-point_a.x, point_b.y-point_a.y, point_b.z-point_a.z);
	};

	// std::cout << "p1 (" << p1.x << " , " << p1.y << " , " << p1.z << ") " <<
	// 			 "p2 (" << p2.x << " , " << p2.y << " , " << p2.z << ") " << std::endl;

	D_vec e1 = ray_vector(triFace.vertex1, triFace.vertex2);
	D_vec e2 = ray_vector(triFace.vertex1, triFace.vertex3);
	D_vec d = (p2 - startPoint).norm();
	D_vec s = D_vec(startPoint.x - triFace.vertex1.x, startPoint.y - triFace.vertex1.y, startPoint.z - triFace.vertex1.z);
	D_vec s1 = cross_product(d, e2);
	D_vec s2 = cross_product(s, e1);

	D_real s1e1 = dot_product(s1, e1);
	D_real t = dot_product(s2, e2)/s1e1;
	D_real b1 = dot_product(s1, s)/s1e1;
	D_real b2 = dot_product(s2, d)/s1e1;

	D_vec intersectPoint;

	if (t >= C_eps && b1 >= C_eps && b2 >= C_eps && (1-b1-b2)>=C_eps)
	{
		intersectPoint.x = (1-b1-b2) * triFace.vertex1.x + b1 * triFace.vertex2.x + b2 * triFace.vertex3.x;
		intersectPoint.y = (1-b1-b2) * triFace.vertex1.y + b1 * triFace.vertex2.y + b2 * triFace.vertex3.y;
		intersectPoint.z = (1-b1-b2) * triFace.vertex1.z + b1 * triFace.vertex2.z + b2 * triFace.vertex3.z;
		
		dis_to_startPoint = two_points_length(startPoint, intersectPoint);
		return true;
	}
	return false;
}

bool Shape::judge_point_within_Shape(D_vec point, Solid_Face triFace, D_vec &intersectPoint)
{
	Solid_Node iscp_temp;
	if (!intersect_line_with_triangle(point, 
		 							  point+D_vec(triFace.faceNorm.at(0), triFace.faceNorm.at(1), triFace.faceNorm.at(2)),
									  triFace,
									  iscp_temp))
		return false;
	intersectPoint.x = iscp_temp.x;
	intersectPoint.y = iscp_temp.y;
	intersectPoint.z = iscp_temp.z;

	D_vec iscp2point = (point - intersectPoint).norm();
	D_vec triFace_norm = D_vec(triFace.faceNorm.at(0), triFace.faceNorm.at(1), triFace.faceNorm.at(2)).norm();
	if (dot_product(iscp2point, triFace_norm) == 1.)
		return true;
	else if (dot_product(iscp2point, triFace_norm) == -1.)
		return false;
	else
		unreachable();
}

bool Shape::judge_point_within_Shape(D_vec point, Solid_Face triFace)
{
	Solid_Node iscp_temp;

	if (!intersect_line_with_triangle(point, 
		 							  point+D_vec(triFace.faceNorm.at(0), triFace.faceNorm.at(1), triFace.faceNorm.at(2)),
									  triFace,
									  iscp_temp))
		return false;

	D_vec iscp2point = (point - D_vec(iscp_temp.x,iscp_temp.y,iscp_temp.z)).norm();
	D_vec triFace_norm = D_vec(triFace.faceNorm.at(0), triFace.faceNorm.at(1), triFace.faceNorm.at(2)).norm();
	if (dot_product(iscp2point, triFace_norm) == 1.)
		return true;
	else if (dot_product(iscp2point, triFace_norm) == -1.)
		return false;
	else
		unreachable();
}

D_real Shape::two_points_length(D_vec a, D_vec b)
{
	return (D_real)sqrtf64((b.x-a.x)*(b.x-a.x) + (b.y-a.y)*(b.y-a.y) + (b.z-a.z)*(b.z-a.z));
}

D_real Shape::point_triangle_distance(D_vec p, Solid_Face triFace)
{
	// Ericson, Real-Time Collision Detection: ClosestPtPointTriangle.
	// Returns the shortest distance from p to the triangle (a,b,c).
	D_vec a(triFace.vertex1.x, triFace.vertex1.y, triFace.vertex1.z);
	D_vec b(triFace.vertex2.x, triFace.vertex2.y, triFace.vertex2.z);
	D_vec c(triFace.vertex3.x, triFace.vertex3.y, triFace.vertex3.z);

	D_vec ab = b - a;
	D_vec ac = c - a;
	D_vec ap = p - a;

	D_real d1 = dot_product(ab, ap);
	D_real d2 = dot_product(ac, ap);
	if (d1 <= 0.0 && d2 <= 0.0) return two_points_length(p, a);  // vertex region A

	D_vec bp = p - b;
	D_real d3 = dot_product(ab, bp);
	D_real d4 = dot_product(ac, bp);
	if (d3 >= 0.0 && d4 <= d3) return two_points_length(p, b);    // vertex region B

	D_real vc = d1 * d4 - d3 * d2;
	if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
		D_real v = d1 / (d1 - d3);
		D_vec proj(a.x + ab.x * v, a.y + ab.y * v, a.z + ab.z * v);
		return two_points_length(p, proj);  // edge AB
	}

	D_vec cp = p - c;
	D_real d5 = dot_product(ab, cp);
	D_real d6 = dot_product(ac, cp);
	if (d6 >= 0.0 && d5 <= d6) return two_points_length(p, c);    // vertex region C

	D_real vb = d5 * d2 - d1 * d6;
	if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
		D_real w = d2 / (d2 - d6);
		D_vec proj(a.x + ac.x * w, a.y + ac.y * w, a.z + ac.z * w);
		return two_points_length(p, proj);  // edge AC
	}

	D_real va = d3 * d6 - d5 * d4;
	if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
		D_real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		D_vec proj(b.x + (c.x - b.x) * w, b.y + (c.y - b.y) * w, b.z + (c.z - b.z) * w);
		return two_points_length(p, proj);  // edge BC
	}

	// face region — project p onto the triangle plane
	D_real denom = 1.0 / (va + vb + vc);
	D_real v = vb * denom;
	D_real w = vc * denom;
	D_vec proj(a.x + ab.x * v + ac.x * w,
	           a.y + ab.y * v + ac.y * w,
	           a.z + ab.z * v + ac.z * w);
	return two_points_length(p, proj);
}

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;

int planeBoxOverlap(double normal[3], double vert[3], double maxbox[3])	// -NJMP-
{
  int q;
  double vmin[3],vmax[3],v;
  for(q=X;q<=Z;q++)
  {
    v=vert[q];					// -NJMP-
    if(normal[q]>0.0f)
    {
      vmin[q]=-maxbox[q] - v;	// -NJMP-
      vmax[q]= maxbox[q] - v;	// -NJMP-
    }
    else
    {
      vmin[q]= maxbox[q] - v;	// -NJMP-
      vmax[q]=-maxbox[q] - v;	// -NJMP-
    }
  }
  if(DOT(normal,vmin)>0.0f) return 0;	// -NJMP-
  if(DOT(normal,vmax)>=0.0f) return 1;	// -NJMP-

  return 0;
}

/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			       	   \
	p2 = a*v2[Y] - b*v2[Z];			       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_X2(a, b, fa, fb)			   \
	p0 = a*v0[Y] - b*v0[Z];			           \
	p1 = a*v1[Y] - b*v1[Z];			       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return false;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p2 = -a*v2[X] + b*v2[Z];	       	       	   \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Y1(a, b, fa, fb)			   \
	p0 = -a*v0[X] + b*v0[Z];		      	   \
	p1 = -a*v1[X] + b*v1[Z];	     	       	   \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
	if(min>rad || max<-rad) return false;

/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return false;

#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				   \
	p1 = a*v1[X] - b*v1[Y];			           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return false;

bool Shape::triBoxOverlap(double boxcenter[3],double boxhalfsize[3],double triverts[3][3])
{
  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */
   double v0[3],v1[3],v2[3];
//   double axis[3];
   double min,max,p0,p1,p2,rad,fex,fey,fez;		// -NJMP- "d" local variable removed
   double normal[3],e0[3],e1[3],e2[3];

   /* This is the fastest branch on Sun */
   /* move everything so that the boxcenter is in (0,0,0) */
   SUB(v0,triverts[0],boxcenter);
   SUB(v1,triverts[1],boxcenter);
   SUB(v2,triverts[2],boxcenter);

   /* compute triangle edges */
   SUB(e0,v1,v0);      /* tri edge 0 */
   SUB(e1,v2,v1);      /* tri edge 1 */
   SUB(e2,v0,v2);      /* tri edge 2 */

   /* Bullet 3:  */
   /*  test the 9 tests first (this was faster) */
   fex = fabsf(e0[X]);
   fey = fabsf(e0[Y]);
   fez = fabsf(e0[Z]);
   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
   AXISTEST_Z12(e0[Y], e0[X], fey, fex);

   fex = fabsf(e1[X]);
   fey = fabsf(e1[Y]);
   fez = fabsf(e1[Z]);
   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
   AXISTEST_Z0(e1[Y], e1[X], fey, fex);

   fex = fabsf(e2[X]);
   fey = fabsf(e2[Y]);
   fez = fabsf(e2[Z]);
   AXISTEST_X2(e2[Z], e2[Y], fez, fey);
   AXISTEST_Y1(e2[Z], e2[X], fez, fex);
   AXISTEST_Z12(e2[Y], e2[X], fey, fex);

   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */

   /* test in X-direction */
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return false;

   /* test in Y-direction */
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return false;

   /* test in Z-direction */
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return false;

   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */
   CROSS(normal,e0,e1);
   // -NJMP- (line removed here)
   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return false;	// -NJMP-

   return true;   /* box and triangle overlaps */
}

Solid_Face Solid_Face::triFace_offset(D_vec dir, D_real offset)
{
	Solid_Face triFace1 = *this;

	D_vec norm = dir.norm();
	D_vec offset_norm = {offset * norm.x, offset * norm.y, offset * norm.z};

	triFace1.vertex1.x += offset_norm.x;
	triFace1.vertex1.y += offset_norm.y;
	triFace1.vertex1.z += offset_norm.z;

	triFace1.vertex2.x += offset_norm.x;
	triFace1.vertex2.y += offset_norm.y;
	triFace1.vertex2.z += offset_norm.z;

	triFace1.vertex3.x += offset_norm.x;
	triFace1.vertex3.y += offset_norm.y;
	triFace1.vertex3.z += offset_norm.z;

	return triFace1;
}

VertexBox Solid_Face::create_triFace_AABB()
{
	D_vec leftSouthBot(min_of_three(this->vertex1.x, this->vertex2.x, this->vertex3.x),
	                   min_of_three(this->vertex1.y, this->vertex2.y, this->vertex3.y),
	                   min_of_three(this->vertex1.z, this->vertex2.z, this->vertex3.z));

	D_vec rightNorthTop(max_of_three(this->vertex1.x, this->vertex2.x, this->vertex3.x),
	                    max_of_three(this->vertex1.y, this->vertex2.y, this->vertex3.y),
	                    max_of_three(this->vertex1.z, this->vertex2.z, this->vertex3.z));

	// std::cout << "leftSouthBot (" << leftSouthBot.x << ", " << leftSouthBot.y << ", " << leftSouthBot.z << ")" << std::endl;
	// std::cout << "rightNorthTop (" << rightNorthTop.x << ", " << rightNorthTop.y << ", " << rightNorthTop.z << ")" << std::endl;

	D_real ddx = C_dx / two_power_n(C_max_level);

		// Expand AABB by 1 cell in each direction to capture all surface-facing cells
	// box_p1 = floor(max/ddx) + 2. The previous code used +1 on the high
	// side, which made the halo asymmetric (-1 low, 0 high) and skipped every
	// cell on the solid-interior side of each triangle — that's why
	// nSameDirFalse was 0 and no SURFACE cells were ever produced.
	D_int box_p0_x = static_cast<D_int>(leftSouthBot.x / ddx + C_eps) - 1;
	D_int box_p0_y = static_cast<D_int>(leftSouthBot.y / ddx + C_eps) - 1;
	D_int box_p0_z = static_cast<D_int>(leftSouthBot.z / ddx + C_eps) - 1;

	D_int box_p1_x = static_cast<D_int>(rightNorthTop.x / ddx + C_eps) + 2;
	D_int box_p1_y = static_cast<D_int>(rightNorthTop.y / ddx + C_eps) + 2;
	D_int box_p1_z = static_cast<D_int>(rightNorthTop.z / ddx + C_eps) + 2;

	D_morton leftSouthBot_code = Morton_Assist::pointer_me->morton_encode(box_p0_x, box_p0_y, box_p0_z);
	D_morton rightNorthTop_code = Morton_Assist::pointer_me->morton_encode(box_p1_x, box_p1_y, box_p1_z);

	// double xxx,yyy,zzz;
	// Morton_Assist::pointer_me->compute_coordinate(leftSouthBot_code, C_max_level, xxx, yyy, zzz);
	// std::cout << "leftSouthBot coordinate  xxx " << xxx << " yyy " << yyy << " zzz " << zzz << std::endl;
	// Morton_Assist::pointer_me->compute_coordinate(rightNorthTop_code, C_max_level, xxx, yyy, zzz);
	// std::cout << "rightNorthTop coordinate  xxx " << xxx << " yyy " << yyy << " zzz " << zzz << std::endl;

	return VertexBox(leftSouthBot, rightNorthTop, 
					 std::array<std::array<int,3>,2>{{{box_p0_x, box_p0_y, box_p0_z},{box_p1_x, box_p1_y, box_p1_z}}},
					 std::array<D_morton,2>{leftSouthBot_code, rightNorthTop_code});

}