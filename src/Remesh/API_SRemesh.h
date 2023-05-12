#ifndef API_REMESHALGORITHM_H
#define API_REMESHALGORITHM_H

#include <iostream>
#include <functional>
#include <vector>

/*
* Do remesh
*/
int API_remesh_non_manifold(
	// ----- input -----//
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	const std::string size_vtk_file_name,
	// ----- output -----//
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	// ----- parameter -----//
	const bool b_useAABBtree = false,
	const bool b_split = false,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk"
);

int API_remesh_non_manifold(
	// ----- input -----//
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	const std::function<double(double, double, double)>& size_function,
	// ----- output -----//
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	// ----- parameter -----//
	const bool b_useAABBtree = false,
	const bool b_split = false,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk"
);

#endif // API_REMESHALGORITHM_H