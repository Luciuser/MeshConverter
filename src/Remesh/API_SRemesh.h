#ifndef API_REMESHALGORITHM_H
#define API_REMESHALGORITHM_H

#include <iostream>
#include <functional>
#include <vector>

// for one mesh which is non-manifold, do remesh with sizefunciton
// input:
//     points: the point coordinates. The length has to be n*3.
//     triangles: the triangular topology. The length has to be m*3;
//     surfaceID: the surface id which belong to every triangle. The length has to be m;
// output:
//     points_out: the point coordinates.
//     triangles_out: the triangular topology.
//     surfaceID_out: the surface id which belong to every triangle.
// parameter:
//     iterations: the iter number. This API will try to cost more time to get better result with a bigger iterations value.
//     size_vtk_file_name: one size function file with point value which is written in vtk format. if size_vtk_file_name = "NO_USE", This API will remesh with average egde length.
//     size_target: try to make all edge length is size_target. if size_target <= 0, This API will remesh with average egde length.
//     b_useAABBtree: build one aabbtree to avoid the intersection of triangles while remeshing.
//     b_split: to split origin face to severals small faces. This option can remesh non-manifold input and try to save geometric features with surface IDs.
//     dulplicate_point_eps: dulplicate the result with this eps.
//     b_write_vtk: output remesh result in vtk format.
//     write_vtk_file_name: the name while outputing remesh result in vtk format. 
//     b_write_split_vtk: output all small faces remesh result in vtk format.
// Return:
//     0  if success;
//     otherwise fail;
int API_remesh_non_manifold_file(
	// ----- input -----//
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	// ----- output -----//
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	// ----- parameter -----//
	const int iterations = 10,
	const std::string size_vtk_file_name = "NO_USE",
	const double size_target = -1,
	const bool b_useAABBtree = false,
	const bool b_split = false,
	const double dulplicate_point_eps = 1e-4,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk",
	const bool b_write_split_vtk = false
);

// for one mesh which is non-manifold, do remesh with sizefunciton
// input:
//     points: the point coordinates. The length has to be n*3.
//     triangles: the triangular topology. The length has to be m*3;
//     surfaceID: the surface id which belong to every triangle. The length has to be m;
// output:
//     points_out: the point coordinates.
//     triangles_out: the triangular topology.
//     surfaceID_out: the surface id which belong to every triangle.
// parameter:
//     iterations: the iter number. This API will try to cost more time to get better result with a bigger iterations value.
//     size_function: one size function. if size_function = nullptr, This API will remesh with average egde length.
//     size_target: try to make all edge length is size_target. if size_target <= 0, This API will remesh with average egde length.
//     b_useAABBtree: build one aabbtree to avoid the intersection of triangles while remeshing.
//     b_split: to split origin face to severals small faces. This option can remesh non-manifold input and try to save geometric features with surface IDs.
//     dulplicate_point_eps: dulplicate the result with this eps.
//     b_write_vtk: output remesh result in vtk format.
//     write_vtk_file_name: the name while outputing remesh result in vtk format. 
//     b_write_split_vtk: output all small faces remesh result in vtk format.
// Return:
//     0  if success;
//     otherwise fail;
int API_remesh_non_manifold(
	// ----- input -----//
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	// ----- output -----//
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	// ----- parameter -----//
	const int iterations = 10,
	const std::function<double(double, double, double)>& size_function = nullptr,
	const double size_target = -1,
	const bool b_useAABBtree = false,
	const bool b_split = false,
	const double dulplicate_point_eps = 1e-4,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk",
	const bool b_write_split_vtk = false
);

#endif // API_REMESHALGORITHM_H