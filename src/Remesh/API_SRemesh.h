#ifndef API_REMESHALGORITHM_H
#define API_REMESHALGORITHM_H

#include <functional>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>

#include "meshIO.h"

#include <Eigen/Dense>

#include "Halfedge/AABB_Tree.h"
#include "Halfedge/PolyMesh.h"
#include "Halfedge/PolyMesh_Base.h"

using namespace MESHIO::polymesh;

namespace MESHIO
{  

/*
* Do remesh
*/
int API_remesh_non_manifold(
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	const std::function<double(double, double, double)>& size_function,
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	const bool b_split = true,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk"
);

int API_remesh_non_manifold(
	const std::vector<double>& points,
	const std::vector<int>& triangles,
	const std::vector<int>& surfaceID,
	const std::string size_vtk_file_name,
	std::vector<double>& points_out,
	std::vector<int>& triangles_out,
	std::vector<int>& surfaceID_out,
	const bool b_split = true,
	const bool b_write_vtk = true,
	const std::string write_vtk_file_name = "remesh_debug.vtk"
);

/*
* parameters
*/
class RemeshParameter {
public:
	RemeshParameter();
	~RemeshParameter();

	// size
	double getVertexSizeOrTargetLow(const Eigen::Vector3d& vertex);
	double getVertexSizeOrTargetHigh(const Eigen::Vector3d& vertex);
	double getVertexSize(const Eigen::Vector3d& vertex);
	
	void setSizeFunction(std::function<double(double, double, double)> size_function) { size_function_ = size_function; b_use_size_function_ = true; };
	void deleteSizeFunction() { b_use_size_function_ = false; };

	double hmin() { return hmin_; }
	double hmax() { return hmax_; }
	void setHmin(double hmin) { hmin_ = hmin; };
	void setHmax(double hmax) { hmax_ = hmax; };
	void setTargetLength(double target_length) { target_length_ = target_length; };

private:
	// size
	bool b_use_size_function_ = false;
	std::function<double(double, double, double)> size_function_;
	double hmin_ = 0.00001; // the min size
	double hmax_ = 100000; // the max size

	double target_length_ = 1;
	double low_ratio_ = 4.0 / 5.0;
	double high_ratio_ = 4.0 / 3.0;
	
};

/*
* manager
*/
class RemeshManager {
public:
	RemeshManager();
	~RemeshManager();

	// Check the valid of input data. It is necessary while using this class.
	// Input: 
	//     points: the point coordinates. The length has to be n*3.
	//     triangles: the triangular topology. The length has to be m*3;
	//     surfaceID: the surface id which belong to every triangle. The length has to be m;
	// Return:
	//     0  if success;
	//     1  otherwise;
	int checkInput(const std::vector<double>& points, const std::vector<int>& triangles, const std::vector<int>& surfaceID);
	
	// Get size function from size vtk file
	// Input: 
	//     filename: the size file in vtk format.
	// Output: 
	//     size_function: the size function from size file
	// Return:
	//     0  if success;
	//     1  otherwise;
	int getSizeFunction(const std::string filename, std::function<double(double, double, double)>& size_function);

	// initial parameter from mesh
	// Input:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     1  otherwise;
	int initial(const Mesh& mesh);

	int addSizeFunction(std::function<double(double, double, double)> size_function);
	int deleteSizeFunction();

	// get mesh data from PolyMesh class
	// Output:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     1  otherwise;
	int getMesh(Mesh& out);

	// repir mesh
	// Output:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     1  otherwise;
	int repair(Mesh& mesh, double eps = 1e-6);

	int remesh(); //  core function
	int remesh(const Mesh& mesh, Mesh& out);
	int remesh(const Mesh &mesh, std::function<double(double, double, double)>& size_function, Mesh &out);
		
private:
	// remesh operation
	void split_long_edges(PolyMesh* mesh, RemeshParameter &parameter);
	void collapse_short_edges(PolyMesh* mesh, RemeshParameter &parameter);
	void delete_lowdegree(PolyMesh* mesh);
	void equalize_valences(PolyMesh* mesh);
	void tangential_relaxation(PolyMesh* mesh);
	void project_to_surface(PolyMesh* mesh, AABB_Tree* abtree);
	void get_AABB_tree(PolyMesh* mesh, AABB_Tree*& abtree);

	double calculateTargetEdgeLength(PolyMesh* mesh);

private:
	//Mesh* mesh_ = nullptr;
	polymesh::PolyMesh half_mesh_;
	AABB_Tree* abtree_ = nullptr;
	RemeshParameter parameter_;
};

}

#endif // API_REMESHALGORITHM_H