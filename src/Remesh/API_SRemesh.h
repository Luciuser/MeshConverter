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

private:
	bool b_use_size_function_;
	std::function<double(double, double, double)> size_function_;
	double hmin_;
	double hmax_;
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
	
	//int initial(Mesh& mesh);
	int addSizeFunction(std::function<double(double, double, double)> size_function);
	int deleteSizeFunction();

	int remeshNonManifold(); //  core function
	int remeshNonManifold(Mesh &mesh); //  core function
		
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
	Mesh* mesh_;
	polymesh::PolyMesh half_mesh_;
	AABB_Tree* abtree_;
	RemeshParameter parameter_;
};

}

#endif // API_REMESHALGORITHM_H