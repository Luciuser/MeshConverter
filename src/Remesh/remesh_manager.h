#ifndef API_REMSHMANAGER_H
#define API_REMSHMANAGER_H

#include <functional>
#include <vector>

#include "meshIO.h"
#include "remesh_parameter.h"

#include <Eigen/Dense>

#include "Halfedge/AABB_Tree.h"
#include "Halfedge/PolyMesh.h"
#include "Halfedge/PolyMesh_Base.h"

namespace MESHIO
{  
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
	//     otherwise fail;
	int checkInput(const std::vector<double>& points, const std::vector<int>& triangles, const std::vector<int>& surfaceID);
	
	// Get size function from size vtk file
	// Input: 
	//     filename: the size file in vtk format.
	// Output: 
	//     size_function: the size function from size file
	// Return:
	//     0  if success;
	//     otherwise fail;
	int getSizeFunction(const std::string filename, std::function<double(double, double, double)>& size_function);

	// initial parameter from mesh
	// Input:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     otherwise fail;
	int initial(const Mesh& mesh);

	// initial global/user target length from mesh
	// Input:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     otherwise fail;
	int initialUserTargetLength(const Mesh& mesh);

	int addSizeFunction(std::function<double(double, double, double)> size_function);
	int deleteSizeFunction();

	// get mesh data from PolyMesh class
	// Output:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     otherwise fail;
	int getMesh(Mesh& out);

	// repir mesh
	// Output:
	//     mesh: point coordinates and triangular topology
	// Return:
	//     0  if success;
	//     otherwise fail;
	int repair(Mesh& mesh, double eps = 1e-6);

	// build AABB tree for check swap operation valid
	// Input:
	//     mesh: point coordinates and triangular topology
	// Output:
	//     tree: AABB tree
	// Return:
	//     0  if success;
	//     otherwise fail;
	int buildAABBTree(Mesh& mesh, aabb::Tree **tree);

	// build AABB tree for check swap operation valid
	// Input:
	//     mesh: point coordinates and triangular topology
	// Output:
	//     tree: AABB tree
	// Return:
	//     0  if success;
	//     otherwise fail;
	int buildAABBTree(Mesh& mesh);

	// update aabb tree by one new polygen
	// Input:
	//     polyface: one new polyface
	// Output:
	//     tree: AABB tree
	// Return:
	//     0  if success;
	//     otherwise fail;
	int updateAABBTree(MESHIO::polymesh::MPolyFace* polyface, aabb::Tree* tree);

	// for one element, get the AABB bounding box.
	// input:
	//     p1/p2/p3: the point coordinates of a trangle
	// output:
	//     lower: the lower position of AABB bounding box
	//     upper: the upper position of AABB bounding box
	// Return:
	//     0  if success;
	//     otherwise fail;
	int calculateBoundingBoxForOneElement(const std::vector<Eigen::Vector3d>& point, std::vector<double>& lower, std::vector<double>& upper);

	// for two triangles, check if exist intersection
	// input:
	//     a/b/c: the point coordinates of triangle 1
	//     o/p/q: the point coordinates of triangle 2
	// Return:
	//     0  if disjoint;
	//     otherwise disjoint;
	int triangleIntersectionCheck(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& c, Eigen::Vector3d& o, Eigen::Vector3d& p, Eigen::Vector3d& q);

	// do remesh
	// Return:
	//     0  if success;
	//     otherwise fail;
	int remesh(); //  core function
	int remesh(const Mesh& mesh, Mesh& out);
	int remesh(const Mesh &mesh, std::function<double(double, double, double)>& size_function, Mesh &out);
		
	aabb::Tree* getAABBTree() { return parameter_.aabbTree(); };
	MESHIO::polymesh::MPolyFace* getPolyFace(int index) { return half_mesh_.polyface(index); };

	void setIteration(int iteration) { parameter_.setIteration(iteration); };
	int iteration() { return parameter_.iteration(); };

	void setUserTargetLength(double target_length) { parameter_.setUserTargetLength(target_length); };

private:
	// remesh operation
	void split_long_edges(MESHIO::polymesh::PolyMesh* mesh, RemeshParameter &parameter);
	void collapse_short_edges(MESHIO::polymesh::PolyMesh* mesh, RemeshParameter &parameter);
	void delete_lowdegree(MESHIO::polymesh::PolyMesh* mesh);
	void equalize_valences(MESHIO::polymesh::PolyMesh* mesh, RemeshParameter& parameter); // TODO(jin) : can only use with triangles
	void tangential_relaxation(MESHIO::polymesh::PolyMesh* mesh);
	void project_to_surface(MESHIO::polymesh::PolyMesh* mesh, AABB_Tree* abtree);
	void get_AABB_tree(MESHIO::polymesh::PolyMesh* mesh, AABB_Tree*& abtree);

	double calculateTargetEdgeLength(MESHIO::polymesh::PolyMesh* mesh);

private:
	//Mesh* mesh_ = nullptr;
	polymesh::PolyMesh half_mesh_;
	AABB_Tree* abtree_ = nullptr; // polygen mesh use
	RemeshParameter parameter_;
};

}

#endif // API_REMSHMANAGER_H