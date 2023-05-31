#ifndef API_REMESHPARAMETER_H
#define API_REMESHPARAMETER_H

#include <functional>
#include <vector>

#include "meshIO.h"

#include <Eigen/Dense>

#include "AABB.h" // use to solve intersection problem

namespace MESHIO
{  

/*
* parameters
*/
class RemeshParameter {
public:
	RemeshParameter();
	~RemeshParameter();

	// ----- AABB tree -----//

	// build AABB tree for check swap operation valid
	// Input:
	//     mesh: point coordinates and triangular topology
	// Output:
	//     tree: AABB tree
	// Return:
	//     0  if success; 
	//     otherwise fail;
	int buildAABBTree(Mesh& mesh);

	int setAABBTreeActive();
	int setAABBTreeDead();
	bool getAABBTreeActive() { return b_use_aabb_tree_; };

	// ----- size -----//
	double getVertexSizeOrTargetLow(const Eigen::Vector3d& vertex);
	double getVertexSizeOrTargetHigh(const Eigen::Vector3d& vertex);
	double getVertexSize(const Eigen::Vector3d& vertex);
	
	void setSizeFunction(std::function<double(double, double, double)> size_function) { size_function_ = size_function; b_use_size_function_ = true; };
	void deleteSizeFunction() { b_use_size_function_ = false; };

	double hmin() { return hmin_; }
	double hmax() { return hmax_; }
	aabb::Tree* aabbTree() { return aabbTree_; };
	void setHmin(double hmin) { hmin_ = hmin; };
	void setHmax(double hmax) { hmax_ = hmax; };
	void setTargetLength(double target_length) { target_length_ = target_length; };

	void setUserTargetLength(double target_length) { target_length_ = target_length; b_use_user_target_length_ = true; };
	bool getUserTargetLengthActive() { return b_use_user_target_length_; };

	// ----- iteration ----- //
	int iteration() { return iteration_; };
	void setIteration(int iteration) { iteration_ = iteration; };

private:
	int iteration_ = 10;

	// AABB tree
	bool b_use_aabb_tree_ = false;
	aabb::Tree* aabbTree_ = nullptr; // intersection use

	// size
	bool b_use_size_function_ = false;
	bool b_use_user_target_length_ = false;
	std::function<double(double, double, double)> size_function_;
	double hmin_ = 0.00001; // the min size
	double hmax_ = 100000; // the max size

	double target_length_ = 1;
	double low_ratio_ = 4.0 / 5.0;
	double high_ratio_ = 4.0 / 3.0;
};


}

#endif // API_REMESHPARAMETER_H