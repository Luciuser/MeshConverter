#ifndef MESHALGORITHM_H
#define MESHALGORITHM_H

#include <Eigen/Dense>
#include <map>
#include <unordered_map>
#include <string>
#include <vector>
#include "meshIO.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace MESHIO{

	bool rotatePoint(std::vector<double> rotateVec, Mesh &mesh);
	bool addBox(std::vector<double> boxVec, Mesh &mesh);
	bool reverseOrient(Eigen::MatrixXi &T);
	bool repair(Mesh &mesh);
	void checkOrientation(Mesh & mesh);
	bool resetOrientation(Mesh &mesh,bool reset_mask=false);
	bool removeDulplicatePoint(Eigen::MatrixXd& V, Eigen::MatrixXi& T, double eps);
	bool createBox(std::vector<double> create_box, Mesh &mesh);
	bool removeBox(Mesh &mesh);
	bool buildTuttleParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);
	bool buildHarmonicParameter(Eigen::MatrixXd& V_3d, Eigen::MatrixXi& T_3d, Eigen::MatrixXd& V_uv);
	bool shuffleSurfaceid(int num, Mesh& mesh);
	bool topoFillHole(Mesh& mesh);
	bool Normalize(Mesh & mesh);

  // called by the algorithm above.
	void dfs_get_loop2(
  int cur, int pre, 
  std::vector<bool>& vis, 
  std::vector<std::vector<int>>& G, 
  std::vector<int>& path, 
  std::vector<std::vector<int>>& loop_lst);
	void boundary_loop_by_dfs2(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::VectorXi &bnd);
	template <
			typename DerivedF,
			typename Derivedb,
			typename VectorIndex,
			typename DerivedF_filled>
	void topological_hole_fill(
			const Eigen::MatrixBase<DerivedF> &F,
			const Eigen::MatrixBase<Derivedb> &b,
			const std::vector<VectorIndex> &holes,
			Eigen::PlainObjectBase<DerivedF_filled> &F_filled,
			Eigen::MatrixXd &V);


	// split origin mesh to several faces by face id and no-manifold feature
// input:
//     mesh: complicate mesh
// Output:
//     meshlist: several meshes
// Return:
//     0  if success;
//     1  otherwise;
	int splitDifferentFaces(const Mesh& mesh, std::vector<Mesh>& meshlist);

	// remove all hanging points in mesh
	// input:
	//     mesh: mesh which has hanging point
	// Return:
	//     0  if success;
	//     1  otherwise;
	int removeHangingPoint(Mesh& mesh);

	// remove all dulplicate topo in mesh
	// input:
	//     mesh: mesh which has dulplicate topo
	// Return:
	//     0  if success;
	//     1  otherwise;
	int removeDulplicateTopo(Mesh& mesh);

	// remove all degradation topo in mesh
	// input:
	//     mesh: mesh which has degradation topo
	// Return:
	//     0  if success;
	//     1  otherwise;
	int removeDegradationTopo(Mesh& mesh);

	// calculate the longest and shortest edges in mesh
	// input:
	//     mesh: 
	// output:
	//     hmax: the longest length
	//     hmin: the shortest length
	//     average: the average length
	// Return:
	//     0  if success;
	//     1  otherwise;
	int calculateEdgesLength(const Mesh& mesh, double &hmax, double &hmin, double &average);

}

#endif
 