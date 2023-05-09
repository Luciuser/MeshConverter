/*
 * @Author: Kejie Fu
 * @Date: 2021-08-19 14:16:09
 * @LastEditTime: 2023-05-05 17:43:42
 * @LastEditors: Kejie Fu
 * @Description: 
 * @FilePath: /EMeshGenSizingFunctionNew/src/tiger_sizingfunction.cpp
 */


#include "tiger_sizingfunction.h"
#include "adaMesh.h"
#include "mmgser.h"
#include "backgroundMesh.h"
#include "geosurface.h"
#include <iostream>
#include <unordered_map>
#include <map>
#include <vector>
#include <mutex>
#include "meshRefiner.h"
#include <string>
#include <fstream>
#include "autogridMesh.h"
namespace SIZING_FUNCTION{
	std::map<int, BackgroundMesh> mapBGMesh;
	static int sizingFunctionObj = 0;
	std::mutex mutex;
	static double PI = 3.14159265358979323846264338327950288419716939937510582;
	static double maximumSizingFactors[9]={0.02, 0.035, 0.055, 0.08, 0.1, 0.15, 0.2, 0.3, 0.5};
	static double minimumSizingFactors[9]={0.002, 0.0015, 0.004, 0.01, 0.018, 0.028, 0.04, 0.05, 0.07};
	static double expandRatios[9]={1.3, 1.35, 1.4, 1.45, 1.5, 1.6, 1.7, 1.85,2};
	static double curvatureFactors[9]={0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9,1};
	static double proximityFactors[9]={1, 0.85, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1};
}

 int API_Create_SurfBKG_SF(double globalSize, int *sfObjID){
	(*sfObjID) = ++SIZING_FUNCTION::sizingFunctionObj;

	SIZING_FUNCTION::mapBGMesh[(*sfObjID)].init(globalSize);	
	return 1;
}


int API_Create_SurfBKG_SF(
	int cvNums, int cvPtNums[], double cvPts[], 
	int fcNums, int fcPtNums[], double fcPts[],
	int lpCvNums[], int lpCvs[],
	int bndPtNum, double bndPts[], 
	int bndFctNum, int bndFcts[],int fctFaces[],
	int bdNum, int bdFcNums[], int bdFcs[],
	int spSettingNum[],
	double spGlSettings[], double spPtSettings[],
	double spCvSettings[], double spFcSettings[],
	double volmGrowthRatio, 
	int *sfObjID
){

	(*sfObjID) = ++SIZING_FUNCTION::sizingFunctionObj;	
	double global_spacing = spGlSettings[0];
	double minimum_spacing = spGlSettings[1];
	double expandRatio = spGlSettings[2];//全局尺寸增长率
	double curvatureFactor = spGlSettings[3];//全局曲率角
	double nelm_proxi = spGlSettings[4];//全局点邻近因子

	std::cout << "new size function API" << std::endl;

	if(spSettingNum[0]==8 && spGlSettings[7]<0){
		std::cout << "Try to define sizing function based on discrete geometry" << std::endl;
		AutoGridMesh amesh0;
		amesh0.init(1, bndPtNum, bndPts, bndFctNum, bndFcts, fctFaces,
            global_spacing, minimum_spacing, expandRatio, curvatureFactor, nelm_proxi);



		amesh0.calculateCurvatureAdaptedSizingValues();
		amesh0.calculateProximityAdaptedSizingValues();
		std::cout << "Finish calculate sizing values" << std::endl;

		amesh0.optimizeSizingValues(1000);
		std::cout << "Finish optimization" << std::endl;
		
		SIZING_FUNCTION::mapBGMesh[(*sfObjID)].init(amesh0.pnts, amesh0.tris, amesh0.sizingValues);


		std::cout << "Finish building sizing function!" << std::endl;

	}
	else{

		double tmpx = bndPts[0], tmpy = bndPts[1], tmpz = bndPts[2];
		double xmax = tmpx, ymax=tmpy, zmax = tmpz;
		double xmin = tmpx, ymin = tmpy, zmin = tmpz;
		for(int i=1; i<bndPtNum; i++){
			int base =3*i;
			double xx = bndPts[base];
			double yy = bndPts[base+1];
			double zz = bndPts[base+2];
			xmax = fmax(xmax, xx);
			xmin = fmin(xmin, xx);
			ymax = fmax(ymax, yy);
			ymin = fmin(ymin, yy);
			zmax = fmax(zmax, zz);
			zmin = fmin(zmin, zz);
		}
		double L = fmax(xmax-xmin, fmax((ymax-ymin), zmax-zmin));

		std::cout << "Start init base" << std::endl;
		AdaMesh amesh0, amesh1, amesh2;
		MMGSer remesher;
		amesh0.tolerantFactor = 0.1;
		amesh0.initBase(1, bndPtNum, bndPts,
				bndFctNum, bndFcts, fctFaces,				
				global_spacing,
				minimum_spacing,
				expandRatio,
				curvatureFactor,
				nelm_proxi,
				cvNums, cvPtNums, cvPts,
				fcNums, fcPtNums, fcPts,
				lpCvNums, lpCvs,
				bdNum, bdFcNums, bdFcs);
		std::cout << "Finish init base" << std::endl;
		std::vector<int> edges;
		amesh0.getBaseFeatureEdges(1, edges);
		//amesh0.exportBaseFeatureEdges("/home/kjfu/research/EMeshGen/test/test_sizefield_edges.vtk");
		std::vector<double> tmpSols(bndPtNum);
		Eigen::MatrixXd sv(bndPtNum, 1);
		double minS = std::numeric_limits<double>::max();
		for(int i=0; i<bndPtNum; i++){
			tmpSols[i] = fmin(amesh0.base_curvatureSizingValues[i], amesh0.base_curveProximitySizingValues[i]);
			sv(i, 0) = tmpSols[i];
			minS = fmin(tmpSols[i], minS);
		}

		// amesh0.calculateBaseTriangleQualities();
		// amesh0.exportBase2VTK("/home/kjfu/research/EMeshGen/test/test_sizefield_base.vtk");
		std::vector<double> startPoints;
		std::vector<int> startElements;
		std::vector<int> startRefs;
		amesh0.exportBaseMesh(1, startPoints, startElements, startRefs);
		remesher.initMesh(startPoints.size()/3, startPoints.data(), startElements.size()/3, startElements.data(), startRefs.data(), edges.size()/2, edges.data());

		remesher.setSize(0.001*L, 0.1*L);
		remesher.setHausdorffDistance(0.001*L);	
		remesher.setGradation(3);
		remesher.setOptimization();



		remesher.run();
		
		std::cout << "Finish 1st refine mesh " << std::endl;
		//remesher.exportVTK("/home/kjfu/research/EMeshGen/test/test_sizefield_before1.vtk");

		Eigen::MatrixXd mid_pnts;
		Eigen::MatrixXi mid_tris;
		Eigen::MatrixXi mid_triRefs;
		Eigen::MatrixXd mid_sizingValues;
		remesher.exportMesh(0, mid_pnts, mid_tris, mid_triRefs);



		amesh0.calculateSizingValues(mid_pnts, mid_tris, mid_triRefs, mid_sizingValues);

		remesher.freeSol();
		remesher.initSol(mid_sizingValues.rows(), mid_sizingValues.data());
		remesher.setOptimization(false);
		remesher.setGradation(expandRatio);
		remesher.setHausdorffDistance(0.001*L);	
		
		remesher.setSize(0.001*L, 0.1*L);

		
		remesher.run();


		std::cout << "Finish 2nd refine mesh: remove banding" << std::endl;
		std::vector<double> newPnts;
		std::vector<int> newTris;
		std::vector<int> newRefs;
		remesher.exportMesh(0, newPnts, newTris, newRefs);

		amesh0.initMesh(0, newPnts.size()/3, newPnts.data(), newTris.size()/3, newTris.data(), newRefs.data());

		amesh0.calculateSizingValues(true);
		std::cout << "Finish calculate sizing values" << std::endl;
		amesh0.optimizeSizingValues(1000);
		// amesh0.exportVTK("/home/kjfu/research/EMeshGen/paper_test_revise/background/test_sizefield.vtk");
		double tolerantSize = minimum_spacing;
		std::cout << "Finish optimization" << std::endl;
		for(int i=0; i<amesh0.sizingValues.size(); i++){
			amesh0.sizingValues(i, 0) = fmax(amesh0.sizingValues(i, 0), tolerantSize);
		}

		SIZING_FUNCTION::mapBGMesh[(*sfObjID)].init(amesh0.pnts, amesh0.tris, amesh0.sizingValues);

		std::cout << "Finish building sizing function!" << std::endl;
	}
	return 1;

}

int API_Create_SurfBKG_SF(
	const char *vtkFile,
	double volmGrowthRatio,
	int *sfObjID
){
	(*sfObjID) = ++SIZING_FUNCTION::sizingFunctionObj;

	SIZING_FUNCTION::mapBGMesh[(*sfObjID)].init(vtkFile);	

	return 1;	
}

double API_Sizing_Query( double x, double y, double z) {
	

	double ret=std::numeric_limits<double>::max();

	for(auto &kv: SIZING_FUNCTION::mapBGMesh){
		if (kv.second.isConst){
			ret=kv.second.m_constSize;
			break;
		}
		double s= kv.second.getSize(x,y,z);
		ret = std::min(ret, s);	
	}
	return ret;
}

void API_Save_Sizefield(int obj, const char *filePath){
	std::string ffilePath(filePath);
	SIZING_FUNCTION::mapBGMesh[obj].exportSizeFieldToVTK(ffilePath);
}


TigerAdasize_API int API_Create_SurfBKG_SF(
	int cvNums, int cvPtNums[], double cvPts[], 
	int fcNums, int fcPtNums[], double fcPts[],
	int lpCvNums[], int lpCvs[],
	int bndPtNum, double bndPts[], 
	int bndFctNum, int bndFcts[],int fctFaces[],
	int bdNum, int bdFcNums[], int bdFcs[],
	SIZE_LEVEL level,
	int *sfObjID,
	bool calculateSurfaceProximity,
	double  tolerantFactor
){	
	(*sfObjID) = ++SIZING_FUNCTION::sizingFunctionObj;
	double tmpx = bndPts[0], tmpy = bndPts[1], tmpz = bndPts[2];
	double xmax = tmpx, ymax=tmpy, zmax = tmpz;
	double xmin = tmpx, ymin = tmpy, zmin = tmpz;
	

	for(int i=1; i<bndPtNum; i++){

		int base =3*i;
		double xx = bndPts[base];
		double yy = bndPts[base+1];
		double zz = bndPts[base+2];
		xmax = fmax(xmax, xx);
		xmin = fmin(xmin, xx);
		ymax = fmax(ymax, yy);
		ymin = fmin(ymin, yy);
		zmax = fmax(zmax, zz);
		zmin = fmin(zmin, zz);
	}
	double L = fmax(xmax-xmin, fmax((ymax-ymin), zmax-zmin));


	double global_spacing = SIZING_FUNCTION::maximumSizingFactors[level]*L;
	double minimum_spacing = SIZING_FUNCTION::minimumSizingFactors[level]*L;
	double expandRatio = SIZING_FUNCTION::expandRatios[level];//全局尺寸增长率
	double curvatureFactor = SIZING_FUNCTION::curvatureFactors[level];//全局曲率角
	double nelm_proxi = SIZING_FUNCTION::proximityFactors[level];//全局点邻近因子
	double tolerantSize = minimum_spacing * tolerantFactor;

	std::cout << "Start init base" << std::endl;
	AdaMesh amesh0, amesh1, amesh2;
	MMGSer remesher;
	amesh0.tolerantFactor = tolerantFactor;
	amesh0.initBase(1, bndPtNum, bndPts,
			bndFctNum, bndFcts, fctFaces,				
			global_spacing,
			minimum_spacing,
			expandRatio,
			curvatureFactor,
			nelm_proxi,
			cvNums, cvPtNums, cvPts,
			fcNums, fcPtNums, fcPts,
			lpCvNums, lpCvs,
			bdNum, bdFcNums, bdFcs);
	std::cout << "Finish init base" << std::endl;
	std::vector<int> edges;
	amesh0.getBaseFeatureEdges(1, edges);
	//amesh0.exportBaseFeatureEdges("/home/kjfu/research/EMeshGen/test/test_sizefield_edges.vtk");
	std::vector<double> tmpSols(bndPtNum);
	Eigen::MatrixXd sv(bndPtNum, 1);

	double minS = std::numeric_limits<double>::max();
	for(int i=0; i<bndPtNum; i++){
		tmpSols[i] = fmin(amesh0.base_curvatureSizingValues[i], amesh0.base_curveProximitySizingValues[i]);
		tmpSols[i] = fmax(tolerantSize, tmpSols[i]);
		sv(i, 0) = tmpSols[i];
		minS = fmin(tmpSols[i], minS);
	}

	// amesh0.calculateBaseTriangleQualities();
	// amesh0.exportBase2VTK("/home/kjfu/research/EMeshGen/test/test_sizefield_base.vtk");
	std::vector<double> startPoints;
	std::vector<int> startElements;
	std::vector<int> startRefs;
	amesh0.exportBaseMesh(1, startPoints, startElements, startRefs);
	remesher.initMesh(startPoints.size()/3, startPoints.data(), startElements.size()/3, startElements.data(), startRefs.data(), edges.size()/2, edges.data());

	remesher.setSize(0.001*L, 0.1*L);
	remesher.setHausdorffDistance(0.001*L);	
	remesher.setGradation(3);
	remesher.setOptimization();



	remesher.run();
	
	std::cout << "Finish 1st refine mesh " << std::endl;
	//remesher.exportVTK("/home/kjfu/research/EMeshGen/test/test_sizefield_before1.vtk");

	Eigen::MatrixXd mid_pnts;
	Eigen::MatrixXi mid_tris;
	Eigen::MatrixXi mid_triRefs;
	Eigen::MatrixXd mid_sizingValues;
	remesher.exportMesh(0, mid_pnts, mid_tris, mid_triRefs);



	amesh0.calculateSizingValues(mid_pnts, mid_tris, mid_triRefs, mid_sizingValues);


	remesher.freeSol();
	remesher.initSol(mid_sizingValues.rows(), mid_sizingValues.data());
	remesher.setOptimization(false);
	remesher.setGradation(expandRatio);
	remesher.setHausdorffDistance(0.001*L);	
	
	remesher.setSize(0.001*L, 0.1*L);

	
	remesher.run();


	std::cout << "Finish 2nd refine mesh: remove banding" << std::endl;
	std::vector<double> newPnts;
	std::vector<int> newTris;
	std::vector<int> newRefs;
	remesher.exportMesh(0, newPnts, newTris, newRefs);

	amesh0.initMesh(0, newPnts.size()/3, newPnts.data(), newTris.size()/3, newTris.data(), newRefs.data());

	amesh0.calculateSizingValues(calculateSurfaceProximity);

	std::cout << "Finish calculate sizing values" << std::endl;
	amesh0.optimizeSizingValues(1000);
	// amesh0.exportVTK("/home/kjfu/research/EMeshGen/paper_test_revise/background/test_sizefield.vtk");
	std::cout << "Finish optimization" << std::endl;
	for(int i=0; i<amesh0.sizingValues.size(); i++){
		amesh0.sizingValues(i, 0) = fmax(amesh0.sizingValues(i, 0), tolerantSize);
	}
	SIZING_FUNCTION::mapBGMesh[(*sfObjID)].init(amesh0.pnts, amesh0.tris, amesh0.sizingValues);

	std::cout << "Finish building sizing function!" << std::endl;
	
	return 1;



}