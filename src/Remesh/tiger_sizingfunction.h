/** ----------------------------------------------------------------------------
	* TIGER, a Trustable Intelligent GridER
	*
	* TIGER的非结构网格单元尺寸控制函数接口定义 (版本号：1.0)
	* Definitions of APIs for Sizing Functions of Unstructured Grid Generation (Version 1.0) 
	*
	* 陈建军 中国 浙江大学航空航天学院
	* 版权所有	  2020年3月5日
	* Chen Jianjun  School of Aeronautics & Astronautics,
	* Zhejiang University, P. R. China
	* Copyright reserved, Mar 5th, 2020
	* 
	* 联系方式
	*   电话：+86-571-87951883
	*   传真：+86-571-87953167
	*   邮箱：chenjj@zju.edu.cn
	* For further information, please conctact
	*  Tel: +86-571-87951883
	*  Fax: +86-571-87953167
	*  Mail: chenjj@zju.edu.cn
	*
	* ------------------------------------------------------------------------------*/
#ifndef __tiger_sizingfunction_h__
#define __tiger_sizingfunction_h__

#ifdef __linux__
#define __declspec(v)
#endif


#ifdef TigerAdasize_EXPORTS
#define TigerAdasize_API __declspec(dllexport)  
#else
#define TigerAdasize_API __declspec(dllimport)  
#endif
// #include "tiger_sizingsmoother.h"

enum SIZE_LEVEL{
    SIZE_EXTREMELY_FINE = 0,
    SIZE_EXTRA_FINE = 1,
    SIZE_FINER = 2,
    SIZE_FINE = 3,
    SIZE_NORMAL =4,
    SIZE_COARSE = 5,
    SIZE_COARSER = 6,
    SIZE_EXTRA_COARSE  = 7,
    SIZE_EXTREMELY_COARSE = 8
};


enum SF_OBJ_TYPE
{
	SF_OBJ_GLOB_VALUE = 0,
	SF_OBJ_CURVE_MESH,
	SF_OBJ_SURF_MESH,
	SF_OBJ_VOLM_MESH,
	SF_OBJ_SIMPLEX_SOURCES
};

/** -------------------------------------------------------------------------------
	* Given a geometry and a set of user defined element spacing paramters attached to
	* the geometry entities, try to create a sizing function defined by a surface background
	* mesh in which the sizing values are adapted to the geometry features.
	*
	* The inputs include three parts:
	* The first part is a simplified B-rep of the geometry. It defines the geometry of 
	* curves & surfaces by Furgeson interpolations. See Chapter 3 of the following references 
	* for more details.
	* 
	*	郑耀，陈建军，著. 非结构网格生成：理论、算法和应用. 科学出版社, 2016. ISBN: 978-7-03-047283-0.
    *	Yao Zheng, Jianjun Chen, Unstructured Mesh Generation: Theories, Algorithms and Applications (in Chinese),
	*	the 1st Edition, ISBN 978-7-03-047283-0, Science Press, Beijing, China, March, 2016 (273 Pages).
	*	surfaces 
	*
	* For each face, the group of curves trimming the surface is recorded as well.
	* 
	* The second part of the input is the tesselation of the geometry (if any).
	*
	* The third part of the input includes user parameters for element spacing. In the present settings, 
	* these parameters are separated into global and local ones.
	*
	* The global parameters include a set of 7 parameters, organized in the order as below,
	* hmax  	--	the maximally allowed element spacing value.
	* hmin  	--	the minimally allowed element spacing value.
	* beta  	--	the default gradation value of neighbouring elements' spacings
	* theta		-- 	the user parameter for local spacing settings adapted to curvatures (unit: degree).
	* nsegPt	--	the user parameter for local spacing settings adapted to point-point proximities 
	* nsegCv	--  the user parameter for local spacing settings adapted to curve-curve proximities
	* nsegFc	--  the user parameter for local spacing settings adapted to face-face proximities
	*
	* The local parameters are classified by which type of entities the sample points are classified on.
	* 
	*	For a sample point classified on an ending point of a curve, the local parameters include a set of 
	*	6 parameters, organized in the order as below,
	* hmax  	--	the maximally allowed element spacing value.
	* hmin  	--	the minimally allowed element spacing value.
	* theta		-- 	the user parameter for local spacing settings adapted to curvatures (unit: degree).
	* nsegPt	--	the user parameter for local spacing settings adapted to point-point proximities 
	* nsegCv	--  the user parameter for local spacing settings adapted to curve-curve proximities
	* nsegFc	--  the user parameter for local spacing settings adapted to face-face proximities
	*
	* For a sample point classified on a curve, the local parameters include a set of 5 parameters, 
	* organized in the order as below,
	* hmax  	--	the maximally allowed element spacing value.
	* hmin  	--	the minimally allowed element spacing value.
	* theta		-- 	the user parameter for local spacing settings adapted to curvatures (unit: degree).
	* nsegCv	--  the user parameter for local spacing settings adapted to curve-curve proximities
	* nsegFc	--  the user parameter for local spacing settings adapted to face-face proximities
	*
	* For a sample point classified on a surface, the local parameters include a set of 5 parameters, 
	* organized in the order as below,
	* hmax  	--	the maximally allowed element spacing value.
	* hmin  	--	the minimally allowed element spacing value.
	* beta  	--	the default gradation value of neighbouring elements' spacings
	* theta		-- 	the user parameter for local spacing settings adapted to curvatures (unit: degree).
	* nsegFc	--  the user parameter for local spacing settings adapted to face-face proximities
	*
	* The outputs are an ID that could be later used to query the created sizing function object
	* 
	* Input Paramerters:
	*	cvNums			--- # of curves
	*	cvPtNums		--- # of points interpolating a curve. dim = cvNums
	*	cvPts			---	the coordinates of curve points. dim = 3*(cvPtNums[0] + ....+cvPtNums[cvNums-1])
	*	fcNums			---	# of surfaces
	*	fcPtNums		---	# of points along U & V directions interpolating a face. dim = 2 * fcNums
	*	fcPts			---	the coordinates of curve points. dim = 3*(fcPtNums[0]*fcPtNum[1] + ....+fcPtNums[2*fcNums-1]*fcPtNums[2*fcNums-2]
	*	lpCvNums		---	# of curves trimming a surface. dim = fcNums
	*	lpCvs			---	the indices of curves trimming each surface. dim = lpCvNums[0] + .... +lpCvNums[fcNums-1]
	* 	bndPtNum		---	# of boundary points of the tessellation 
	*								= 0 				if no such a tessellation is given
	*								otherwise		a valid tessellation is given 
	*	bndPts			---	the coordinates of boundary points. dim = 3*bndPtNum
	*	bndFctNum		---	# of boundary facets of the tessellation 				
	*   bndFcts			---	the connections of boundary faces. dim = 3*bndFctNum
	* 	spSettingNum	--- # of space settings for different dimensional entities. dim = 3
	*								spSettingsNum[0] -- # of settings for points
	*								spSettingsNum[1] -- # of settings for curves
	* 							    spSettingsNum[2] -- # of settings for faces
	*	spGlSettings	---user parameters for global spacing settings. dim = 8
	* 	spPtSettings	---user parameters for spacing settings on points. dim = 7 * spSettingNum. For each group of settings, 
	*								the first parameter points to the point id,  the remaining 6 parameters are organized as mentioned before.
	* 	spCvSettings	---user parameters for spacing settings on curves. dim = 6 * spSettingNum. For each group of settings, 
	*								the first parameter points to the curve id,  the remaining 6 parameters are organized as mentioned before.
	* 	spFcSettings	---user parameters for spacing settings on curves. dim = 6 * spSettingNum. For each group of settings, 
	*								the first parameter points to the face id,  the remaining 6 parameters are organized as mentioned before.
	* 	volmGrowthRatio	---the growth ratio used to compute sizing values at volume points. 
	* Output Parameters
	*	sfObjID				the ID of the sizing function object
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
TigerAdasize_API int API_Create_SurfBKG_SF(
	int cvNums, int cvPtNums[], double cvPts[], 
	int fcNums, int fcPtNums[], double fcPts[],
	int lpCvNums[], int lpCvs[],
	int bndPtNum, double bndPts[], 
	int bndFctNum, int bndFcts[],int fctFaces[],
	int spSettingNum[],
	double spGlSettings[], double spPtSettings[],
	double spCvSettings[], double spFcSettings[],
	double volmGrowthRatio, 
	int *sfObjID
	);


TigerAdasize_API int API_Create_SurfBKG_SF(
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
	);

TigerAdasize_API int API_Create_SurfBKG_SF(
	int cvNums, int cvPtNums[], double cvPts[], 
	int fcNums, int fcPtNums[], double fcPts[],
	int lpCvNums[], int lpCvs[],
	int bndPtNum, double bndPts[], 
	int bndFctNum, int bndFcts[],int fctFaces[],
	int bdNum, int bdFcNums[], int bdFcs[],
	SIZE_LEVEL level,
	int *sfObjID,
	bool calculateSurfaceProximity = true,	
	double  tolerantFactor = 0.1
	);



TigerAdasize_API int API_Create_SurfBKG_SF(
	int ptNum, double ptCoords[], double ptValues[],
	int elmNum, int elmConns[],
	double volmGrowthRatio,
	int *sfObjID
);

TigerAdasize_API int API_Create_SurfBKG_SF(
	const char *vtkFile,
	double volmGrowthRatio,
	int *sfObjID
);
/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a tetrahedral mesh, try to create a sizing 
	* function to record it
	*	NOTE: 
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	*	mshPtNums		# of mesh points
	*	mshPts			the coordinates of mesh points. dim = 3*mshPtNums
	* ptValues		sizing values at mesh points
	*	volElmNum		# of volume elements
	*	volElms			the connections of volume elements. dim = 4*volElmNum
	* Output Parameters
	*	sfObjID			the ID of the sizing function object
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
TigerAdasize_API int API_Create_VolmBKG_SF(
	int mshPtNum, double mshPts[], double ptValues[],
	int volElmNum, int volElms[],
	int *sfObjID
	);

/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a triangular mesh, try to create a sizing 
	* function to record it
	*
	*	NOTE: no smoothing is done
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	* bndPtNum			# of boundary points of the tessellation 
	*	bndPts				the coordinates of boundary points. dim = 3*bndPtNum
	*	ptValues			sizing values at mesh points
	*	bndFctNum			# of boundary facets of the tessellation 				
	* bndFcts				the connections of boundary faces. dim = 3*bndFctNum
	* volmGrowthRatio	the growth ratio used to compute sizing values at volume points.
	* Output Parameters
	*	sfObjID				the ID of the sizing function object
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
TigerAdasize_API int API_Create_SurfBKG_SF(
	int bndPtNum, double bndPts[], double ptValues[], 
	int bndFctNum, int bndFcts[],
	double volmGrowthRatio,
	int *sfObjID
	);
	
/** -------------------------------------------------------------------------------
	* Given a sizing function defined by a linear mesh, try to create a sizing 
	* function to record it
	*
	*	NOTE: no smoothing is done
	* 1. no smoothing is done
	* 2. linear function is used in default
	*
	* Input Paramerters:
	* bndPtNum			# of boundary points of the tessellation 
	*	bndPts				the coordinates of boundary points. dim = 3*bndPtNum
	*	ptValues			sizing values at mesh points
	*	bndSegNum			# of linear segments of the tessellation 				
	* bndSegs				the connections of linear segements. dim = 2*bndFctNum
	* volmGrowthRatio	the growth ratio used to compute sizing values at volume points.
	* Output Parameters
	*	sfObjID				the ID of the sizing function object
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
TigerAdasize_API int API_Create_CurveBKG_SF(
	int bndPtNum, double bndPts[], double ptValues[], 
	int bndSegNum, int bndSegs[], double growFactor,
	double volmGrowthRatio,
	int *sfObjID
	);

/** -------------------------------------------------------------------------------
	* Given a sizing function defined by simplex sources, try to create a sizing 
	* function to record it
	*
	*	NOTE: no smoothing is done
	* 1. no smoothing is done
	* 2. exponential function is used in default
	*
	* Input Paramerters:
	* sourceNums		# of different types of simplex sources
	*								sourceNums[0]			# of point sources
	*								sourceNums[1]			# of line sources
	*								sourceNums[2]			# of triange sources
	*	centers				the coordinates of source centers. 
	*								dim = sourceNums[0]	+ 2*sourceNums[1]	+ 3*sourceNums[2]	
	*	innerRadii		the innner radii at center positions 
	*								dim = sourceNums[0]	+ 2*sourceNums[1]	+ 3*sourceNums[2]	
	*	outerRadii		the outer radii at center positions 
	*								dim = sourceNums[0]	+ 2*sourceNums[1]	+ 3*sourceNums[2]	
	* ptValues			sizing values at center positions
	*								dim = sourceNums[0]	+ 2*sourceNums[1]	+ 3*sourceNums[2]	
	* Output Parameters
	*	sfObjID				the ID of the sizing function object
	* Returned value:
	*	0				Completely succeed
	*	1				Partially succeed
	*	otherwise		Fail. Error info. will be defined later
	* ------------------------------------------------------------------------------*/
TigerAdasize_API int API_Create_SimplexSources_SF(
	int sourceNums[], double centers[],	double innerRadii[], 
	double outerRadii[], double ptValues[],
	int *sfObjID
	);

TigerAdasize_API int API_Create_SurfBKG_SF(double globalSize, int *sfObjID);

TigerAdasize_API int API_Query_SF_Objects(int *num, int **sfObjIDs);
TigerAdasize_API int API_Enable_SF_Object(int sfObjID);
TigerAdasize_API int API_Disable_SF_Object(int sfObjID);
TigerAdasize_API int API_Query_SF_Object_Validity(int sfObjID, int *valid);
TigerAdasize_API int API_Query_SF_Object_Type(int sfObjID, int *type);

TigerAdasize_API int API_Delete_SF_Object(int sfObjID);

TigerAdasize_API double API_Sizing_Query(double x, double y, double z);
TigerAdasize_API void API_Save_Sizefield(int obj, const char *filePath);
#endif /* __tiger_sizingfunction_h__ */