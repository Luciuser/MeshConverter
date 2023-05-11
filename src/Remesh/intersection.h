#pragma once


//namespace Adaptation {
#define REAL double


	//============================================================================//
	//                                                                            //
	// Labels (enumeration declarations) used by TetGen.                          //
	//                                                                            //
	//============================================================================//

	  // Labels that signify the type of a vertex. 
	enum verttype {
		UNUSEDVERTEX, DUPLICATEDVERTEX, RIDGEVERTEX, /*ACUTEVERTEX,*/
		FACETVERTEX, VOLVERTEX, FREESEGVERTEX, FREEFACETVERTEX,
		FREEVOLVERTEX, NREGULARVERTEX, DEADVERTEX
	};

	// Labels that signify the result of triangle-triangle intersection test.
	enum interresult {
		DISJOINT, INTERSECT, SHAREVERT, SHAREEDGE, SHAREFACE,
		TOUCHEDGE, TOUCHFACE, ACROSSVERT, ACROSSEDGE, ACROSSFACE,
		SELF_INTERSECT
	};

	// Labels that signify the result of point location.
	enum locateresult {
		UNKNOWN, OUTSIDE, INTETRAHEDRON, ONFACE, ONEDGE, ONVERTEX,
		ENCVERTEX, ENCSEGMENT, ENCSUBFACE, NEARVERTEX, NONREGULAR,
		INSTAR, BADELEMENT, NULLCAVITY, SHARPCORNER, FENSEDIN,
		NONCOPLANAR, SELF_ENCROACH
	};

#define SETVECTOR3(V, a0, a1, a2) (V)[0] = (a0); (V)[1] = (a1); (V)[2] = (a2)

#define SWAP2(a0, a1, tmp) (tmp) = (a0); (a0) = (a1); (a1) = (tmp)

	// dot() returns the dot product: v1 dot v2.
	REAL dot(REAL* v1, REAL* v2);


	// cross() computes the cross product: n = v1 cross v2.
	void cross(REAL* v1, REAL* v2, REAL* n);

	// distance() computes the Euclidean distance between two points.
	inline REAL distance(REAL* p1, REAL* p2);

	REAL distance2(REAL* p1, REAL* p2);

	REAL norm2(REAL x, REAL y, REAL z);


	void facenormal(REAL *pa, REAL *pb, REAL *pc, REAL *n, int pivot, REAL* lav);


	REAL orient3d(REAL *pa, REAL *pb, REAL *pc, REAL *pd);


// Triangle-edge intersection test (robust)
	int tri_edge_2d(REAL *A, REAL *B, REAL *C, REAL *P, REAL *Q, REAL *R, int level, int *types, int *pos);
	int tri_edge_tail(REAL *A, REAL *B, REAL *C, REAL *P, REAL *Q, REAL *R, REAL sP, REAL sQ, int level, int *types, int *pos);
	int tri_edge_inter(REAL *A, REAL *B, REAL *C, REAL *P, REAL *Q);



// Triangle-triangle intersection test (robust)
	int tri_tri_inter(REAL *A, REAL *B, REAL *C, REAL *O, REAL *P, REAL *Q);
	int tri_edge_inter_tail(REAL *A, REAL *B, REAL *C, REAL *P, REAL *Q, REAL s_p, REAL s_q);












//}