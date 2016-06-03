#ifndef CERESREPROJECTIONERRORS
#define CERESREPROJECTIONERRORS

#include "ceres/ceres.h"
#include "ceres/rotation.h"

struct ReprojectionError {
  ReprojectionError(double observed_x, double observed_y, double focal, double principalPoint_x, double principalPoint_y)
      : observed_x(observed_x), observed_y(observed_y), focal(focal), princPoint_x(principalPoint_x), princPoint_y(principalPoint_y) {}

protected:

  double observed_x;
  double observed_y;
  double focal;
  double princPoint_x;
  double princPoint_y;

  template <typename T>
    void projectInto(const T* const point, const T* const params, T* predicted) const {

 	  T p[3];
 	  T rotparams[3]; rotparams[0] = params[0]; rotparams[1] = params[1]; rotparams[2] = params[2];
 	  ceres::AngleAxisRotatePoint(rotparams, point, p); //takes only the first three elements
 	  p[0] += params[3]; p[1] += params[4]; p[2] += params[5];

 	//p[0] *= camera[6];
 	//p[1] *= camera[6];
 	  p[0] *= T(focal);
 	  p[1] *= T(focal);

 	//No: Compute the center of distortion. The sign change comes from
 	// the camera model that Noah Snavely's Bundler assumes, whereby
 	// the camera coordinate system has a negative z axis.
 	//T xp =  p[0] / p[2];
 	//T yp =  p[1] / p[2];

 	//No: Apply second and fourth order radial distortion.
 	//const T& l1 = camera[7];
 	//const T& l2 = camera[8];
 	//T r2 = xp*xp + yp*yp;
 	//T distortion = T(1.0) + r2  * (l1 + l2  * r2);

 	//add principal point
 	//xp += camera[7];
 	//yp += camera[8];
 	//p[0] += principalPoint[0];
 	//p[1] += principalPoint[1];

 	 predicted[0] = p[0]/p[2] + T(princPoint_x);
     predicted[1] = p[1]/p[2] + T(princPoint_y);

   }

public:

  template <typename T>
  bool operator()(const T* const camera,
                  const T* const point,
                  T* residuals) const {
    // camera[0,1,2] are the angle-axis rotation.
    // camera[3,4,5] are the translation.

	T predicted[2];
	projectInto(point,camera,predicted);

	residuals[0] = predicted[0] - T(observed_x);
    residuals[1] = predicted[1] - T(observed_y);

	return true;
  }

 // Factory to hide the construction of the CostFunction object from
   // the client code.
   static ceres::CostFunction* Create(const double observed_x,
                                      const double observed_y,
                                      const double focal,
									  const double principalPoint_x,
									  const double principalPoint_y) {
     return (new ceres::AutoDiffCostFunction<ReprojectionError, 2, 6, 3>( //size of residual, size of cam, size of point X
                 new ReprojectionError(observed_x, observed_y, focal, principalPoint_x,principalPoint_y)));
   }

};

// TODO Remove?
//from grid point to pivot reprojection error
//the observed_point is the pivot's projection
struct LatticeGridToPivotReprojectionError : ReprojectionError {
	LatticeGridToPivotReprojectionError(double observed_x, double observed_y, double focal,
			double principalPoint_x, double principalPoint_y, int a1, int a2)
	    : ReprojectionError(observed_x, observed_y, focal, principalPoint_x,principalPoint_y),
		a1(a1), a2(a2){};

  template <typename T>
  bool operator()(const T* const camera,
		  	  	  const T* const params,
                  const T* const point,
                  T* residuals) const {

	//camera[0,1,2]:   orientation in (x,y,z) axes
	//camera[3,4,5]:   position in (X,Y,Z) of camera center

	//params[0,1,2]: basisVector1 3D (width,a1)
	//params[3,4,5]:basisVector2 3D (height,a2)

	// observed_p is the pivot's 2D coordinates

	//bring 3D point to pivot's position
	T point1[3];
	point1[0] = point[0] - T(a1)*params[0] - T(a2)*params[3];
	point1[1] = point[1] - T(a1)*params[1] - T(a2)*params[4];
	point1[2] = point[2] - T(a1)*params[2] - T(a2)*params[5];

	T predicted[2];
	projectInto(point1,camera,predicted);

	residuals[0] = predicted[0] - T(observed_x);
    residuals[1] = predicted[1] - T(observed_y);

    return true;
  }

   static ceres::CostFunction* Create(const double observed_x,
                                      const double observed_y,
                                      const double focal,
									  const double principalPoint_x,
									  const double principalPoint_y,
									  const int a1,
									  const int a2){
	   return (new ceres::AutoDiffCostFunction<LatticeGridToPivotReprojectionError, 2, 6, 6, 3>( //size of residual, size of cam, sizeof params, size of point X
                 new LatticeGridToPivotReprojectionError(observed_x, observed_y, focal, principalPoint_x,principalPoint_y, a1, a2)));
   }

  int a1;
  int a2;

};

//TODO Remove?
struct LatticePivotToGridReprojectionError : ReprojectionError {
	LatticePivotToGridReprojectionError(double observed_x, double observed_y, double focal, double principalPoint_x, double principalPoint_y,int a1, int a2)
    : ReprojectionError(observed_x, observed_y, focal, principalPoint_x, principalPoint_y),
	a1(a1), a2(a2){};

  template <typename T>
  bool operator()(const T* const camera,
		  	  	  const T* const params,
		  	  	  const T* const pivotPoint,
                  T* residuals) const {

	//camera[0,1,2]:   orientation in (x,y,z) axes
	//camera[3,4,5]:   position in (X,Y,Z) of camera center,
	//params[0,1,2]: basisVector1 3D (width,a1)
	//params[3,4,5]: basisVector2 3D (height,a2)

	// observed_p is the grid point's 2D coordinates

	T pivot[3]; pivot[0] = pivotPoint[0]; pivot[1] = pivotPoint[1]; pivot[2] = pivotPoint[2];

	//bring pivot to 3d point
	pivot[0] += (T(a1)*params[0] + T(a2)*params[3]);
	pivot[1] += (T(a1)*params[1] + T(a2)*params[4]);
	pivot[2] += (T(a1)*params[2] + T(a2)*params[5]);

	T predicted[2];
	projectInto(pivot,camera,predicted);

	residuals[0] = predicted[0] - T(observed_x);
	residuals[1] = predicted[1] - T(observed_y);

    return true;

  }

   static ceres::CostFunction* Create(const double observed_x,
                                      const double observed_y,
                                      const double focal,
                                      const double principalPoint_x,
                                      const double principalPoint_y,
                                      const int a1,
                                      const int a2) {
     return (new ceres::AutoDiffCostFunction<LatticePivotToGridReprojectionError, 2, 6, 6, 3>( //size of residual, size of cam, size of params, size of point X
                 new LatticePivotToGridReprojectionError(observed_x, observed_y, focal, principalPoint_x, principalPoint_y, a1, a2)));
   }

  int a1;
  int a2;

};


struct LatticePairwiseReprojectionError : ReprojectionError {
	LatticePairwiseReprojectionError(double observed_x, double observed_y, double focal, double principalPoint_x, double principalPoint_y,int a1, int a2, int b1, int b2)
    : ReprojectionError(observed_x, observed_y, focal, principalPoint_x, principalPoint_y),
	a1(a1), a2(a2), b1(b1), b2(b2){};

  template <typename T>
  bool operator()(const T* const camera,
		  	  	  const T* const params,
		  	  	  const T* const pointA,
                  T* residuals) const {

	//camera[0,1,2]:   orientation in (x,y,z) axes
	//camera[3,4,5]:   position in (X,Y,Z) of camera center,
	//params[0,1,2]: 	basisVector1 3D (width,a1,b1)
	//params[3,4,5]:	basisVector2 3D (height,a2,b2)

	// observed_p is the grid point's 2D coordinates

	T p[3]; p[0] = pointA[0]; p[1] = pointA[1]; p[2] = pointA[2];

	//bring pointA to lower left
	p[0] -= (T(a1)*params[0] + T(a2)*params[3]);
	p[1] -= (T(a1)*params[1] + T(a2)*params[4]);
	p[2] -= (T(a1)*params[2] + T(a2)*params[5]);


	//bring now bring pointA in the place of the other point (whose observation is available)
	p[0] += (T(b1)*params[0] + T(b2)*params[3]);
	p[1] += (T(b1)*params[1] + T(b2)*params[4]);
	p[2] += (T(b1)*params[2] + T(b2)*params[5]);

	T predicted[2];
	projectInto(p,camera,predicted);

	residuals[0] = predicted[0] - T(observed_x);
	residuals[1] = predicted[1] - T(observed_y);

    return true;

  }

   static ceres::CostFunction* Create(const double observed_x,
                                      const double observed_y,
                                      const double focal,
                                      const double principalPoint_x,
                                      const double principalPoint_y,
                                      const int a1,
                                      const int a2,
                                      const int b1,
                                      const int b2) {
     return (new ceres::AutoDiffCostFunction<LatticePairwiseReprojectionError, 2, 6, 6, 3>( //size of residual, size of cam, size of params, size of point X
                 new LatticePairwiseReprojectionError(observed_x, observed_y, focal, principalPoint_x, principalPoint_y, a1, a2,b1,b2)));
   }

  int a1;
  int a2;
  int b1;
  int b2;
};

/*!
 * Helper function to transform basis vectors according to a specified transformation.
 *
 * @param[in] transformation	The transformation to apply to the basis vectors.
 * @param[in] basisVector0		The first basis vector to be transformed.
 * @param[in] basisVector1  	The second basis vector to be transformed.
 * @param[out] basisVector0Transformed	The first basis vector after transformation (Note that a renaming can have taken place, such that this vector
 * 										might also originate from basisVector1).
 * @param[out] basisVector1Transformed	The second basis vector after transformation (Note that a renaming can have taken place, such that this vector
 * 										might also originate from basisVector0).
 */
template <typename T>
static int transformBasisVectors(int transformation, const T* const basisVector0, const T* const basisVector1, T* basisVector0Transformed, T* basisVector1Transformed){

		switch(transformation){

			case 0: basisVector0Transformed[0] = basisVector0[0];
					basisVector0Transformed[1] = basisVector0[1];
					basisVector0Transformed[2] = basisVector0[2];

					basisVector1Transformed[0] = basisVector1[0];
					basisVector1Transformed[1] = basisVector1[1];
					basisVector1Transformed[2] = basisVector1[2];

					return 0;

			case 1: basisVector0Transformed[0] = - basisVector0[0];
					basisVector0Transformed[1] = - basisVector0[1];
					basisVector0Transformed[2] = - basisVector0[2];

					basisVector1Transformed[0] = basisVector1[0];
					basisVector1Transformed[1] = basisVector1[1];
					basisVector1Transformed[2] = basisVector1[2];

					return 0;

			case 2: basisVector0Transformed[0] = basisVector0[0];
					basisVector0Transformed[1] = basisVector0[1];
					basisVector0Transformed[2] = basisVector0[2];

					basisVector1Transformed[0] = - basisVector1[0];
					basisVector1Transformed[1] = - basisVector1[1];
					basisVector1Transformed[2] = - basisVector1[2];

					return 0;

			case 3: basisVector0Transformed[0] = - basisVector0[0];
					basisVector0Transformed[1] = - basisVector0[1];
					basisVector0Transformed[2] = - basisVector0[2];

					basisVector1Transformed[0] = - basisVector1[0];
					basisVector1Transformed[1] = - basisVector1[1];
					basisVector1Transformed[2] = - basisVector1[2];
					return 0;

			case 4: basisVector0Transformed[0] = basisVector1[0];
					basisVector0Transformed[1] = basisVector1[1];
					basisVector0Transformed[2] = basisVector1[2];

					basisVector1Transformed[0] = basisVector0[0];
					basisVector1Transformed[1] = basisVector0[1];
					basisVector1Transformed[2] = basisVector0[2];
					return 0;

			case 5: basisVector0Transformed[0] = - basisVector1[0];
					basisVector0Transformed[1] = - basisVector1[1];
					basisVector0Transformed[2] = - basisVector1[2];

					basisVector1Transformed[0] = basisVector0[0];
					basisVector1Transformed[1] = basisVector0[1];
					basisVector1Transformed[2] = basisVector0[2];
					return 0;

			case 6: basisVector0Transformed[0] = basisVector1[0];
					basisVector0Transformed[1] = basisVector1[1];
					basisVector0Transformed[2] = basisVector1[2];

					basisVector1Transformed[0] = - basisVector0[0];
					basisVector1Transformed[1] = - basisVector0[1];
					basisVector1Transformed[2] = - basisVector0[2];
					return 0;

			case 7: basisVector0Transformed[0] = - basisVector1[0];
					basisVector0Transformed[1] = - basisVector1[1];
					basisVector0Transformed[2] = - basisVector1[2];

					basisVector1Transformed[0] = - basisVector0[0];
					basisVector1Transformed[1] = - basisVector0[1];
					basisVector1Transformed[2] = - basisVector0[2];
					return 0;

			default: return -1;
		}
}

/*!
 * Struct for the error term for basis vector difference.
 */
struct VectorDifferenceError {

	// consolidation transformation parameter of the first basisvector set
	int cTransformation1;

	// consolidation transformation parameter of the second basisvector set
	int cTransformation2;

	// if true, returns the residual for the difference in basisvector 0 (transformed),
	// if false, returns the residual for difference in basisvector 1 (transformed)
	bool selectVector0;

	VectorDifferenceError(int aCTransformation1, int aCTransformation2, bool aSelectVector0):
		cTransformation1(aCTransformation1), cTransformation2(aCTransformation2), selectVector0(aSelectVector0){};

	template <typename T>
	bool operator()(	const T* const basisVectorSet1,
						const T* const basisVectorSet2,
						T* residuals) const {

		//basisVectorSet[0,1,2]: 	basisVector0
		//basisVectorSet[3,4,5]:	basisVector1

		// observed_p is the grid point's 2D coordinates

		T  basisVector10[3], basisVector11[3], basisVector20[3], basisVector21[3], finalBasisVector10[3], finalBasisVector11[3], finalBasisVector20[3], finalBasisVector21[3];

		basisVector10[0] = basisVectorSet1[0];
		basisVector10[1] = basisVectorSet1[1];
		basisVector10[2] = basisVectorSet1[2];

		basisVector11[0] = basisVectorSet1[3];
		basisVector11[1] = basisVectorSet1[4];
		basisVector11[2] = basisVectorSet1[5];

		basisVector20[0] = basisVectorSet2[0];
		basisVector20[1] = basisVectorSet2[1];
		basisVector20[2] = basisVectorSet2[2];

		basisVector21[0] = basisVectorSet2[3];
		basisVector21[1] = basisVectorSet2[4];
		basisVector21[2] = basisVectorSet2[5];

		// Transformation of both basis vector sets to common reference system
		transformBasisVectors(cTransformation1, basisVector10, basisVector11, finalBasisVector10, finalBasisVector11);
		transformBasisVectors(cTransformation2, basisVector20, basisVector21, finalBasisVector20, finalBasisVector21);

		T diff0[3], diff1[3];

		diff0[0] = finalBasisVector10[0] - finalBasisVector20[0];
		diff0[1] = finalBasisVector10[1] - finalBasisVector20[1];
		diff0[2] = finalBasisVector10[2] - finalBasisVector20[2];

		diff1[0] = finalBasisVector11[0] - finalBasisVector21[0];
		diff1[1] = finalBasisVector11[1] - finalBasisVector21[1];
		diff1[2] = finalBasisVector11[2] - finalBasisVector21[2];

		if(selectVector0){
			residuals[0] = diff0[0];
			residuals[1] = diff0[1];
			residuals[2] = diff0[2];
		}
		else{
			residuals[0] = diff1[0];
			residuals[1] = diff1[1];
			residuals[2] = diff1[2];
		}

		return true;

	}

	static ceres::CostFunction* Create(const int aCTransformation1,
									  const int aCTransformation2,
									  bool aSelectVector0) {
		return (new ceres::AutoDiffCostFunction<VectorDifferenceError, 3, 6, 6>( //size of residual, size of basisVectorSet1, size of basisVectorSet2
				 new VectorDifferenceError(aCTransformation1, aCTransformation2, aSelectVector0)));
	}

};


#endif
