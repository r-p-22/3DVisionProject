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

	//params[0,1,2]:   Lowerleft position 3D -- DEPRECATED/NOT USED
	//params[3,4,5]: basisVector1 3D (width,a1)
	//params[6,7,8]:basisVector2 3D (height,a2)

	// observed_p is the pivot's 2D coordinates

	//bring 3D point to pivot's position
	T point1[3];
	point1[0] = point[0] - T(a1)*params[3] - T(a2)*params[6];
	point1[1] = point[1] - T(a1)*params[4] - T(a2)*params[7];
	point1[2] = point[2] - T(a1)*params[5] - T(a2)*params[8];

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
	   return (new ceres::AutoDiffCostFunction<LatticeGridToPivotReprojectionError, 2, 6, 9, 3>( //size of residual, size of cam, sizeof params, size of point X
                 new LatticeGridToPivotReprojectionError(observed_x, observed_y, focal, principalPoint_x,principalPoint_y, a1, a2)));
   }

  int a1;
  int a2;

};

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
	//params[0,1,2]:   Lowerleft position 3D
	//params[3,4,5]: basisVector1 3D (width,a1)
	//params[6,7,8]:basisVector2 3D (height,a2)

	// observed_p is the grid point's 2D coordinates

	T pivot[3]; pivot[0] = pivotPoint[0]; pivot[1] = pivotPoint[1]; pivot[2] = pivotPoint[2];

	//bring pivot to 3d point
	pivot[0] += (T(a1)*params[3] + T(a2)*params[6]);
	pivot[1] += (T(a1)*params[4] + T(a2)*params[7]);
	pivot[2] += (T(a1)*params[5] + T(a2)*params[8]);

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
     return (new ceres::AutoDiffCostFunction<LatticePivotToGridReprojectionError, 2, 6, 9,3>( //size of residual, size of cam, size of point X
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
	//params[0,1,2]:   Lowerleft position 3D
	//params[3,4,5]: basisVector1 3D (width,a1,b1)
	//params[6,7,8]:basisVector2 3D (height,a2,b2)

	// observed_p is the grid point's 2D coordinates

	T p[3]; p[0] = pointA[0]; p[1] = pointA[1]; p[2] = pointA[2];

	//bring pointA to lower left
	p[0] -= (T(a1)*params[3] + T(a2)*params[6]);
	p[1] -= (T(a1)*params[4] + T(a2)*params[7]);
	p[2] -= (T(a1)*params[5] + T(a2)*params[8]);


	//bring now bring pointA in the place of the other point (whose observation is available)
	p[0] += (T(b1)*params[3] + T(b2)*params[6]);
	p[1] += (T(b1)*params[4] + T(b2)*params[7]);
	p[2] += (T(b1)*params[5] + T(b2)*params[8]);

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
     return (new ceres::AutoDiffCostFunction<LatticePairwiseReprojectionError, 2, 6, 9,3>( //size of residual, size of cam, size of point X
                 new LatticePairwiseReprojectionError(observed_x, observed_y, focal, principalPoint_x, principalPoint_y, a1, a2,b1,b2)));
   }

  int a1;
  int a2;
  int b1;
  int b2;
};


#endif
