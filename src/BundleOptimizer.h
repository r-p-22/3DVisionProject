
#ifndef BUNDLEOPTIMIZER_H_
#define BUNDLEOPTIMIZER_H_

#include "ceresReprojectionErrors.h"

#include "3dtools.h"
//#include "latticeStruct.h"
#include "inputManager.h"
#include "latticeClass.h"

#include <Eigen/Dense>


/* Data Layout:
 *
 * Lattice[i]------- parameters
 *			|
 *			--- TriangulatedPoints(idcs)------3D position __V3d
 *							|
 *							--pointmeasurements: [i]--- viewId
 *										|
 *										--- 2d image position __V2d
 */

class BundleOptimizer {

private:

	/* CeresLattCamModel is a struct to hold the camera pose information.
	 * To optimize with Ceres, the camera parameters must be given in a double array.
	 * Contents of model[15]:
	 * 		0-2:   orientation in (x,y,z) axes,
	 * 		3-5:   pose in (X,Y,Z) of camera center,
	 * 		6-8:   Lowerleft position 3D
	 * 		9-11:  basisVector1 3D
	 * 		12-14: basisVector2 3D
	 * Note that the lattice information is simply ignored in normal bundle adjustment.
	 */
	struct CeresLattCamModel{
		double model[15];
	};

	struct FLAGS{
		bool robust = false;
	};

	vector<LatticeClass> Lattices;
	vector<TriangulatedPoint> allPoints;
	vector<int> camViewIndices; //mapping from camera pose index to image

	size_t number_of_cams;
	size_t num_lattices;

	CeresLattCamModel* CameraModel;

	double principalPoint[2];
	double focalLength;

	ceres::Problem problem;
	FLAGS FLAGS;


public:
	BundleOptimizer(vector<LatticeClass> &Lattices, inputManager &inpM);
	virtual ~BundleOptimizer();

	void setupNormalOptimizer();
	void setupGridOptimizer();
	void solve();
	vector< Eigen::Matrix<double,3,4> > getOptimizedCameras();




	/*
	 * Extract the angles (in all x,y,z axes) and position and place it in an array
	 */
	static void PmatrixToCeres(double* ceres_pose, const Eigen::Matrix<double,3,4> R){

		double xang = atan2(R(2,1), R(2,2));
		double yang = atan2(-R(2,0), sqrt(R(2,1)*R(2,1) + R(2,2)*R(2,2)) );
		double zang = atan2(R(1,0), R(0,0));

		ceres_pose[0] = xang;   ceres_pose[1] = yang;   ceres_pose[2] = zang;
		ceres_pose[3] = R(0,3); ceres_pose[4] = R(1,3); ceres_pose[5] = R(2,3);

		return;
	}

	/*
	 * Compose the P camera matrix given the angles and the position.
	 */
	static Eigen::Matrix<double,3,4> PmatrixFromCeres(const double* cam_pose){

	    double x = cam_pose[0];
	    double y = cam_pose[1];
	    double z = cam_pose[2];

	    Eigen::Matrix3d X =Eigen::Matrix3d::Identity();
	    Eigen::Matrix3d Y =Eigen::Matrix3d::Identity();
	    Eigen::Matrix3d Z =Eigen::Matrix3d::Identity();

	    X(2,2) = cos(x);
	    X(2,3) = -sin(x);
	    X(3,2) = sin(x);
	    X(3,3) = cos(x);

	    Y(1,1) = cos(y);
	    Y(1,3) = sin(y);
	    Y(3,1) = -sin(y);
	    Y(3,3) = cos(y);

	    Z(1,1) = cos(z);
	    Z(1,2) = -sin(z);
	    Z(2,1) = sin(z);
	    Z(2,2) = cos(z);

	    Eigen::Matrix3d R;

	    R = Z*Y;
	    R = R*X;

	    Eigen::Matrix<double,3,4> RT;
	    RT.block<3,3>(0,0) = R;
	    RT(0,3) = cam_pose[3];
	    RT(1,3) = cam_pose[4];
	    RT(2,3) = cam_pose[5];

	    return RT;
	}
};

#endif /* BUNDLEOPTIMIZER_H_ */
