
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
	 *
	 * Contents of camera[6]:
	 * 		0-2:   orientation in (x,y,z) axes,
	 * 		3-5:   pose in (X,Y,Z) of camera center,
	 *
	 * Contents of model[6]:
	 * 		0-2: basisVector1 3D
	 * 		3-5: basisVector2 3D
	 */

	struct CeresCamModel{
		double model[6];
	};
	struct CeresLattModel{
		double model[6];
	};

	struct FLAGS{
		bool robust = true;
	};

	vector<LatticeClass> Lattices;
	list<list<LatticeClass>> consolidatedLattices;
	vector<TriangulatedPoint> allPoints;
	vector<int> camViewIndices; //mapping from camera pose index to image

	size_t number_of_cams;
	size_t num_lattices;
	size_t numConsolidatedLattices;

	CeresCamModel* CameraModel;
	CeresLattModel* LattModel;
	CeresLattModel* consolidatedLatticeModels;
	CeresLattModel* advancedConsolidatedLatticeModels;

	double principalPoint[2];
	double focalLength;

	ceres::Problem problem;
	FLAGS FLAGS;


public:
	BundleOptimizer(vector<LatticeClass> &Lattices, inputManager &inpM);
	virtual ~BundleOptimizer();

	void initLatticeModels();
	void initConsolidatedLatticeModels();
	void initAdvancedConsolidatedLatticeModels();

	void setupNormalOptimizer();
	void setupAllNormalOptimizer();
	void setupGridOptimizer();
	void setupPairwiseLatticeOptimizer();
	void setupPairwiseConsolidatedLatticeOptimizer();
	void setupAdvancedPairwiseConsolidatedLatticeOptimizer();


	void solve();
	vector< Eigen::Matrix<double,3,4> > getOptimizedCameras();

	void setLatticeParameters(vector<LatticeClass> &);
	void setConsolidatedLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices);


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

	    X(1,1) = cos(x);
	    X(1,2) = -sin(x);
	    X(2,1) = sin(x);
	    X(2,2) = cos(x);

	    Y(0,0) = cos(y);
	    Y(0,2) = sin(y);
	    Y(2,0) = -sin(y);
	    Y(2,2) = cos(y);

	    Z(0,0) = cos(z);
	    Z(0,1) = -sin(z);
	    Z(1,0) = sin(z);
	    Z(1,1) = cos(z);

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