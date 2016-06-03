
#ifndef BUNDLEOPTIMIZER_H_
#define BUNDLEOPTIMIZER_H_

#include "ceresReprojectionErrors.h"

#include "3dtools.h"
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

	list<list<LatticeClass>> consolidatedLattices; 	/*!< List of consolidated lattices */
	vector<TriangulatedPoint>* allPoints;
	vector<int> camViewIndices; //mapping from camera pose index to image

	vector<ceres::ResidualBlockId> pointReprojectionResiduals; 		/*!< IDs of all point reprojection residuals added to the model */
	vector<ceres::ResidualBlockId> gridTransformationResiduals; 	/*!< IDs of all grid transformation residuals added to the model */
	vector<ceres::ResidualBlockId> basisVectorResiduals; 			/*!< IDs of all basis vector residuals added to the model */

	size_t number_of_cams;			/*!< Number of cameras. */
	size_t num_lattices;			/*!< Number of lattices. */
	size_t numConsolidatedLattices;	/*!< Number of consolidation groups */

	CeresCamModel* CameraModel;
	CeresLattModel* consolidatedLatticeModel;		/*!< Block of lattice basis vector variables that are optimized, for every lattice */
	CeresLattModel* rigidConsolidatedLatticeModel;	/*!< Block of lattice basis vector variables that are optimized, for every consolidation group */

	double principalPoint[2];
	double focalLength;

	ceres::Problem problem;
	FLAGS FLAGS;

	/*!
	 * Initializes num_lattices.
	 */
	void initNumLattices();

	/*!
	 * Initializes consolidatedLatticeModel.
	 */
	void initConsolidatedLatticeModel();

	/*!
	 * Initialize rigidConsolidatedLatticeModel.
	 */
	void initRigidConsolidatedLatticeModel();

	void addPointReprojectionResiduals();

	/*!
	 * Adds grid transformation residuals for every lattice and basis vector residuals for every pair of lattices in a consolidation group,
	 * according to the specified weights.
	 *
	 * @param[in] gridTransformationWeight	The weight for the grid transformation residuals.
	 * @param[in] basisVectorWeight			The weight for the basis vector residuals.
	 */
	void addGridTransformationAndBasisVectorResiduals(double gridTransformationWeight, double basisVectorWeight);

	/*!
	 * Adds grid transformation residuals for every lattice, in the rigid model (lattices in the same consolidation group share the same
	 * basis vector set), according to the specified weight.
	 *
	 * @param[in] rigidGridTransformationWeight	The weight for the rigid grid transformation residuals.
	 */
	void addRigidGridTransformationResiduals(double rigidGridTransformationWeight);


public:

	/*!
	 * Different cost types for different types of residuals, as well as one additional type for total costs
	 */
	enum CostType{POINT_REPROJECTION, GRID_TRANSFORMATION, BASIS_VECTORS, TOTAL};

	/*!
	 * Constructor
	 *
	 * @param[in] aConsolidatedLattices	A list of lattice consolidation groups
	 * @param[in] inpM					A pointer to the input manager that was already used to find the lattice consolidation groups
	 */
	BundleOptimizer(list<list<LatticeClass> > &aConsolidatedLattices, inputManager &inpM);

	/*
	 * Destructor
	 */
	virtual ~BundleOptimizer();

	/*!
	 * Set up a standard optimizer that takes into account only point reprojection errors.
	 */
	void setupStandardOptimizer();

	/*!
	 * Set up an optimizer that takes into account point reprojection errors, as well as additional pairwise constraints to on-grid points of lattices
	 * ("grid transformation residuals").These are defined in a way that every lattice in a consolidation group has its own set of basis vectors,
	 * but there is another set of constraints penalizing deviations between the basis vector sets of lattices in the same consolidation group ("basis
	 * vector residuals").
	 *
	 * @param[in] gridTransformationWeight	The weight for the additional grid transformation residuals.
	 * @param[in] basisVectorWeight			The weight for the additional basis vector residuals.
	 */
	void setupConsolidatedLatticeOptimizer(double gridTransformationWeight, double basisVectorWeight);

	/*!
	 * Set up an optimizer that takes into account point reprojection errors, as well as additional pairwise constraints to on-grid points of lattices
	 * ("rigid grid transformation residuals"). These are defined in a way that all lattices of a consolidation group share the same set of basis vectors.
	 *
	 * @param[in] rigidGridTransformationWeight	The weight for the additional grid transformation residuals.
	 */
	void setupRigidConsolidatedLatticeOptimizer(double rigidGridTransformationWeight);

	void solve();
	vector< Eigen::Matrix<double,3,4> > getOptimizedCameras();

	/*!
	 * Reads out the optimized lattice parameters (basisvectors and the correspondingly updated corner) and writes their new values to
	 * the specified consolidated lattices list. Use this method only if the BundleOptimizer was set up with setupConsolidatedLatticeOptimizer,
	 * otherwise the behavior is undefined.
	 *
	 * @param[in,out] aConsolidatedLattices	The consolidatedLattices object that was already given when initializing the BundleOptimizer.
	 */
	void readoutLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices);

	/*!
	 * Reads out the optimized lattice parameters (basisvectors and the correspondingly updated corner) and writes their new values to
	 * the specified consolidated lattices list. Use this method only if the BundleOptimizer was set up with setupRigidConsolidatedLatticeOptimizer,
	 * otherwise the behavios is undefined.
	 *
	 * @param[in,out] aConsolidatedLattices	The consolidatedLattices object that was already given when initializing the BundleOptimizer.
	 */
	void readoutRigidLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices);

	/*!
	 * Calculate the summed costs of a specified cost type.
	 *
	 * @param[in] type 	The type of costs that should be summed.
	 * @return	The summed costs of the specified cost type.
	 */
	double calculateCost(CostType type);

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
