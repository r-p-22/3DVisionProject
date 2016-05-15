/*
 * BundleOptimizer.cpp
 *
 *  Created on: May 14, 2016
 *      Author: lianos91
 */

#include "BundleOptimizer.h"

#include "ceres/ceres.h"

BundleOptimizer::BundleOptimizer(vector<LatticeClass> &Lattices, inputManager &inpM) {

	allPoints = inpM.getPointModel();
	camViewIndices = inpM.getViewIds();

	number_of_cams = inpM.getCamPoses().size();
	num_lattices = Lattices.size();

	focalLength = inpM.getK()(0,0);
	principalPoint[0] = inpM.getK()(0,2);
	principalPoint[1] = inpM.getK()(1,2);

	CameraModel = new CeresLattCamModel[number_of_cams];

	//change camera representation: From [R|T] to an array of angles (roll, pitch, yaw) and positions.
	for (size_t k=0; k<number_of_cams; k++){
		PmatrixToCeres(CameraModel[k].model, inpM.getCamPoses()[k]);
	}
}

BundleOptimizer::~BundleOptimizer() {
	delete [] CameraModel;
}

/* This is the normal reprojection error
 *
 * for the 3D points, the optimization is performed in-place: The allPoints will contain the new positions
 * for the cameras and Lattices, the corresponding entry of CameraModel is changed.
*/
void BundleOptimizer::setupNormalOptimizer() {

	ceres::LossFunction* loss_function = FLAGS.robust ? new ceres::HuberLoss(1.0) : NULL;

	//create a residual term for each observation of each 3D point of each lattice. ( sum_L{sum_3Dp{sum_2d_obs{}}} )
	for (size_t l = 0; l < num_lattices; ++l) {
		for (size_t i=0; i < Lattices[l].groupPointsIdx.size(); i++){
			//we take all the points in the group
			TriangulatedPoint Tp = allPoints[Lattices[l].groupPointsIdx[i]];
			for (size_t j = 0; j < Tp.measurements.size(); j++ ){
				//cam pose selection
				int imgview = Tp.measurements[j].view;
				size_t k = 0;
				for (k = 0; k < number_of_cams; k++){
					if (camViewIndices[k] == imgview)
						break;
				}
				//at this point, CameraModel[k] contains the camera to optimize for this 2d observation.

				//declare cost function
				ceres::CostFunction* cost_function =
				ReprojectionError::Create(
				   (double)Tp.measurements[j].pos[0], //observation.x
				   (double)Tp.measurements[j].pos[1], //observation.y
				   focalLength,						  //focal length and principal point are assumed constant
	  			   principalPoint[0],
	  			   principalPoint[1]
				);
				//residual error block
				problem.AddResidualBlock(cost_function,
						   loss_function, //if NULL then squared loss
						   CameraModel[k].model,
						   allPoints[Lattices[l].groupPointsIdx[i]].pos.data());
						   //Lattices[l].pointsInGroup[i].data());
			}
		}
	}

}

void BundleOptimizer::solve() {

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = true;

	ceres::Solve(options, &problem, &summary);

	std::cout << summary.FullReport() << "\n";
}

vector< Eigen::Matrix<double,3,4> > BundleOptimizer::getOptimizedCameras(){

	vector< Eigen::Matrix<double,3,4> > camPoses;
	Eigen::Matrix<double,3,4> newPose;

	for(size_t i=0; i < number_of_cams; i++){
		newPose = PmatrixFromCeres(CameraModel[i].model);
		camPoses.push_back(newPose);
	}

	return camPoses;
}

