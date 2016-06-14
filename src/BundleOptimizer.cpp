/*
 * BundleOptimizer.cpp
 *
 *  Created on: May 14, 2016
 *      Author: lianos91
 */

#include "BundleOptimizer.h"

#include "ceres/ceres.h"


BundleOptimizer::BundleOptimizer(list<list<LatticeClass> > &aConsolidatedLattices, inputManager &inpM) {

	allPoints = &(inpM.pointModel);
	camViewIndices = inpM.getViewIds();

	number_of_cams = inpM.getCamPoses().size();

	consolidatedLattices = list<list<LatticeClass> >(aConsolidatedLattices);
	numConsolidatedLattices = consolidatedLattices.size();

	initNumLattices();

	focalLength = inpM.getK()(0,0);
	principalPoint[0] = inpM.getK()(0,2);
	principalPoint[1] = inpM.getK()(1,2);

	CameraModel = new CeresCamModel[number_of_cams];

	//change camera representation: From [R|T] to an array of angles (roll, pitch, yaw) and positions.
	for (size_t k=0; k<number_of_cams; k++){
		PmatrixToCeres(CameraModel[k].model, inpM.getCamPoses()[k]);
	}

	initRigidConsolidatedLatticeModel();
	initConsolidatedLatticeModel();
}

BundleOptimizer::~BundleOptimizer() {
	delete [] CameraModel;
	delete [] rigidConsolidatedLatticeModel;
	delete [] consolidatedLatticeModel;
}

void BundleOptimizer::initNumLattices(){

	list<list<LatticeClass>>::iterator consolidatedIt;

	size_t latticeCount = 0;
	for(consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		latticeCount += (*consolidatedIt).size();
	}

	num_lattices = latticeCount;
}

void BundleOptimizer::initRigidConsolidatedLatticeModel(){

	rigidConsolidatedLatticeModel = new CeresLattModel[numConsolidatedLattices];

	list<list<LatticeClass>>::iterator consolidatedIt;

	int consolidatedGroupID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		// iterate over all lattices in the group to find one with consolidationTransformation == 0 (for simplicity)

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			if ((*latticeIt).consolidationTransformation == 0){

				rigidConsolidatedLatticeModel[consolidatedGroupID].model[0] = (*latticeIt).LattStructure.basisVectors[0].x();
				rigidConsolidatedLatticeModel[consolidatedGroupID].model[1] = (*latticeIt).LattStructure.basisVectors[0].y();
				rigidConsolidatedLatticeModel[consolidatedGroupID].model[2] = (*latticeIt).LattStructure.basisVectors[0].z();
				rigidConsolidatedLatticeModel[consolidatedGroupID].model[3] = (*latticeIt).LattStructure.basisVectors[1].x();
				rigidConsolidatedLatticeModel[consolidatedGroupID].model[4] = (*latticeIt).LattStructure.basisVectors[1].y();
				rigidConsolidatedLatticeModel[consolidatedGroupID].model[5] = (*latticeIt).LattStructure.basisVectors[1].z();

				break;
			}
		}
		consolidatedGroupID++;
	}
}

void BundleOptimizer::initConsolidatedLatticeModel(){


	consolidatedLatticeModel = new CeresLattModel[num_lattices];

	list<list<LatticeClass>>::iterator consolidatedIt;

	int latticeID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			consolidatedLatticeModel[latticeID].model[0] = (*latticeIt).LattStructure.basisVectors[0].x();
			consolidatedLatticeModel[latticeID].model[1] = (*latticeIt).LattStructure.basisVectors[0].y();
			consolidatedLatticeModel[latticeID].model[2] = (*latticeIt).LattStructure.basisVectors[0].z();
			consolidatedLatticeModel[latticeID].model[3] = (*latticeIt).LattStructure.basisVectors[1].x();
			consolidatedLatticeModel[latticeID].model[4] = (*latticeIt).LattStructure.basisVectors[1].y();
			consolidatedLatticeModel[latticeID].model[5] = (*latticeIt).LattStructure.basisVectors[1].z();

			latticeID++;
		}
	}
}

void BundleOptimizer::setupStandardOptimizer(){
	addPointReprojectionResiduals();
}

void BundleOptimizer::setupConsolidatedLatticeOptimizer(double gridTransformationWeight, double basisVectorWeight){
	addPointReprojectionResiduals();
	addGridTransformationAndBasisVectorResiduals(gridTransformationWeight, basisVectorWeight);
}

void BundleOptimizer::setupRigidConsolidatedLatticeOptimizer(double rigidGridTransformationWeight){
	addPointReprojectionResiduals();
	addRigidGridTransformationResiduals(rigidGridTransformationWeight);
}

//Bundle adjustment for all points available
void BundleOptimizer::addPointReprojectionResiduals() {

	//create a residual term for each observation of each 3D point of each lattice. ( sum_L{sum_p3D{sum_2dobs{}}} )

	for (int p3d = 0; p3d < (*allPoints).size(); p3d++){
		TriangulatedPoint Tp = (*allPoints)[p3d];

		for (size_t view_id = 0; view_id < Tp.measurements.size(); view_id++ ){ //iteration over image observations of this 3d point

			//cam pose selection
			int imgview = Tp.measurements[view_id].view;
			size_t cam_id = 0;
			for (cam_id = 0; cam_id < number_of_cams; cam_id++){
				if (camViewIndices[cam_id] == imgview)
					break;
			}
			//at this point, CameraModel[cam_id] contains the camera to optimize for this 2d observation.

			//declare cost function for normal reprojection
			ceres::CostFunction* cost_function =
				ReprojectionError::Create(
				   (double)Tp.measurements[view_id].pos[0], //observation.x
				   (double)Tp.measurements[view_id].pos[1], //observation.y
				   focalLength,						  //focal length and principal point are assumed constant
				   principalPoint[0],
				   principalPoint[1]
			);
			//residual error block
			ceres::ResidualBlockId residualID = problem.AddResidualBlock(cost_function,
					   NULL, //if NULL then squared loss
					   CameraModel[cam_id].model,
					   (*allPoints)[p3d].pos.data());

			pointReprojectionResiduals.push_back(residualID);
		}


	}
}

void BundleOptimizer::addRigidGridTransformationResiduals(double rigidGridTransformationWeight){

	ceres::LossFunction* loss_function = new ceres::ScaledLoss(NULL,rigidGridTransformationWeight,ceres::TAKE_OWNERSHIP);;
	ceres::CostFunction* cost_function;

	//create a residual term for each observation of each pair of 3D point on every lattice.

	list<list<LatticeClass> >::iterator consolidatedIt;
	list<LatticeClass>::iterator latticeIt;

	int consolidatedGroupID = 0;

	for(consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			int numLatticeGridPoints = (*latticeIt).latticeGridIndices.size();
			//iterate over all pairs of 3d points: Tp,Tp2

			for (size_t p3d_id=0; p3d_id < numLatticeGridPoints; p3d_id++){ //iteration over 3d points in that lattice

				TriangulatedPoint Tp = (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id].first];
				int a1 = (*latticeIt).latticeGridIndices[p3d_id].second[0];
				int a2 = (*latticeIt).latticeGridIndices[p3d_id].second[1];

				for (size_t p3d_id2=0; p3d_id2 < numLatticeGridPoints; p3d_id2++){ //iteration over 3d points in that lattice

					if (p3d_id == p3d_id2){
						continue;
					}

					TriangulatedPoint Tp2 = (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id2].first];
					int b1 = (*latticeIt).latticeGridIndices[p3d_id2].second[0];
					int b2 = (*latticeIt).latticeGridIndices[p3d_id2].second[1];

					// Incorporate proper basis vector transformation

					int swap;

					switch((*latticeIt).consolidationTransformation){

						case 1: a1 = -a1;
								b1 = -b1;
								break;

						case 2: a2 = -a2;
								b2 = -b2;
								break;

						case 3: a1 = -a1;
								a2 = -a2;
								b1 = -b1;
								b2 = -b2;
								break;

						case 4: swap = a1;
								a1 = a2;
								a2 = swap;
								swap = b1;
								b1 = b2;
								b2 = swap;
								break;

						case 5: swap = a1;
								a1 = a2;
								a2 = swap;
								a1 = -a1;
								swap = b1;
								b1 = b2;
								b2 = swap;
								b1 = -b1;
								break;

						case 6: swap = a1;
								a1 = a2;
								a2 = swap;
								a2 = -a2;
								swap = b1;
								b1 = b2;
								b2 = swap;
								b2 = -b2;
								break;

						case 7: swap = a1;
								a1 = a2;
								a2 = swap;
								a1 = -a1;
								a2 = -a2;
								swap = b1;
								b1 = b2;
								b2 = swap;
								b1 = -b1;
								b2 = -b2;
								break;

						default: break;
					}

					for (size_t view_id = 0; view_id < Tp2.measurements.size(); view_id++ ){ //iteration over image observations of this 3d point

						//cam pose selection
						int imgview = Tp2.measurements[view_id].view;
						size_t cam_id = 0;
						for (cam_id = 0; cam_id < number_of_cams; cam_id++){
							if (camViewIndices[cam_id] == imgview)
								break;
						}
						//at this point, CameraModel[cam_id] contains the camera to optimize for this 2d observation.

						//a1,a2 are the coordinates of the 3D point under optimization, Tp
						//b1,b2 are the coordinates of the other point, where Tp will be transformed to
						//declare cost function for normal reprojection
						cost_function =
							LatticePairwiseReprojectionError::Create(
							   (double)Tp2.measurements[view_id].pos[0], //observation.x
							   (double)Tp2.measurements[view_id].pos[1], //observation.y
							   focalLength,						  //focal length and principal point are assumed constant
							   principalPoint[0],
							   principalPoint[1],
							   a1,
							   a2,
							   b1,
							   b2
						);
						//residual error block
						;

						ceres::ResidualBlockId residualID = problem.AddResidualBlock(cost_function,
								   loss_function, //if NULL then squared loss
								   CameraModel[cam_id].model,
								   rigidConsolidatedLatticeModel[consolidatedGroupID].model,
								   (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id].first].pos.data());

						gridTransformationResiduals.push_back(residualID);
					}
				}
			}
		}
		consolidatedGroupID++;
	}
}

void BundleOptimizer::addGridTransformationAndBasisVectorResiduals(double gridTransformationWeight, double basisVectorWeight){

	ceres::LossFunction* loss_function  = new ceres::ScaledLoss(NULL,gridTransformationWeight,ceres::TAKE_OWNERSHIP);
	ceres::LossFunction* loss_function_bvecs =  new ceres::ScaledLoss(NULL,basisVectorWeight,ceres::TAKE_OWNERSHIP);

	ceres::CostFunction* cost_function;

	list<list<LatticeClass> >::iterator consolidatedIt;
	list<LatticeClass>::iterator latticeIt, latticeIt1, latticeIt2;

	int consolidatedGroupID = 0;
	int latticeID = 0;

	for(consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		// Add penalty functions for differences of corresponding basis vectors between every lattice pair in a consolidated group

		// copy the ID of the first lattice in the consolidated group
		int lattice1ID = latticeID;

		// Iterate over all lattices in the consolidated group
		for(latticeIt1 = (*consolidatedIt).begin(); latticeIt1 != (*consolidatedIt).end(); ++latticeIt1){

			LatticeClass lattice1 = (*latticeIt1);

			// copy and increase the ID of the current lattice by one
			int lattice2ID = lattice1ID+1;

			latticeIt2 = list<LatticeClass>::iterator(latticeIt1);
			latticeIt2++;

			// Iterate over all lattices with higher index in the consolidated group
			for(; latticeIt2 != (*consolidatedIt).end(); ++latticeIt2){

				LatticeClass lattice2 = (*latticeIt2);

				int cTransformation1 = (*latticeIt1).consolidationTransformation;
				int cTransformation2 = (*latticeIt2).consolidationTransformation;

				//add cost function (lattice1, lattice2) for basis vector 0
				cost_function = VectorDifferenceError::Create(cTransformation1, cTransformation2, true);

				ceres::ResidualBlockId residualID = problem.AddResidualBlock(cost_function,
										loss_function_bvecs, //if NULL then squared loss
										consolidatedLatticeModel[lattice1ID].model,
										consolidatedLatticeModel[lattice2ID].model);

				basisVectorResiduals.push_back(residualID);

				//add cost function (lattice1, lattice2) for basis vector 1
				cost_function = VectorDifferenceError::Create(cTransformation1, cTransformation2, false);

				residualID = problem.AddResidualBlock(cost_function,
										loss_function_bvecs, //if NULL then squared loss
										consolidatedLatticeModel[lattice1ID].model,
										consolidatedLatticeModel[lattice2ID].model);

				basisVectorResiduals.push_back(residualID);

				lattice2ID++;
			}
			lattice1ID++;
		}

		// Create a residual term for each observation of each pair of 3D points on the lattice.

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			int numLatticeGridPoints = (*latticeIt).latticeGridIndices.size();
			//iterate over all pairs of 3d points: Tp,Tp2

			for (size_t p3d_id=0; p3d_id < numLatticeGridPoints; p3d_id++){ //iteration over 3d points in that lattice

				TriangulatedPoint Tp = (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id].first];
				int a1 = (*latticeIt).latticeGridIndices[p3d_id].second[0];
				int a2 = (*latticeIt).latticeGridIndices[p3d_id].second[1];


				for (size_t p3d_id2=0; p3d_id2 < numLatticeGridPoints; p3d_id2++){ //iteration over 3d points in that lattice

					if (p3d_id == p3d_id2){
						continue;
					}

					TriangulatedPoint Tp2 = (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id2].first];
					int b1 = (*latticeIt).latticeGridIndices[p3d_id2].second[0];
					int b2 = (*latticeIt).latticeGridIndices[p3d_id2].second[1];

					for (size_t view_id = 0; view_id < Tp2.measurements.size(); view_id++ ){ //iteration over image observations of this 3d point

						//cam pose selection
						int imgview = Tp2.measurements[view_id].view;
						size_t cam_id = 0;
						for (cam_id = 0; cam_id < number_of_cams; cam_id++){
							if (camViewIndices[cam_id] == imgview)
								break;
						}
						//at this point, CameraModel[cam_id] contains the camera to optimize for this 2d observation.

						//a1,a2 are the coordinates of the 3D point under optimization, Tp
						//b1,b2 are the coordinates of the other point, where Tp will be transformed to
						//declare cost function for normal reprojection
						cost_function =
							LatticePairwiseReprojectionError::Create(
							   (double)Tp2.measurements[view_id].pos[0], //observation.x
							   (double)Tp2.measurements[view_id].pos[1], //observation.y
							   focalLength,						  //focal length and principal point are assumed constant
							   principalPoint[0],
							   principalPoint[1],
							   a1,
							   a2,
							   b1,
							   b2
						);
						//residual error block
						;
						ceres::ResidualBlockId residualID = problem.AddResidualBlock(cost_function,
								   loss_function,
								   CameraModel[cam_id].model,
								   consolidatedLatticeModel[latticeID].model,
								   (*allPoints)[(*latticeIt).latticeGridIndices[p3d_id].first].pos.data());

						gridTransformationResiduals.push_back(residualID);
					}
				}
			}
			latticeID++;
		}
		consolidatedGroupID++;
	}
}

void BundleOptimizer::solve() {

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 200;
	options.function_tolerance = 1.e-6;
	options.parameter_tolerance = 1.e-6;
	ceres::Solve(options, &problem, &summary);

	std::cout << summary.FullReport() << "\n";

}

double BundleOptimizer::calculateCost(CostType type){

	ceres::Problem::EvaluateOptions options;

	switch(type){

		case POINT_REPROJECTION: 	if (pointReprojectionResiduals.size() == 0){
										return 0;
									}
									else{
										options.residual_blocks = pointReprojectionResiduals;
									}
									break;

		case GRID_TRANSFORMATION: 	if (gridTransformationResiduals.size() == 0){
										return 0;
									}
									else{
										options.residual_blocks = gridTransformationResiduals;
									}
									break;

		case BASIS_VECTORS:			if (basisVectorResiduals.size() == 0){
										return 0;
									}
									else{
										options.residual_blocks = basisVectorResiduals;
									}
									break;

		case TOTAL: 				break;

		default: 					return -1;

	}

	double totalCost = 0;
	problem.Evaluate(options, &totalCost, nullptr, nullptr, nullptr);

	return totalCost;
}

vector< Eigen::Matrix<double,3,4> > BundleOptimizer::getOptimizedCameras(){

	vector< Eigen::Matrix<double,3,4> > camPoses;
	Eigen::Matrix<double,3,4>* newPose;

	for(size_t i=0; i < number_of_cams; i++){
		newPose = new Eigen::Matrix<double,3,4>(PmatrixFromCeres(CameraModel[i].model));
		camPoses.push_back((*newPose));
	}

	return camPoses;
}

void BundleOptimizer::readoutLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices){

	list<list<LatticeClass>>::iterator consolidatedIt;

	int latticeID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = aConsolidatedLattices.begin(); consolidatedIt != aConsolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			(*latticeIt).LattStructure.basisVectors[0].x() = consolidatedLatticeModel[latticeID].model[0];
			(*latticeIt).LattStructure.basisVectors[0].y() = consolidatedLatticeModel[latticeID].model[1];
			(*latticeIt).LattStructure.basisVectors[0].z() = consolidatedLatticeModel[latticeID].model[2];
			(*latticeIt).LattStructure.basisVectors[1].x() = consolidatedLatticeModel[latticeID].model[3];
			(*latticeIt).LattStructure.basisVectors[1].y() = consolidatedLatticeModel[latticeID].model[4];
			(*latticeIt).LattStructure.basisVectors[1].z() = consolidatedLatticeModel[latticeID].model[5];

			// adjust corner
			int a1 = (*latticeIt).latticeGridIndices[0].second[0];
			int a2 = (*latticeIt).latticeGridIndices[0].second[1];

			(*latticeIt).LattStructure.corner = (*allPoints)[(*latticeIt).latticeGridIndices[0].first].pos -
						a1*(*latticeIt).LattStructure.basisVectors[0] - a2*(*latticeIt).LattStructure.basisVectors[1];

			latticeID++;
		}
	}
}

void BundleOptimizer::readoutRigidLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices){

	int consolidationGroupID = 0;

	list<list<LatticeClass>>::iterator consolidatedIt;

	for (consolidatedIt = aConsolidatedLattices.begin(); consolidatedIt != aConsolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			// write back the basis vectors with respect to the proper transformation

			switch((*latticeIt).consolidationTransformation){

			case 0:

				(*latticeIt).LattStructure.basisVectors[0][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[5]; break;

			case 1:

				(*latticeIt).LattStructure.basisVectors[0][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[5]; break;

			case 2:

				(*latticeIt).LattStructure.basisVectors[0][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[5]; break;

			case 3:

				(*latticeIt).LattStructure.basisVectors[0][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[5]; break;

			case 4:

				(*latticeIt).LattStructure.basisVectors[0][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[2]; break;

			case 5:

				(*latticeIt).LattStructure.basisVectors[0][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[2]; break;

			case 6:

				(*latticeIt).LattStructure.basisVectors[0][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = rigidConsolidatedLatticeModel[consolidationGroupID].model[2]; break;

			case 7:

				(*latticeIt).LattStructure.basisVectors[0][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = - rigidConsolidatedLatticeModel[consolidationGroupID].model[2]; break;

			default: break;

			}

		// adjust corner

		int a1 = (*latticeIt).latticeGridIndices[0].second[0];
		int a2 = (*latticeIt).latticeGridIndices[0].second[1];
		(*latticeIt).LattStructure.corner = (*allPoints)[(*latticeIt).latticeGridIndices[0].first].pos -
				a1*(*latticeIt).LattStructure.basisVectors[0] - a2*(*latticeIt).LattStructure.basisVectors[1];

		}

		consolidationGroupID++;
	}
}


