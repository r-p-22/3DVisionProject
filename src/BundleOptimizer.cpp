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

	consolidatedLattices = LatticeClass::consolidateLattices(Lattices);
	numConsolidatedLattices = consolidatedLattices.size();

	this->Lattices = vector<LatticeClass>(Lattices);

	focalLength = inpM.getK()(0,0);
	principalPoint[0] = inpM.getK()(0,2);
	principalPoint[1] = inpM.getK()(1,2);

	CameraModel = new CeresCamModel[number_of_cams];

	//change camera representation: From [R|T] to an array of angles (roll, pitch, yaw) and positions.
	for (size_t k=0; k<number_of_cams; k++){
		PmatrixToCeres(CameraModel[k].model, inpM.getCamPoses()[k]);
	}

	initLatticeModels();
	initConsolidatedLatticeModels();
	initAdvancedConsolidatedLatticeModels();
}

BundleOptimizer::~BundleOptimizer() {
	delete [] CameraModel;
	delete [] LattModel;
	delete [] consolidatedLatticeModels;
	delete [] advancedConsolidatedLatticeModels;
}

void BundleOptimizer::initLatticeModels(){

	/* Contents of LattModel.model[6]:
	 * 		0-2: basisVector1 3D
	 * 		3-5: basisVector2 3D
	 */
	LattModel = new CeresLattModel[num_lattices];
	for (size_t k=0; k<num_lattices; k++){

		LattModel[k].model[0] = Lattices[k].LattStructure.basisVectors[0].x();
		LattModel[k].model[1] = Lattices[k].LattStructure.basisVectors[0].y();
		LattModel[k].model[2] = Lattices[k].LattStructure.basisVectors[0].z();
		LattModel[k].model[3] = Lattices[k].LattStructure.basisVectors[1].x();
		LattModel[k].model[4] = Lattices[k].LattStructure.basisVectors[1].y();
		LattModel[k].model[5] = Lattices[k].LattStructure.basisVectors[1].z();
	}
}

void BundleOptimizer::initConsolidatedLatticeModels(){

	consolidatedLatticeModels = new CeresLattModel[numConsolidatedLattices];

	list<list<LatticeClass>>::iterator consolidatedIt;

	int consolidatedGroupID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		// iterate over all lattices in the group to find one with consolidationTransformation == 0 (for simplicity)

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			if ((*latticeIt).consolidationTransformation == 0){

				consolidatedLatticeModels[consolidatedGroupID].model[0] = (*latticeIt).LattStructure.basisVectors[0].x();
				consolidatedLatticeModels[consolidatedGroupID].model[1] = (*latticeIt).LattStructure.basisVectors[0].y();
				consolidatedLatticeModels[consolidatedGroupID].model[2] = (*latticeIt).LattStructure.basisVectors[0].z();
				consolidatedLatticeModels[consolidatedGroupID].model[3] = (*latticeIt).LattStructure.basisVectors[1].x();
				consolidatedLatticeModels[consolidatedGroupID].model[4] = (*latticeIt).LattStructure.basisVectors[1].y();
				consolidatedLatticeModels[consolidatedGroupID].model[5] = (*latticeIt).LattStructure.basisVectors[1].z();

				break;
			}
		}
		consolidatedGroupID++;
	}
}

void BundleOptimizer::initAdvancedConsolidatedLatticeModels(){


	advancedConsolidatedLatticeModels = new CeresLattModel[num_lattices];

	list<list<LatticeClass>>::iterator consolidatedIt;

	int latticeID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			advancedConsolidatedLatticeModels[latticeID].model[0] = (*latticeIt).LattStructure.basisVectors[0].x();
			advancedConsolidatedLatticeModels[latticeID].model[1] = (*latticeIt).LattStructure.basisVectors[0].y();
			advancedConsolidatedLatticeModels[latticeID].model[2] = (*latticeIt).LattStructure.basisVectors[0].z();
			advancedConsolidatedLatticeModels[latticeID].model[3] = (*latticeIt).LattStructure.basisVectors[1].x();
			advancedConsolidatedLatticeModels[latticeID].model[4] = (*latticeIt).LattStructure.basisVectors[1].y();
			advancedConsolidatedLatticeModels[latticeID].model[5] = (*latticeIt).LattStructure.basisVectors[1].z();

			latticeID++;
		}
	}
}



/* This is the normal reprojection error
 *
 * for the 3D points, the optimization is performed in-place: The allPoints will contain the new positions
 * for the cameras and Lattices, the corresponding entry of CameraModel is changed.
*/
void BundleOptimizer::setupNormalOptimizer() {

	ceres::LossFunction* loss_function = FLAGS.robust ? new ceres::HuberLoss(1.0) : NULL;

	//create a residual term for each observation of each 3D point of each lattice. ( sum_L{sum_p3D{sum_2dobs{}}} )
	for (size_t latt_id = 0; latt_id < num_lattices; ++latt_id) { //iteration over lattices

		//for (size_t p3d_id=0; p3d_id < Lattices[latt_id].groupPointsIdx.size(); p3d_id++){ //iteration over 3d points in that lattice
		for (size_t p3d_id=0; p3d_id < Lattices[latt_id].latticeGridIndices.size(); p3d_id++){ //iteration over 3d points in that lattice

			TriangulatedPoint Tp = allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first];

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
				problem.AddResidualBlock(cost_function,
						   loss_function, //if NULL then squared loss
						   CameraModel[cam_id].model,
						   allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first].pos.data());
			}

		}
	}
}

void BundleOptimizer::setupAllNormalOptimizer() {

	//ceres::LossFunction* loss_function= NULL;//new ceres::HuberLoss(1.00);
	ceres::LossFunction* loss_function= new ceres::HuberLoss(15.00);

	//create a residual term for each observation of each 3D point of each lattice. ( sum_L{sum_p3D{sum_2dobs{}}} )

	for (int p3d = 0; p3d < allPoints.size(); p3d++){
		TriangulatedPoint Tp = allPoints[p3d];

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
			problem.AddResidualBlock(cost_function,
					   loss_function, //if NULL then squared loss
					   CameraModel[cam_id].model,
					   allPoints[p3d].pos.data());
		}


	}
}



/*void BundleOptimizer::setupGridOptimizer(){

	ceres::LossFunction* loss_function =  new ceres::ScaledLoss(NULL,10,ceres::TAKE_OWNERSHIP);
	ceres::CostFunction* cost_function;
	//create a residual term for each observation of each 3D point of each lattice. ( sum_L{sum_p3D{sum_2dobs{}}} )
	for (size_t latt_id = 0; latt_id < num_lattices; ++latt_id) { //iteration over lattices

		//must first find pivot index
		vector<pair<int, vector<int> >>::iterator pivot = std::find_if( Lattices[latt_id].latticeGridIndices.begin(),
				Lattices[latt_id].latticeGridIndices.end(), find_pivot());
		int pivotIndex = pivot->first;
		TriangulatedPoint Tp_pivot = allPoints[pivotIndex];

		for (size_t view_id = 0; view_id < Tp_pivot.measurements.size(); view_id++ ){ //iteration over image observations of pivot				cout << "view id = " << view_id << endl;

			//cam pose selection
			int imgview = Tp_pivot.measurements[view_id].view;
			size_t cam_id = 0;
			for (cam_id = 0; cam_id < number_of_cams; cam_id++){
				if (camViewIndices[cam_id] == imgview)
					break;
			}

			for (size_t p3d_id=0; p3d_id < Lattices[latt_id].latticeGridIndices.size(); p3d_id++){ //iteration over 3d points in that lattice
				if (Lattices[latt_id].latticeGridIndices[p3d_id].first == pivotIndex)
					continue;
				int a1 = Lattices[latt_id].latticeGridIndices[p3d_id].second[0];
				int a2 = Lattices[latt_id].latticeGridIndices[p3d_id].second[1];

				//at this point, CameraModel[cam_id] contains the camera to optimize for this 2d observation.

				//declare cost function for lattice reprojection from grid point to pivot
				//ceres::CostFunction*
				cost_function =
					LatticeGridToPivotReprojectionError::Create(
					   (double)Tp_pivot.measurements[view_id].pos[0], //observation.x
					   (double)Tp_pivot.measurements[view_id].pos[1], //observation.y
					   focalLength,						  //focal length and principal point are assumed constant
					   principalPoint[0],
					   principalPoint[1],
					   a1,
					   a2
						);
				//residual error block
				problem.AddResidualBlock(cost_function,
					   loss_function, //if NULL then squared loss
					   CameraModel[cam_id].model,
					   LattModel[latt_id].model,
					   allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first].pos.data());

			}

		}

		//for (size_t p3d_id=0; p3d_id < Lattices[latt_id].groupPointsIdx.size(); p3d_id++){ //iteration over 3d points in that lattice
		for (size_t p3d_id=0; p3d_id < Lattices[latt_id].latticeGridIndices.size(); p3d_id++){ //iteration over 3d points in that lattice

			TriangulatedPoint Tp = allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first];
			int a1 = Lattices[latt_id].latticeGridIndices[p3d_id].second[0];
			int a2 = Lattices[latt_id].latticeGridIndices[p3d_id].second[1];


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
				//ceres::CostFunction*
				cost_function =
					ReprojectionError::Create(
					   (double)Tp.measurements[view_id].pos[0], //observation.x
					   (double)Tp.measurements[view_id].pos[1], //observation.y
					   focalLength,						  //focal length and principal point are assumed constant
					   principalPoint[0],
					   principalPoint[1]
				);
				//residual error block
				problem.AddResidualBlock(cost_function,
						loss_function, //if NULL then squared loss
						   CameraModel[cam_id].model,
						   allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first].pos.data());

				//declare cost function for lattice reprojection from pivot to grid point
				//ceres::CostFunction*
				cost_function =
						LatticePivotToGridReprojectionError::Create(
						   (double)Tp.measurements[view_id].pos[0], //observation.x
						   (double)Tp.measurements[view_id].pos[1], //observation.y
						   focalLength,						  //focal length and principal point are assumed constant
						   principalPoint[0],
						   principalPoint[1],
						   a1,
						   a2
							);
				//residual error block
				problem.AddResidualBlock(cost_function,
						   loss_function, //if NULL then squared loss
						   CameraModel[cam_id].model,
						   LattModel[latt_id].model,
						   allPoints[pivotIndex].pos.data());

			}

		}

	}

}*/



void BundleOptimizer::setupPairwiseLatticeOptimizer(){

	ceres::LossFunction* loss_function =  new ceres::ScaledLoss(new ceres::HuberLoss(100),100,ceres::TAKE_OWNERSHIP);
	ceres::CostFunction* cost_function;
	//create a residual term for each observation of each pair of 3D point on every lattice.
	for (size_t latt_id = 0; latt_id < num_lattices; ++latt_id) { //iteration over lattices

		//iterate over all pairs of 3d points: Tp,Tp2

		for (size_t p3d_id=0; p3d_id < Lattices[latt_id].latticeGridIndices.size(); p3d_id++){ //iteration over 3d points in that lattice

			TriangulatedPoint Tp = allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first];
			int a1 = Lattices[latt_id].latticeGridIndices[p3d_id].second[0];
			int a2 = Lattices[latt_id].latticeGridIndices[p3d_id].second[1];


			for (size_t p3d_id2=0; p3d_id2 < Lattices[latt_id].latticeGridIndices.size(); p3d_id2++){ //iteration over 3d points in that lattice

				if (p3d_id == p3d_id2){
					continue;
				}

				TriangulatedPoint Tp2 = allPoints[Lattices[latt_id].latticeGridIndices[p3d_id2].first];
				int b1 = Lattices[latt_id].latticeGridIndices[p3d_id2].second[0];
				int b2 = Lattices[latt_id].latticeGridIndices[p3d_id2].second[1];

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
					problem.AddResidualBlock(cost_function,
							loss_function, //if NULL then squared loss
							   CameraModel[cam_id].model,
							   LattModel[latt_id].model,
							   allPoints[Lattices[latt_id].latticeGridIndices[p3d_id].first].pos.data());
				}

			}

		}

	}

}

void BundleOptimizer::setupPairwiseConsolidatedLatticeOptimizer(){

	ceres::LossFunction* loss_function = FLAGS.robust ? new ceres::HuberLoss(1.50) : NULL;
	ceres::CostFunction* cost_function;

	//create a residual term for each observation of each pair of 3D point on every lattice.

	list<list<LatticeClass> >::iterator consolidatedIt;
	list<LatticeClass>::iterator latticeIt;

	int consolidatedGroupID = 0;

	for(consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			cout << "*** new lattice" << endl;

			int numLatticeGridPoints = (*latticeIt).latticeGridIndices.size();
			//iterate over all pairs of 3d points: Tp,Tp2

			for (size_t p3d_id=0; p3d_id < numLatticeGridPoints; p3d_id++){ //iteration over 3d points in that lattice

				TriangulatedPoint Tp = allPoints[(*latticeIt).latticeGridIndices[p3d_id].first];
				int a1 = (*latticeIt).latticeGridIndices[p3d_id].second[0];
				int a2 = (*latticeIt).latticeGridIndices[p3d_id].second[1];


				for (size_t p3d_id2=0; p3d_id2 < numLatticeGridPoints; p3d_id2++){ //iteration over 3d points in that lattice

					if (p3d_id == p3d_id2){
						continue;
					}

					TriangulatedPoint Tp2 = allPoints[(*latticeIt).latticeGridIndices[p3d_id2].first];
					int b1 = (*latticeIt).latticeGridIndices[p3d_id2].second[0];
					int b2 = (*latticeIt).latticeGridIndices[p3d_id2].second[1];

					// Incorporate proper basis vector transformation

					//TODO Remove these
					int a1old = a1;
					int a2old = a2;
					int b1old = b1;
					int b2old = b2;

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


					/* Debugging outputs

					Vector3d basisVector0 = (*latticeIt).LattStructure.basisVectors[0];
					Vector3d basisVector1 = (*latticeIt).LattStructure.basisVectors[1];
					Vector3d consolidatedBasisVector0 = Vector3d(consolidatedLatticeModels[consolidatedGroupID].model[0], consolidatedLatticeModels[consolidatedGroupID].model[1],consolidatedLatticeModels[consolidatedGroupID].model[2]);
					Vector3d consolidatedBasisVector1 = Vector3d(consolidatedLatticeModels[consolidatedGroupID].model[3], consolidatedLatticeModels[consolidatedGroupID].model[4],consolidatedLatticeModels[consolidatedGroupID].model[5]);

					Vector3d aOld = a1old*basisVector0 + a2old*basisVector1;
					Vector3d aNew = a1*consolidatedBasisVector0 + a2*consolidatedBasisVector1;
					Vector3d bOld = b1old*basisVector0 + b2old*basisVector1;
					Vector3d bNew = b1*consolidatedBasisVector0 + b2*consolidatedBasisVector1;
					Vector3d diffA = aOld - aNew;
					Vector3d diffB = bOld - bNew;

					if (diffA.norm() != 0){
						cout << "***" << endl;
						cout << "diffA: " << endl;
						cout << diffA << endl;
						cout << "diffA norm: " << diffA.norm() << endl;
					}
					if (diffB.norm() != 0){
						cout << "***" << endl;
						cout << "diffB: " << endl;
						cout << diffB << endl;
						cout << "diffB norm: " << diffB.norm() << endl;
					}*/



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
						problem.AddResidualBlock(cost_function,
								   NULL, //if NULL then squared loss
								   CameraModel[cam_id].model,
								   consolidatedLatticeModels[consolidatedGroupID].model,
								   allPoints[(*latticeIt).latticeGridIndices[p3d_id].first].pos.data());
					}
				}
			}
		}
		consolidatedGroupID++;
	}
}

void BundleOptimizer::setupAdvancedPairwiseConsolidatedLatticeOptimizer(){

	ceres::LossFunction* loss_function = FLAGS.robust ? new ceres::HuberLoss(1.50) : NULL;
	ceres::CostFunction* cost_function;

	//create a residual term for each observation of each pair of 3D point on every lattice.

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

				// DEBUGGING OUTPUTS
				/*cout << "===================" << endl;
				cout << lattice1ID << endl;
				cout << lattice2ID << endl;

				cout << "lattice " << lattice1ID << " v0:" << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[0] << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[1] << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[2] << endl;
				cout << endl;
				cout << (*latticeIt1).LattStructure.basisVectors[0] << endl;
				cout << endl;
				cout << "lattice " << lattice1ID << " v1:" << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[3] << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[4] << endl;
				cout << advancedConsolidatedLatticeModels[lattice1ID].model[5] << endl;
				cout << endl;
				cout << (*latticeIt1).LattStructure.basisVectors[1] << endl;
				cout << endl;

				cout << "lattice " << lattice2ID << " v0:" << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[0] << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[1] << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[2] << endl;
				cout << endl;
				cout << (*latticeIt2).LattStructure.basisVectors[0] << endl;
				cout << endl;
				cout << "lattice " << lattice2ID << " v1:" << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[3] << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[4] << endl;
				cout << advancedConsolidatedLatticeModels[lattice2ID].model[5] << endl;
				cout << endl;
				cout << (*latticeIt2).LattStructure.basisVectors[1] << endl;
				cout << endl;*/


				LatticeClass lattice2 = (*latticeIt2);

				int cTransformation1 = (*latticeIt1).consolidationTransformation;
				int cTransformation2 = (*latticeIt2).consolidationTransformation;

				//add cost function (lattice1, lattice2)
				cost_function = VectorDifferenceError::Create(cTransformation1, cTransformation2);

				// *** Debugging
				/*VectorDifferenceError myerr = VectorDifferenceError(cTransformation1, cTransformation2);

				double res[2];

				myerr.operator ()(advancedConsolidatedLatticeModels[lattice1ID].model, advancedConsolidatedLatticeModels[lattice2ID].model, res);

				cout << "res 0: " << res[0] << endl;
				cout << "res 1: " << res[1] << endl;*/
				// ***

				problem.AddResidualBlock(cost_function,
											NULL, //if NULL then squared loss
											advancedConsolidatedLatticeModels[lattice1ID].model,
											advancedConsolidatedLatticeModels[lattice2ID].model);

				lattice2ID++;
			}
			lattice1ID++;
		}

		// Add pairwise cost functions for grid points

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			cout << "*** new lattice" << endl;

			int numLatticeGridPoints = (*latticeIt).latticeGridIndices.size();
			//iterate over all pairs of 3d points: Tp,Tp2

			for (size_t p3d_id=0; p3d_id < numLatticeGridPoints; p3d_id++){ //iteration over 3d points in that lattice

				TriangulatedPoint Tp = allPoints[(*latticeIt).latticeGridIndices[p3d_id].first];
				int a1 = (*latticeIt).latticeGridIndices[p3d_id].second[0];
				int a2 = (*latticeIt).latticeGridIndices[p3d_id].second[1];


				for (size_t p3d_id2=0; p3d_id2 < numLatticeGridPoints; p3d_id2++){ //iteration over 3d points in that lattice

					if (p3d_id == p3d_id2){
						continue;
					}

					TriangulatedPoint Tp2 = allPoints[(*latticeIt).latticeGridIndices[p3d_id2].first];
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
						problem.AddResidualBlock(cost_function,
								   NULL, //if NULL then squared loss
								   CameraModel[cam_id].model,
								   advancedConsolidatedLatticeModels[latticeID].model,
								   allPoints[(*latticeIt).latticeGridIndices[p3d_id].first].pos.data());
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

vector< Eigen::Matrix<double,3,4> > BundleOptimizer::getOptimizedCameras(){

	vector< Eigen::Matrix<double,3,4> > camPoses;
	Eigen::Matrix<double,3,4>* newPose;

	for(size_t i=0; i < number_of_cams; i++){
		newPose = new Eigen::Matrix<double,3,4>(PmatrixFromCeres(CameraModel[i].model));
		camPoses.push_back((*newPose));
	}

	return camPoses;
}

void BundleOptimizer::setLatticeParameters(vector<LatticeClass> &Lattices){

	for (int k=0; k<num_lattices; k++){

		//take a random point and subtract the coordinates

		//Lattices[k].LattStructure.lowerLeftCorner = allPoints[LattModel[k].pivotIndex].pos;

		Lattices[k].LattStructure.basisVectors[0][0] = LattModel[k].model[0];
		Lattices[k].LattStructure.basisVectors[0][1] = LattModel[k].model[1];
		Lattices[k].LattStructure.basisVectors[0][2] = LattModel[k].model[2];
		Lattices[k].LattStructure.basisVectors[1][0] = LattModel[k].model[3];
		Lattices[k].LattStructure.basisVectors[1][1] = LattModel[k].model[4];
		Lattices[k].LattStructure.basisVectors[1][2] = LattModel[k].model[5];

		int a1 = Lattices[k].latticeGridIndices[0].second[0];
		int a2 = Lattices[k].latticeGridIndices[0].second[1];
		Lattices[k].LattStructure.lowerLeftCorner = allPoints[Lattices[k].latticeGridIndices[0].first].pos -
				a1*Lattices[k].LattStructure.basisVectors[0] - a2*Lattices[k].LattStructure.basisVectors[1];

	}
}

void BundleOptimizer::setAdvancedConsolidatedLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices){

	list<list<LatticeClass>>::iterator consolidatedIt;

	int latticeID = 0;

	// iterate over all consolidated groups
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for(latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); latticeIt++){

			(*latticeIt).LattStructure.basisVectors[0].x() = advancedConsolidatedLatticeModels[latticeID].model[0];
			(*latticeIt).LattStructure.basisVectors[0].y() = advancedConsolidatedLatticeModels[latticeID].model[1];
			(*latticeIt).LattStructure.basisVectors[0].z() = advancedConsolidatedLatticeModels[latticeID].model[2];
			(*latticeIt).LattStructure.basisVectors[1].x() = advancedConsolidatedLatticeModels[latticeID].model[3];
			(*latticeIt).LattStructure.basisVectors[1].y() = advancedConsolidatedLatticeModels[latticeID].model[4];
			(*latticeIt).LattStructure.basisVectors[1].z() = advancedConsolidatedLatticeModels[latticeID].model[5];

			int a1 = (*latticeIt).latticeGridIndices[0].second[0];
			int a2 = (*latticeIt).latticeGridIndices[0].second[1];
			(*latticeIt).LattStructure.lowerLeftCorner = allPoints[(*latticeIt).latticeGridIndices[0].first].pos -
						a1*(*latticeIt).LattStructure.basisVectors[0] - a2*(*latticeIt).LattStructure.basisVectors[1];

			latticeID++;
		}
	}
}

void BundleOptimizer::setConsolidatedLatticeParameters(list<list<LatticeClass> > &aConsolidatedLattices){

	int consolidationGroupID = 0;

	list<list<LatticeClass>>::iterator consolidatedIt;

	for (consolidatedIt = aConsolidatedLattices.begin(); consolidatedIt != aConsolidatedLattices.end(); ++consolidatedIt){

		list<LatticeClass>::iterator latticeIt;

		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			// write back the basis vectors with respect to the proper transformation

			switch((*latticeIt).consolidationTransformation){

			case 0:

				(*latticeIt).LattStructure.basisVectors[0][0] = consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = consolidatedLatticeModels[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = consolidatedLatticeModels[consolidationGroupID].model[5]; break;

			case 1:

				(*latticeIt).LattStructure.basisVectors[0][0] = - consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = - consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = - consolidatedLatticeModels[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = consolidatedLatticeModels[consolidationGroupID].model[5]; break;

			case 2:

				(*latticeIt).LattStructure.basisVectors[0][0] = consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = consolidatedLatticeModels[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = - consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = - consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = - consolidatedLatticeModels[consolidationGroupID].model[5]; break;

			case 3:

				(*latticeIt).LattStructure.basisVectors[0][0] = - consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[0][1] = - consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[0][2] = - consolidatedLatticeModels[consolidationGroupID].model[2];
				(*latticeIt).LattStructure.basisVectors[1][0] = - consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[1][1] = - consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[1][2] = - consolidatedLatticeModels[consolidationGroupID].model[5]; break;

			case 4:

				(*latticeIt).LattStructure.basisVectors[0][0] = consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = consolidatedLatticeModels[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = consolidatedLatticeModels[consolidationGroupID].model[2]; break;

			case 5:

				(*latticeIt).LattStructure.basisVectors[0][0] = consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = consolidatedLatticeModels[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = - consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = - consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = - consolidatedLatticeModels[consolidationGroupID].model[2]; break;

			case 6:

				(*latticeIt).LattStructure.basisVectors[0][0] = - consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = - consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = - consolidatedLatticeModels[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = consolidatedLatticeModels[consolidationGroupID].model[2]; break;

			case 7:

				(*latticeIt).LattStructure.basisVectors[0][0] = - consolidatedLatticeModels[consolidationGroupID].model[3];
				(*latticeIt).LattStructure.basisVectors[0][1] = - consolidatedLatticeModels[consolidationGroupID].model[4];
				(*latticeIt).LattStructure.basisVectors[0][2] = - consolidatedLatticeModels[consolidationGroupID].model[5];
				(*latticeIt).LattStructure.basisVectors[1][0] = - consolidatedLatticeModels[consolidationGroupID].model[0];
				(*latticeIt).LattStructure.basisVectors[1][1] = - consolidatedLatticeModels[consolidationGroupID].model[1];
				(*latticeIt).LattStructure.basisVectors[1][2] = - consolidatedLatticeModels[consolidationGroupID].model[2]; break;

			default: break;

			}

		// adjust lowerLeftCorner

		int a1 = (*latticeIt).latticeGridIndices[0].second[0];
		int a2 = (*latticeIt).latticeGridIndices[0].second[1];
		(*latticeIt).LattStructure.lowerLeftCorner = allPoints[(*latticeIt).latticeGridIndices[0].first].pos -
				a1*(*latticeIt).LattStructure.basisVectors[0] - a2*(*latticeIt).LattStructure.basisVectors[1];

		}

		consolidationGroupID++;
	}
}


