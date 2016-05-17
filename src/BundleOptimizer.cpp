/*
 * BundleOptimizer.cpp
 *
 *  Created on: May 14, 2016
 *      Author: lianos91
 */

#include "BundleOptimizer.h"

#include "ceres/ceres.h"


struct find_pivot {
    vector<int> v{0,0};
    find_pivot(vector<int> v) : v(v) {}
    find_pivot(){}
    bool operator () ( const pair<int, vector<int> >& m ) const
    {
        return ((m.second[0] == v[0])&&(m.second[1]==v[1]));
    }
};

BundleOptimizer::BundleOptimizer(vector<LatticeClass> &Lattices, inputManager &inpM) {

	allPoints = inpM.getPointModel();
	camViewIndices = inpM.getViewIds();

	number_of_cams = inpM.getCamPoses().size();
	num_lattices = Lattices.size();

	this->Lattices = vector<LatticeClass>(Lattices);

	focalLength = inpM.getK()(0,0);
	principalPoint[0] = inpM.getK()(0,2);
	principalPoint[1] = inpM.getK()(1,2);

	CameraModel = new CeresCamModel[number_of_cams];


	/* Contents of LattModel.model[9]:
	 * 		0-2:   Lowerleft position 3D -- DEPRECATED
	 * 		3-5:  basisVector1 3D
	 * 		6-8: basisVector2 3D
	 */
	LattModel = new CeresLattModel[num_lattices];
	for (size_t k=0; k<num_lattices; k++){

		//find pivot point and set as lower left
		vector<pair<int, vector<int> >>::iterator pivot = std::find_if( Lattices[k].latticeGridIndices.begin(),
					Lattices[k].latticeGridIndices.end(), find_pivot());
		int pivotIndex = pivot->first;
		cout << "pivotIndex: " << pivotIndex << endl;
		TriangulatedPoint Tp_pivot = allPoints[pivotIndex];

		LattModel[k].pivotIndex = pivotIndex;

		LattModel[k].model[3] = Lattices[k].LattStructure.basisVectors[0].x();
		LattModel[k].model[4] = Lattices[k].LattStructure.basisVectors[0].y();
		LattModel[k].model[5] = Lattices[k].LattStructure.basisVectors[0].z();
		LattModel[k].model[6] = Lattices[k].LattStructure.basisVectors[1].x();
		LattModel[k].model[7] = Lattices[k].LattStructure.basisVectors[1].y();
		LattModel[k].model[8] = Lattices[k].LattStructure.basisVectors[1].z();
	}

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

	ceres::LossFunction* loss_function = NULL;//FLAGS.robust ? new ceres::HuberLoss(2.0) : NULL;

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

void BundleOptimizer::setupGridOptimizer(){

	ceres::LossFunction* loss_function = FLAGS.robust ? new ceres::HuberLoss(.50) : NULL;
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
						   NULL, //if NULL then squared loss
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

}


void BundleOptimizer::solve() {

	ceres::Solver::Options options;
	ceres::Solver::Summary summary;

	options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;
	options.minimizer_progress_to_stdout = true;
	options.max_num_iterations = 200;
	options.function_tolerance = 1.e-12;
	options.parameter_tolerance = 1.e-12;
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

		Lattices[k].LattStructure.lowerLeftCorner = allPoints[LattModel[k].pivotIndex].pos;

		Lattices[k].LattStructure.basisVectors[0][0] = LattModel[k].model[3];
		Lattices[k].LattStructure.basisVectors[0][1] = LattModel[k].model[4];
		Lattices[k].LattStructure.basisVectors[0][2] = LattModel[k].model[5];
		Lattices[k].LattStructure.basisVectors[1][0] = LattModel[k].model[6];
		Lattices[k].LattStructure.basisVectors[1][1] = LattModel[k].model[7];
		Lattices[k].LattStructure.basisVectors[1][2] = LattModel[k].model[8];
	}


}


