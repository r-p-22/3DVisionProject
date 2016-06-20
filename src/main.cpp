//============================================================================
// Name        : main.cpp
// Author      : group 15
// Version     : 6. April 2016
// Copyright   :
// Description : main function to run 3D vision Project
//============================================================================

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include "detectRepPoints.h"

#include <Eigen/Dense>

#include <iostream>
#include <fstream>
#include <string>
#include "3dtools.h"
#include "inputManager.h"
#include "latticeClass.h"

#include "BundleOptimizer.h"

using namespace std;

IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");

template<typename T>
vector<T> concatenate(vector<vector<T> > V){
	vector<T> out(V.at(0));
	//out.push_back();
	for (int i = 1; i< V.size(); i++){
		out.insert(out.end(), V.at(i).begin(), V.at(i).end());
	}

	return out;
}

template<typename T>
void writeToFile(vector<T> vec, const char* filename){
	ofstream file3(filename);
	for (int i=0; i<vec.size(); i++){
		file3 << vec[i].format(CommaInitFmt)<<endl;
	}
	file3.close();
}



template<typename T>
void readFromFile(vector<T>& resvec, int fieldsize, const char* filename){
	ifstream infile(filename);
	Eigen::VectorXd datavec(fieldsize);
	string val;
	getline(infile, val, ',');
	while (true){

		datavec[0] = stod(val,NULL);
		for(int j=1;j<fieldsize-1;j++)
		{
			getline(infile, val, ',');
			datavec[j] = stod(val,NULL);
		}
		getline(infile, val, '\n');
		datavec[fieldsize-1] = stod(val,NULL);

		resvec.push_back(datavec);

		cout << datavec << endl;

		getline(infile, val, ',');

		if (infile.eof())
			break;
	}
	infile.close();

}


float calculateReprojectionError(inputManager inpM){
	Vector2f pa;
	CameraMatrix cam;
	cam.setIntrinsic(inpM.getK());

	float reprojectionError = 0;

	for (int pointidx=0; pointidx < inpM.getPoints().size(); pointidx++){

		for (int j=0; j<inpM.getPointModel()[pointidx].measurements.size();j++){
			int imgview = inpM.getPointModel()[pointidx].measurements[j].view;

			int i=0;
			for (i=0;i<inpM.getCamPoses().size();i++){
				if (inpM.getViewIds()[i] == imgview)
					break;
			}
			Eigen::Matrix<double,3,4> P = inpM.getCamPoses()[i];
			cam.setOrientation(P);
			pa = cam.projectPoint(inpM.getPointModel()[pointidx].pos).cast<float>();

			float err = (pa - inpM.getPointModel()[pointidx].measurements[j].pos).norm();
			//if (err > 50)
			//	err = 50;
			reprojectionError += err;
		}
	}

	return reprojectionError;
}

void outputDistanceVectors(string filename, const vector<Vector3d> &allModelPoints, bool outputWidth){

	vector<Vector3d> widthVecs;
	vector<Vector3d> heightVecs;

	Vector3d vec_1_0_to_2_0 = allModelPoints[1264]-allModelPoints[1050];
	Vector3d vec_2_0_to_3_0 = allModelPoints[858]-allModelPoints[1264];
	Vector3d vec_3_0_to_3_1 = allModelPoints[493]-allModelPoints[858];
	Vector3d vec_3_1_to_3_2 = allModelPoints[1370]-allModelPoints[493];
	Vector3d vec_2_2_to_3_2 = allModelPoints[1370]-allModelPoints[673];
	Vector3d vec_2_0_to_2_2 = (allModelPoints[673]-allModelPoints[1264])*0.5;
	Vector3d vec_0_2_to_2_2 = (allModelPoints[673]-allModelPoints[950])*0.5;

	widthVecs.push_back(vec_1_0_to_2_0);
	widthVecs.push_back(vec_2_0_to_3_0);
	widthVecs.push_back(vec_2_2_to_3_2);
	widthVecs.push_back(vec_0_2_to_2_2);

	heightVecs.push_back(vec_3_0_to_3_1);
	heightVecs.push_back(vec_3_1_to_3_2);
	heightVecs.push_back(vec_2_0_to_2_2);

	ofstream os(filename);

	if(outputWidth){
		for (int i = 0; i < widthVecs.size(); i++){
			os << widthVecs[i][0] << "," << widthVecs[i][1] << "," << widthVecs[i][2] << endl;
		}
	}
	else{
		for (int i = 0; i < heightVecs.size(); i++){
			os << heightVecs[i][0] << "," << heightVecs[i][1] << "," << heightVecs[i][2] << endl;
		}
	}

	os.close();
}

void outputCostsInConsole(BundleOptimizer &bal){

	double totalCosts, pointReprojectionCosts, gridTransformationCosts, basisVectorCosts;

	totalCosts = bal.calculateCost(BundleOptimizer::CostType::TOTAL);
	pointReprojectionCosts = bal.calculateCost(BundleOptimizer::CostType::POINT_REPROJECTION);
	gridTransformationCosts = bal.calculateCost(BundleOptimizer::CostType::GRID_TRANSFORMATION);
	basisVectorCosts = bal.calculateCost(BundleOptimizer::CostType::BASIS_VECTORS);

	cout << "Total costs: " << totalCosts << endl;

	// as control, should be equal to total costs
	cout << "Summed costs: " << (pointReprojectionCosts + gridTransformationCosts + basisVectorCosts) << endl;

	cout << "Points costs: " << pointReprojectionCosts << endl;
	cout << "Grid transformation costs: " << gridTransformationCosts << endl;
	cout << "Basis vector costs: " << basisVectorCosts << endl;
}

int main(int argc, char** argv)
{
	cv::initModule_nonfree();
    // check argc
    if(argc != 5)
    {
        cout << "Usage: ./latt_bal images.txt points.txt cams.txt K.txt" << endl;
        return -1;
    }

    // -----------------------------------------------------------------------
    // REPETITIVE POINTS
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Computing groups of repetitive points." << endl << endl;
    detectRepPoints myRepPoints(argv,0);                      // new class, 1: compute from images, 0: take sift features from file
    vector<vector<Eigen::Vector3d> > groupsOfPoints;
    vector<vector<int> > groupsOfPointsIndices;

    cout << "will get groups" << endl;
    // compute groups
    groupsOfPoints = myRepPoints.getGroups();
    groupsOfPointsIndices = myRepPoints.getGroupIndices();

    // print results
    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    // write results to file in grouping folder - automatically

    ////// Import STUFF
    inputManager inpM(argv);

	// -----------------------------------------------------------------------
	// LATTICE DETECTION
	// -----------------------------------------------------------------------

    int point_count_index = inpM.getPoints().size();

    int validlattices[] = {0,8 ,11,13,17,31, 33}; //10?14?26?34? //27,31 has 3points | and 18 ofc

	vector<LatticeClass> allLattices;
	cout << "importing lattices" << endl;

	//points.txt
	// 1st and 2nd group: Too big
	//
	//for (size_t v=0;v<groupsOfPoints.size(); v++)
	for (size_t i=0;i<7; i++) {

		int v = validlattices[i];
		string filename = "./data/savedLattices/lattice"+to_string(v)+".txt";
		LatticeClass mylatt(inpM,groupsOfPoints[v],groupsOfPointsIndices[v],filename.c_str());
		allLattices.push_back(mylatt);

	}

	int size = allLattices.size();

	list<list<LatticeClass> > consolidatedLattices = LatticeClass::consolidateLattices(allLattices);

	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------


	vector<Vector3d> allModelPointsBefore = inpM.getPoints();

	outputDistanceVectors("./data/distanceVectors/width_vectors_before_BA.txt", allModelPointsBefore, true);
	outputDistanceVectors("./data/distanceVectors/height_vectors_before_BA.txt", allModelPointsBefore, false);


	BundleOptimizer bal(consolidatedLattices, inpM);

	cout << "setting up optimization..." << endl;

	double gridTransformationWeight = 10;
	double basisVectorWeight = 100;

	//bal.setupStandardOptimizer();
	bal.setupConsolidatedLatticeOptimizer(gridTransformationWeight,basisVectorWeight);
	//bal.setupRigidConsolidatedLatticeOptimizer(0.01);

	outputCostsInConsole(bal);

	bal.solve();

	cout << "optimization done" << endl;

	// read out optimized values
	inpM.updatePointsWithModel();
	inpM.setCamPoses(bal.getOptimizedCameras());
	bal.readoutLatticeParameters(consolidatedLattices);
	//bal.readoutRigidLatticeParameters(consolidatedLattices);

	outputCostsInConsole(bal);

	vector<Vector3d> allModelPoints = inpM.getPoints();

	outputDistanceVectors("./data/distanceVectors/width_vectors_"+to_string(gridTransformationWeight)+"_"+to_string(basisVectorWeight)+".txt", allModelPoints, true);
	outputDistanceVectors("./data/distanceVectors/height_vectors_"+to_string(gridTransformationWeight)+"_"+to_string(basisVectorWeight)+".txt", allModelPoints, false);

	return 1;
}


