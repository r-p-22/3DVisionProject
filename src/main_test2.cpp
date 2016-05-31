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
#include <string>     // std::string, std::stof
//#include <vector>
//#include "my_v3d_vrmlio.h"
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
        cout << "Usage: ./3DVisionProject images.txt points.txt cams.txt K.txt" << endl;
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
    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    // write results to file in grouping folder - automatically

    ////// Import STUFF
    inputManager inpM(argv);

	Vector3d col(255.,0.,0.);
	writePointsToVRML(inpM.getPoints(),col,"initial3D.wrl",false);

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

		cout << i << endl;
		int v = validlattices[i];
		string filename = "./data/savedLattices/lattice"+to_string(v)+".txt";
		LatticeClass mylatt(inpM,groupsOfPoints[v],groupsOfPointsIndices[v],filename.c_str());
	    //LatticeClass mylatt(inpM,groupsOfPoints[i],groupsOfPointsIndices[i]);
		//mylatt.fitLattice();
		//mylatt.saveLatticeToFile(filename.c_str());

		//mylatt.projectLatticeToImage();
		//mylatt.projectGroupToImage();

        //point_count_index = mylatt.densifyStructure(point_count_index);

	    //writeQuantilePointsToVRML(inpM.getPoints(),"fitted_latts.wrl", .99);
		//mylatt.writeToVRML("fitted_latts.wrl",true);
		allLattices.push_back(mylatt);

	}

	int size = allLattices.size();

	/*
	cout << "***" << endl;

	float reprojError = 0;
	for (int i=0; i < size; i++){

		//cout << i << endl;
		//cout << "v0: " << allLattices[i].LattStructure.basisVectors[0] << endl;
		//cout << "v1: " << allLattices[i].LattStructure.basisVectors[1] << endl;

		//allLattices[i].projectLatticeToImage();

		cout << "Reprojection error before bal: ";
		cout << allLattices[i].calculateReprojectionError() << endl;
		//reprojError += allLattices[i].calculateReprojectionError();
	}*/

	//cout << "Reprojection error before bal: ";
	//cout << calculateReprojectionError(inpM) << endl;

	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------

	/*BundleOptimizer bal(allLattices, inpM);

	cout << "setting up optimization..." << endl;

	//bal.setupAllNormalOptimizer();
	bal.setupNormalOptimizer();
	bal.setupPairwiseLatticeOptimizer();

	bal.solve();

	//Need to substitute the new cameras in the inputManager.
	//the inM.pointModel points have already been updated
	inpM.setCamPoses(bal.getOptimizedCameras());

	bal.setLatticeParameters(allLattices);

	cout << "***" << endl;
	writePointsToVRML(inpM.getPoints(),col,"latticesPairwise.wrl",false);

	for ( int i=0; i < size; i++){

		//cout << i << endl;
		//cout << "v0: " << allLattices[i].LattStructure.basisVectors[0] << endl;
		//cout << "v1: " << allLattices[i].LattStructure.basisVectors[1] << endl;

		cout << "Reprojection error after bal: ";
		cout << allLattices[i].calculateReprojectionError() << endl;
		//reprojError += allLattices[i].calculateReprojectionError();

		allLattices[i].writeToVRML("latticesPairwise.wrl",true);

		//allLattices[i].projectGroupToImage();
		//allLattices[i].projectLatticeToImage();
	}

	cout << "Reprojection error after bal: ";
	cout << calculateReprojectionError(inpM) << endl;

    return 1;*/

	list<list<LatticeClass> > consolidatedLattices = LatticeClass::consolidateLattices(allLattices);

	list<list<LatticeClass> >::iterator consolidatedIt;
	list<LatticeClass>::iterator latticeIt;

	/*for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			//cout << "Reprojection error before bal: ";
			//cout << (*latticeIt).calculateReprojectionError() << endl;
		}
	}*/

	/*for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){

			cout << (*latticeIt).latticeGridIndices.size() << endl;
			cout << (*latticeIt).LattStructure.basisVectors[0] << endl;
			cout << (*latticeIt).LattStructure.basisVectors[1] << endl;
			//(*latticeIt).projectGroupToImage();
			//(*latticeIt).projectLatticeToImage();
		}
	}




	return 1;*/

	/*vector<Vector3d> allModelPoints = inpM.getPoints();

	vector<Vector3d> widthVecs;
		vector<Vector3d> heightVecs;

		Vector3d vec_2686_to_1020 = (allModelPoints[1020]-allModelPoints[2686])*0.5;
		Vector3d vec_1020_to_1002 = allModelPoints[1002]-allModelPoints[1020];
		Vector3d vec_614_to_599 = allModelPoints[599]-allModelPoints[614];
		widthVecs.push_back(vec_2686_to_1020);
		widthVecs.push_back(vec_1020_to_1002);
		widthVecs.push_back(vec_614_to_599);

		Vector3d vec_1020_to_2186 = allModelPoints[2186]-allModelPoints[1020];
		Vector3d vec_614_to_555 = allModelPoints[555]-allModelPoints[614];
		Vector3d vec_1241_to_195 = allModelPoints[195]-allModelPoints[1241];
		heightVecs.push_back(vec_1020_to_2186);
		heightVecs.push_back(vec_614_to_555);
		heightVecs.push_back(vec_1241_to_195);

		for (int i=0; i < widthVecs.size(); i++){
			cout << widthVecs[i] << endl;
			cout << "***" << endl;
		}

		for (int i=0; i < heightVecs.size(); i++){
			cout << heightVecs[i] << endl;
			cout << "***" << endl;
		}

		return 1;*/

	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------

	BundleOptimizer bal(consolidatedLattices, inpM);

	cout << "setting up optimization..." << endl;

	double gridTransformationWeight = 0;
	double basisVectorWeight = 0;

	//bal.setupStandardOptimizer();
	bal.setupConsolidatedLatticeOptimizer(gridTransformationWeight,basisVectorWeight);
	//bal.setupRigidConsolidatedLatticeOptimizer(0.01);

	outputCostsInConsole(bal);

	bal.solve();
	cout << "optimization done" << endl;

	inpM.updatePointsWithModel();
	inpM.setCamPoses(bal.getOptimizedCameras());

	bal.readoutLatticeParameters(consolidatedLattices);
	//bal.readoutRigidLatticeParameters(consolidatedLattices);

	outputCostsInConsole(bal);

	vector<Vector3d> allModelPoints = inpM.getPoints();

	outputDistanceVectors("./data/distanceVectors/width_vectors_"+to_string(gridTransformationWeight)+"_"+to_string(basisVectorWeight)+".txt", allModelPoints, true);
	outputDistanceVectors("./data/distanceVectors/height_vectors_"+to_string(gridTransformationWeight)+"_"+to_string(basisVectorWeight)+".txt", allModelPoints, false);


	/*vector<Vector3d> widthVecs;
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

	for (int i=0; i < widthVecs.size(); i++){
		cout << widthVecs[i] << endl;
		cout << "***" << endl;
	}

	for (int i=0; i < heightVecs.size(); i++){
		cout << heightVecs[i] << endl;
		cout << "***" << endl;
	}

	vector<Vector3d> standardWidthVecs;
	vector<Vector3d> standardHeightVecs;
	standardWidthVecs.push_back(Vector3d(0.441192,0.00210911,0.0975631));
	standardWidthVecs.push_back(Vector3d(0.441767,0.00620341,0.106381));
	standardWidthVecs.push_back(Vector3d(0.44441,0.00683026,0.103692));
	standardWidthVecs.push_back(Vector3d(0.442898,0.0040882,0.0931621));
	standardHeightVecs.push_back(Vector3d(-0.00319963,-0.103002,0.0252757));
	standardHeightVecs.push_back(Vector3d(-0.00518085,-0.105306,0.0222232));
	standardHeightVecs.push_back(Vector3d(-0.00551161,-0.104468,0.0250939));



	Vector3d colors[2];
	colors[0] = Vector3d(255.,0.,0.);
	colors[1] = Vector3d(0.,0.,255.);

	writePointsToVRML(widthVecs, colors[1],"widthVecs.wrl",false);
	writePointsToVRML(heightVecs,colors[1],"heightVecs.wrl",false);
	writePointsToVRML(standardWidthVecs,colors[0],"widthVecs.wrl",true);
	writePointsToVRML(standardHeightVecs,colors[0],"heightVecs.wrl",true);*/

	//writeQuantilePointsToVRML(allModelPoints,col,"test.wrl",false,0.95);

	// LATTICE0
	/*
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

	// get average
	Vector3d averageWidthVec = Vector3d(0,0,0);
	for (int i=0; i < widthVecs.size(); i++){
		averageWidthVec = averageWidthVec + widthVecs[i];
	}
	averageWidthVec = averageWidthVec/((double)widthVecs.size());

	// normalize with average
	for (int i=0; i < widthVecs.size(); i++){
		widthVecs[i] = widthVecs[i]/averageWidthVec.norm();
	}
	averageWidthVec = averageWidthVec/averageWidthVec.norm();

	// subtract average
	for (int i=0; i < widthVecs.size(); i++){
		widthVecs[i] = widthVecs[i]-averageWidthVec;
	}

	// get average
		Vector3d averageHeightVec = Vector3d(0,0,0);
		for (int i=0; i < heightVecs.size(); i++){
			averageHeightVec = averageHeightVec + heightVecs[i];
		}
		averageHeightVec = averageHeightVec/((double)heightVecs.size());

		// normalize with average
		for (int i=0; i < heightVecs.size(); i++){
			heightVecs[i] = heightVecs[i]/averageHeightVec.norm();
		}
		averageHeightVec = averageHeightVec/averageHeightVec.norm();

		// subtract average
		for (int i=0; i < heightVecs.size(); i++){
			heightVecs[i] = heightVecs[i]-averageHeightVec;
		}



	for (int i=0; i < widthVecs.size(); i++){
		cout << widthVecs[i] << endl;
		cout << "***" << endl;
	}

	for (int i=0; i < heightVecs.size(); i++){
		cout << heightVecs[i] << endl;
		cout << "***" << endl;
	}

	vector<Vector3d> standardWidthVecs;
	vector<Vector3d> standardHeightVecs;
	standardWidthVecs.push_back(Vector3d(0.441192,0.00210911,0.0975631));
	standardWidthVecs.push_back(Vector3d(0.441767,0.00620341,0.106381));
	standardWidthVecs.push_back(Vector3d(0.44441,0.00683026,0.103692));
	standardWidthVecs.push_back(Vector3d(0.442898,0.0040882,0.0931621));
	standardHeightVecs.push_back(Vector3d(-0.00319963,-0.103002,0.0252757));
	standardHeightVecs.push_back(Vector3d(-0.00518085,-0.105306,0.0222232));
	standardHeightVecs.push_back(Vector3d(-0.00551161,-0.104468,0.0250939));

	// get average
		Vector3d standardAverageWidthVec = Vector3d(0,0,0);
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardAverageWidthVec = standardAverageWidthVec + standardWidthVecs[i];
		}
		standardAverageWidthVec = standardAverageWidthVec/((double)standardWidthVecs.size());

		// normalize with average
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardWidthVecs[i] = standardWidthVecs[i]/standardAverageWidthVec.norm();
		}
		standardAverageWidthVec = standardAverageWidthVec/standardAverageWidthVec.norm();

		// subtract average
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardWidthVecs[i] = standardWidthVecs[i]-standardAverageWidthVec;
		}

		// get average
				Vector3d standardAverageHeightVec = Vector3d(0,0,0);
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardAverageHeightVec = standardAverageHeightVec + standardHeightVecs[i];
				}
				standardAverageHeightVec = standardAverageHeightVec/((double)standardHeightVecs.size());

				// normalize with average
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardHeightVecs[i] = standardHeightVecs[i]/standardAverageHeightVec.norm();
				}
				standardAverageHeightVec = standardAverageHeightVec/standardAverageHeightVec.norm();

				// subtract average
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardHeightVecs[i] = standardHeightVecs[i]-standardAverageHeightVec;
				}

	Vector3d colors[7];
	colors[0] = Vector3d(255.,0.,0.);
	colors[1] = Vector3d(175.,175.,0.);
	colors[2] = Vector3d(0.,0.,255.);
	colors[3] = Vector3d(175.,0.,175.);
	colors[4] = Vector3d(0.,255.,175.);
	colors[5] = Vector3d(0.,175.,255.);
	colors[6] = Vector3d(175.,0.,255.);

	//writePointsToVRML(standardWidthVecs,colors[0],"widthVecs.wrl",false);
	writePointsToVRML(standardHeightVecs,colors[0],"heightVecs.wrl",false);
	//writePointsToVRML(widthVecs, colors[2],"widthVecs.wrl",true);
	writePointsToVRML(heightVecs,colors[2],"heightVecs.wrl",true);
*/

/*	vector<Vector3d> widthVecs;
	vector<Vector3d> heightVecs;

	Vector3d vec_2686_to_1020 = (allModelPoints[1020]-allModelPoints[2686])*0.5;
	Vector3d vec_1020_to_1002 = allModelPoints[1002]-allModelPoints[1020];
	Vector3d vec_614_to_599 = allModelPoints[599]-allModelPoints[614];
	widthVecs.push_back(vec_2686_to_1020);
	widthVecs.push_back(vec_1020_to_1002);
	widthVecs.push_back(vec_614_to_599);

	Vector3d vec_1020_to_2186 = allModelPoints[2186]-allModelPoints[1020];
	Vector3d vec_614_to_555 = allModelPoints[555]-allModelPoints[614];
	Vector3d vec_1241_to_195 = allModelPoints[195]-allModelPoints[1241];
	heightVecs.push_back(vec_1020_to_2186);
	heightVecs.push_back(vec_614_to_555);
	heightVecs.push_back(vec_1241_to_195);



	// get average
	Vector3d averageWidthVec = Vector3d(0,0,0);
	for (int i=0; i < widthVecs.size(); i++){
		averageWidthVec = averageWidthVec + widthVecs[i];
	}
	averageWidthVec = averageWidthVec/((double)widthVecs.size());

	// normalize with average
	for (int i=0; i < widthVecs.size(); i++){
		widthVecs[i] = widthVecs[i]/averageWidthVec.norm();
	}
	averageWidthVec = averageWidthVec/averageWidthVec.norm();

	// subtract average
	for (int i=0; i < widthVecs.size(); i++){
		widthVecs[i] = widthVecs[i]-averageWidthVec;
	}

	// get average
		Vector3d averageHeightVec = Vector3d(0,0,0);
		for (int i=0; i < heightVecs.size(); i++){
			averageHeightVec = averageHeightVec + heightVecs[i];
		}
		averageHeightVec = averageHeightVec/((double)heightVecs.size());

		// normalize with average
		for (int i=0; i < heightVecs.size(); i++){
			heightVecs[i] = heightVecs[i]/averageHeightVec.norm();
		}
		averageHeightVec = averageHeightVec/averageHeightVec.norm();

		// subtract average
		for (int i=0; i < heightVecs.size(); i++){
			heightVecs[i] = heightVecs[i]-averageHeightVec;
		}

	for (int i=0; i < widthVecs.size(); i++){
		cout << widthVecs[i] << endl;
		cout << "***" << endl;
	}

	for (int i=0; i < heightVecs.size(); i++){
		cout << heightVecs[i] << endl;
		cout << "***" << endl;
	}

	vector<Vector3d> standardWidthVecs;
	vector<Vector3d> standardHeightVecs;
	standardWidthVecs.push_back(Vector3d(0.44683,0.000258562,0.0801904));
	standardWidthVecs.push_back(Vector3d(0.442024,0.00392044,0.0803062));
	standardWidthVecs.push_back(Vector3d(0.443965,0.0013403,0.100773));
	standardHeightVecs.push_back(Vector3d(-0.0247044,-0.562359,0.122301));
	standardHeightVecs.push_back(Vector3d(-0.0370549,-0.540935,0.127764));
	standardHeightVecs.push_back(Vector3d(-0.0252952,-0.560397,0.188671));

	// get average
		Vector3d standardAverageWidthVec = Vector3d(0,0,0);
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardAverageWidthVec = standardAverageWidthVec + standardWidthVecs[i];
		}
		standardAverageWidthVec = standardAverageWidthVec/((double)standardWidthVecs.size());

		// normalize with average
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardWidthVecs[i] = standardWidthVecs[i]/standardAverageWidthVec.norm();
		}
		standardAverageWidthVec = standardAverageWidthVec/standardAverageWidthVec.norm();

		// subtract average
		for (int i=0; i < standardWidthVecs.size(); i++){
			standardWidthVecs[i] = standardWidthVecs[i]-standardAverageWidthVec;
		}

		// get average
				Vector3d standardAverageHeightVec = Vector3d(0,0,0);
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardAverageHeightVec = standardAverageHeightVec + standardHeightVecs[i];
				}
				standardAverageHeightVec = standardAverageHeightVec/((double)standardHeightVecs.size());

				// normalize with average
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardHeightVecs[i] = standardHeightVecs[i]/standardAverageHeightVec.norm();
				}
				standardAverageHeightVec = standardAverageHeightVec/standardAverageHeightVec.norm();

				// subtract average
				for (int i=0; i < standardHeightVecs.size(); i++){
					standardHeightVecs[i] = standardHeightVecs[i]-standardAverageHeightVec;
				}

	Vector3d colors[7];
	colors[0] = Vector3d(255.,0.,0.);
	colors[1] = Vector3d(175.,175.,0.);
	colors[2] = Vector3d(0.,0.,255.);
	colors[3] = Vector3d(175.,0.,175.);
	colors[4] = Vector3d(0.,255.,175.);
	colors[5] = Vector3d(0.,175.,255.);
	colors[6] = Vector3d(175.,0.,255.);

	writePointsToVRML(standardWidthVecs,colors[0],"consWidthVecs_1_100000000.wrl",false);
	writePointsToVRML(standardHeightVecs,colors[0],"consHeightVecs_1_100000000.wrl",false);
	writePointsToVRML(widthVecs, colors[2],"consWidthVecs_1_100000000.wrl",true);
	writePointsToVRML(heightVecs,colors[2],"consHeightVecs_1_100000000.wrl",true);*/






	/*vector<vector<Vector3d> > allLatticePoints = vector<vector<Vector3d> >();
	vector<Vector3d> allModelPoints = inpM.getPoints();

	//Vector3d firstVec = allModelPoints[0];
	vector<pair<int,vector<int> > >::iterator pointsIt;
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			vector<Vector3d> pointsVec;
			for(pointsIt = (*latticeIt).latticeGridIndices.begin(); pointsIt != (*latticeIt).latticeGridIndices.end(); ++pointsIt){
				int index = (*pointsIt).first;
				pointsVec.push_back(allModelPoints[index]);
				cout << allModelPoints[index] << endl;
				cout << "---" << endl;
				//allModelPoints[index] = firstVec;
			}
			cout << "***" << endl;
			allLatticePoints.push_back(pointsVec);
			//cout << (*latticeIt).LattStructure.lowerLeftCorner << endl;
		}
	}

	Vector3d colors[7];
	colors[0] = Vector3d(255.,255.,0.);
	colors[1] = Vector3d(0.,255.,255.);
	colors[2] = Vector3d(255.,0.,255.);
	colors[3] = Vector3d(175.,255.,0.);
	colors[4] = Vector3d(0.,255.,175.);
	colors[5] = Vector3d(0.,175.,255.);
	colors[6] = Vector3d(175.,0.,255.);

	//writePointsToVRML(allModelPoints,col,"test.wrl",false);

	vector<vector<Vector3d> >::iterator allLatticePointsIt;

	int colorID = 0;
	bool append = false;
	for (allLatticePointsIt = allLatticePoints.begin(); allLatticePointsIt != allLatticePoints.end(); ++allLatticePointsIt){
		writePointsToVRML((*allLatticePointsIt),colors[colorID],"rigidOptimizedLattices_001.wrl",append);
		append = true;
		colorID++;
	}

	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			(*latticeIt).writeToVRML("rigidOptimizedLattices_001.wrl",true);
		}
	}*/



	//bal.readoutLatticeParameters(consolidatedLattices);

	//bal.setLatticeParameters(allLattices);
	//bal.setConsolidatedLatticeParameters(consolidatedLattices);
	//bal.readoutAdvancedConsolidatedLatticeParameters(consolidatedLattices);

	//vector<Vector3d> allModelPoints = inpM.getPoints();

	//writeQuantilePointsToVRML(inpM.getPoints(),col,"test.wrl", 0.95, false);

	/*vector<vector<Vector3d> > allLatticePoints = vector<vector<Vector3d> >();
	vector<Vector3d> allModelPoints = inpM.getPoints();

	//Vector3d firstVec = allModelPoints[0];
	vector<pair<int,vector<int> > >::iterator pointsIt;
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			vector<Vector3d> pointsVec;
			for(pointsIt = (*latticeIt).latticeGridIndices.begin(); pointsIt != (*latticeIt).latticeGridIndices.end(); ++pointsIt){
				int index = (*pointsIt).first;
				pointsVec.push_back(allModelPoints[index]);
				cout << allModelPoints[index] << endl;
				cout << "---" << endl;
				//allModelPoints[index] = firstVec;
			}
			cout << "***" << endl;
			pointsVec.push_back((*latticeIt).LattStructure.lowerLeftCorner);
			allLatticePoints.push_back(pointsVec);
			//cout << (*latticeIt).LattStructure.lowerLeftCorner << endl;
		}
	}

	Vector3d colors[7];
	colors[0] = Vector3d(255.,255.,0.);
	colors[1] = Vector3d(0.,255.,255.);
	colors[2] = Vector3d(255.,0.,255.);
	colors[3] = Vector3d(175.,255.,0.);
	colors[4] = Vector3d(0.,255.,175.);
	colors[5] = Vector3d(0.,175.,255.);
	colors[6] = Vector3d(175.,0.,255.);

	//writePointsToVRML(allModelPoints,col,"test.wrl",false);

	vector<vector<Vector3d> >::iterator allLatticePointsIt;

	int colorID = 0;
	bool append = true;
	for (allLatticePointsIt = allLatticePoints.begin(); allLatticePointsIt != allLatticePoints.end(); ++allLatticePointsIt){
		writePointsToVRML((*allLatticePointsIt),colors[colorID],"test.wrl",append);
		append = true;
		colorID++;
	}

	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			(*latticeIt).writeToVRML("test.wrl",true);
		}
	}*/


	/*writePointsToVRML(inpM.getPoints(),col,"latticeAdvConsol.wrl",false);

	cout << "total reproj error: " << calculateReprojectionError(inpM) << endl;

	float groupReprojError = 0;
	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		cout << "New consolidation class" << endl;
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			cout << "Reprojection error after bal: ";
			//float err = (*latticeIt).calculateReprojectionError();
			cout << err << endl;
			groupReprojError += err;
			//(*latticeIt).projectGroupToImage();
			//(*latticeIt).projectLatticeToImage();
			//(*latticeIt).writeToVRML("latticeAdvConsol.wrl",true);

			//cout << "***" << endl;
			//cout << "v0: " << (*latticeIt).LattStructure.basisVectors[0] << endl;
			//cout << "v1: " << (*latticeIt).LattStructure.basisVectors[1] << endl;

		}
	}

	cout << "summed group reproj error: " << groupReprojError << endl;*/

	return 1;
}


