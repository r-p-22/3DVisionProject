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
				if (err > 10)
					err = 0;
				reprojectionError += (pa - inpM.getPointModel()[pointidx].measurements[j].pos).norm();
			}
		}

		return reprojectionError;
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
    //myRepPoints.printGroupMembers();

    // write results to file in grouping folder - automatically

    ////// Import STUFF
    inputManager inpM(argv);

	Vector3d col(255.,0.,0.);
	writePointsToVRML(inpM.getPoints(),col,"initial3D.wrl",false);

	// -----------------------------------------------------------------------
	// LATTICE DETECTION
	// -----------------------------------------------------------------------

    int point_count_index = inpM.getPoints().size();

    //done:0,2,4,8,11,17,18
    //13:ok, after increasing int.comb threshold | 33, needed to lower ransac thresh and increase iters
    //27:wrong plane. 2,4: thrown away
    //BAD: 29:no 2 vecs
    int validlattices[] = {0,8 ,11,13,17,31, 33}; //10?14?26?34? //27,31 has 3points | and 18 ofc
    int corrlattices[] = {8}; //10?14?26?34? //27,31 has 3points | and 18 ofc
	vector<LatticeClass> allLattices;

	/*for (size_t v=0;v<1; v++)
	{
		int i = corrlattices[v];
		cout << i << endl;
		string filename = "./data/savedLattices/lattice"+to_string(i)+".txt";
		//LatticeClass mylatt(inpM,groupsOfPoints[i],groupsOfPointsIndices[i],filename.c_str());
	    LatticeClass mylatt(inpM,groupsOfPoints[i],groupsOfPointsIndices[i]);
		//mylatt.projectGroupToImage();
		mylatt.fitLattice();
		mylatt.saveLatticeToFile(filename.c_str());

		mylatt.projectLatticeToImage();

        //point_count_index = mylatt.densifyStructure(point_count_index);

	    //writeQuantilePointsToVRML(inpM.getPoints(),"fitted_latts.wrl", .99);
		//mylatt.writeToVRML("fitted_latts.wrl",true);
	}*/

	for (size_t i=0; i<7; i++)
		{
			int index = validlattices[i];
			cout << index << endl;
			string filename = "./data/savedLattices/lattice"+to_string(index)+".txt";
			LatticeClass lattice(inpM,groupsOfPoints[index],groupsOfPointsIndices[index],filename.c_str());
			allLattices.push_back(lattice);
		}

	/*list<list<LatticeClass> > consolidation = LatticeClass::consolidateLattices(allLattices);

	list<list<LatticeClass> >::iterator consolidationIt;

	for (consolidationIt=consolidation.begin(); consolidationIt != consolidation.end(); consolidationIt++){

		list<LatticeClass>::iterator latticeIt;

		for (latticeIt = (*consolidationIt).begin(); latticeIt != (*consolidationIt).end(); latticeIt++){

			cout << (*latticeIt).LattStructure.basisVectors[0] << endl;
			cout << (*latticeIt).LattStructure.basisVectors[1] << endl;
			cout << (*latticeIt).consolidationTransformation << endl;
			cout << "---" << endl;
			(*latticeIt).projectLatticeToImage();
		}
		cout << "***" << endl;
	}*/

	/*{	int i = 13;
		string filename = "./data/savedLattices/lattice"+to_string(i)+".txt";
		LatticeClass mylatt(inpM,groupsOfPoints[i],groupsOfPointsIndices[i],filename.c_str());
		mylatt.projectGroupToImage();

		//mylatt.fitLattice();
		//mylatt.projectLatticeToImage();

		writePointsToVRML(inpM.getPoints(),col,"latt13.wrl",false);
		mylatt.writeToVRML("latt13.wrl",true);

		allLattices.push_back(mylatt);

	}*/

	int size = allLattices.size();

	cout << "***" << endl;

	float reprojError = 0;
	for (int i=0; i < size; i++){

		//cout << i << endl;
		/*cout << "v0: " << allLattices[i].LattStructure.basisVectors[0] << endl;
		cout << "v1: " << allLattices[i].LattStructure.basisVectors[1] << endl;
		 */
		//allLattices[i].projectLatticeToImage();

		cout << "Reprojection error before bal: ";
		cout << allLattices[i].calculateReprojectionError() << endl;
		//reprojError += allLattices[i].calculateReprojectionError();
	}
	//cout << "Reprojection error before bal: ";
	/*is known to be 60465.6*/
//	cout << calculateReprojectionError(inpM) << endl;

	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------

	BundleOptimizer bal(allLattices, inpM);

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

    return 1;
/*
	list<list<LatticeClass> > consolidatedLattices = LatticeClass::consolidateLattices(allLattices);

	allLattices[0].projectLatticeToImage();

	list<list<LatticeClass> >::iterator consolidatedIt;
	list<LatticeClass>::iterator latticeIt;

	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		cout << "New consolidation class" << endl;
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			cout << "Reprojection error before bal: ";
			cout << (*latticeIt).calculateReprojectionError() << endl;

			cout << "***" << endl;
			cout << "v0: " << (*latticeIt).LattStructure.basisVectors[0] << endl;
			cout << "v1: " << (*latticeIt).LattStructure.basisVectors[1] << endl;

		}
	}


	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------

	BundleOptimizer bal(allLattices, inpM);

	cout << "setting up optimization..." << endl;

	bal.setupNormalOptimizer();
	//bal.setupPairwiseLatticeOptimizer();
	bal.setupPairwiseConsolidatedLatticeOptimizer();

	bal.solve();

	//Need to substitute the new cameras in the inputManager.
	//the inM.pointModel points have already been updated
	inpM.setCamPoses(bal.getOptimizedCameras());

	//bal.setLatticeParameters(allLattices);
	bal.setConsolidatedLatticeParameters(consolidatedLattices);

	writePointsToVRML(inpM.getPoints(),col,"latt13.wrl",false);

	for (consolidatedIt = consolidatedLattices.begin(); consolidatedIt != consolidatedLattices.end(); ++consolidatedIt){
		cout << "New consolidation class" << endl;
		for (latticeIt = (*consolidatedIt).begin(); latticeIt != (*consolidatedIt).end(); ++latticeIt){
			cout << "Reprojection error after bal: ";
			cout << (*latticeIt).calculateReprojectionError() << endl;

			(*latticeIt).projectGroupToImage();
			(*latticeIt).projectLatticeToImage();
			(*latticeIt).writeToVRML("latt13.wrl",true);

			cout << "***" << endl;
			cout << "v0: " << (*latticeIt).LattStructure.basisVectors[0] << endl;
			cout << "v1: " << (*latticeIt).LattStructure.basisVectors[1] << endl;

		}
	}
*/
	return 1;
}


