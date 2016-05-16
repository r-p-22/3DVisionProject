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
#include "my_v3d_vrmlio.h"
#include "3dtools.h"
#include "inputManager.h"
#include "latticeClass.h"

//#include "BundleOptimizer.h"

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
	writePointsToVRML(groupsOfPoints.at(18),col,"group18.wrl",false);

	// -----------------------------------------------------------------------
	// LATTICE DETECTION
	// -----------------------------------------------------------------------

	vector<LatticeClass> allLattices;
    int point_count_index = inpM.getPoints().size();
	//for (size_t i=0;i<groupsOfPoints.size(); i++)
	{
		int i = 18;
		string filename = "./data/savedLattices/lattice"+to_string(i)+".txt";
		LatticeClass mylatt(inpM,groupsOfPoints[i],groupsOfPointsIndices[i],filename.c_str());
		mylatt.projectGroupToImage();
	    //mylatt.fitLattice();
		//mylatt.saveLatticeToFile(filename.c_str());
		mylatt.projectLatticeToImage();

        point_count_index = mylatt.densifyStructure(point_count_index);

	    //writeQuantilePointsToVRML(inpM.getPoints(),"fitted_latts.wrl", .99);
		//mylatt.writeToVRML("fitted_latts.wrl",true);
		allLattices.push_back(mylatt);
	}


		//TODO Is this line still needed?
		cout << allLattices[0].lattStructure.plane << endl;

	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------
/*
	BundleOptimizer balNormal(allLattices, inpM);

	cout << "setting up optimization..." << endl;

	balNormal.setupNormalOptimizer();

	cout << inpM.getPointModel().size() << endl;

	balNormal.solve();

	//Need to substitute the new cameras in the inputManager.
	//the inM.pointModel points have already been updated
	inpM.setCamPoses(balNormal.getOptimizedCameras());

	cout << "projecting..." << endl;
	allLattices[0].projectGroupToImage();
*/
    return 1;

}


