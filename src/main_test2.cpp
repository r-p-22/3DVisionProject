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
	writePointsToVRML(groupsOfPoints.at(13),col,"group18.wrl",false);

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

	list<list<LatticeClass> > consolidation = LatticeClass::consolidateLattices(allLattices);

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
	}

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


	/*
	// LATTICE CONSOLIDATION TEST

			vector<LatticeClass> lattices = vector<LatticeClass>();

			LatticeClass lattice1 = LatticeClass(Vector3d(0,-1,0), Vector3d(1,0,0), Vector4d(1,1,1,1));
			LatticeClass lattice2 = LatticeClass(Vector3d(0,-1,0), Vector3d(-1,0,0), Vector4d(1,1,1,1));
			LatticeClass lattice3 = LatticeClass(Vector3d(0,1,0), Vector3d(-1,0,0), Vector4d(1,1,1,1));
			LatticeClass lattice4 = LatticeClass(Vector3d(0,1.15,0), Vector3d(1.15,0,0), Vector4d(1,1,1,1));
			LatticeClass lattice5 = LatticeClass(Vector3d(1.15,0,0), Vector3d(0,-1.15,0), Vector4d(1,1,1,1));
			LatticeClass lattice6 = LatticeClass(Vector3d(-1.15,0,0), Vector3d(0,-1.15,0), Vector4d(1,1,1,1));
			LatticeClass lattice7 = LatticeClass(Vector3d(-1,0,0), Vector3d(0,1,0), Vector4d(1,1,1,1));
			LatticeClass lattice8 = LatticeClass(Vector3d(1,0,0), Vector3d(0,1,0), Vector4d(1,1,1,1));

			lattices.push_back(lattice1);
			lattices.push_back(lattice2);
			lattices.push_back(lattice6);
			lattices.push_back(lattice4);
			lattices.push_back(lattice7);
			lattices.push_back(lattice8);
			lattices.push_back(lattice3);
			lattices.push_back(lattice5);


			list<list<LatticeClass> > consolidation = LatticeClass::consolidateLattices(lattices);

			list<list<LatticeClass> >::iterator consolidationIt;

			for (consolidationIt=consolidation.begin(); consolidationIt != consolidation.end(); consolidationIt++){
				list<LatticeClass>::iterator latticeIt;
				for (latticeIt = (*consolidationIt).begin(); latticeIt != (*consolidationIt).end(); latticeIt++){
					cout << (*latticeIt).LattStructure.basisVectors[0] << endl;
					cout << (*latticeIt).LattStructure.basisVectors[1] << endl;
					cout << (*latticeIt).consolidationTransformation << endl;
					cout << "---" << endl;
				}
				cout << "***" << endl;
			}
		 */

	/*allLattices[0].projectLatticeToImage();

	cout << "Reprojection error before bal: ";
	cout << allLattices[0].calculateReprojectionError() << endl;


	// -----------------------------------------------------------------------
	// BUNDLE ADJUSTMENT OPTIMIZATION
	// -----------------------------------------------------------------------

	BundleOptimizer bal(allLattices, inpM);

	cout << "setting up optimization..." << endl;

	bal.setupNormalOptimizer();
	bal.setupPairwiseLatticeOptimizer();

	bal.solve();

	//Need to substitute the new cameras in the inputManager.
	//the inM.pointModel points have already been updated
	inpM.setCamPoses(bal.getOptimizedCameras());

	bal.setLatticeParameters(allLattices);

	cout << "Reprojection error after bal: ";
	cout << allLattices[0].calculateReprojectionError() << endl;

	allLattices[0].projectGroupToImage();
	allLattices[0].projectLatticeToImage();*/

    return 1;

}


