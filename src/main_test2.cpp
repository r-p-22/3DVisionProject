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
//#include "latticeStruct.h"
#include "inputManager.h"
#include "latticeClass.h"

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
    myRepPoints.printGroupMembers();

    // write results to file in grouping folder - automatically

    ////// Import STUFF
    inputManager inpM(argv);

   // myRepPoints.visualiseGroups();

    writeQuantilePointsToVRML(inpM.getPoints(),"group18.wrl", .99);
	Vector3d col(255.,0.,0.);
	cout << groupsOfPoints[18].size() << endl;
	writePointsToVRML(groupsOfPoints[18],col,"group18.wrl",true);

    LatticeClass mylatt(inpM,groupsOfPoints[18],groupsOfPointsIndices[18]);
    mylatt.projectGroupToImage();

    mylatt.fitLattice();

	mylatt.projectLatticeToImage();

	mylatt.writeToVRML("fitted_latts2.wrl",true);

    return 1;

    //==========================================================================================================
    //==========================================================================================================
    //==========================================================================================================
    //==========================================================================================================
    //==========================================================================================================
    //==========================================================================================================
    //==========================================================================================================
    // -----------------------------------------------------------------------
    // PLANE FITTING
    // -----------------------------------------------------------------------
    vector<Eigen::Vector3d> allPoints = inpM.getPoints();

    //groupsOfPoints contains the groups

    vector<vector<Eigen::Vector3d> > projectedGroupsOfPoints;
    vector<vector<Eigen::Vector3d> > validGroupsOfPoints;

    std::vector<Eigen::Vector4d> fittedPlanes;


    PlaneFitter pf;

    vector<int> discarded;

    int maxid = 0;
    int maxx = -1;
    vector<Eigen::Vector3d> best18_projectedGroupsOfPoints;
    vector<int> best18_Indices,ids;
    Eigen::Vector4d best18plane;
    for (int i = 0; i < groupsOfPoints.size(); i++){
    	ids = pf.ransacFit(groupsOfPoints[i],groupsOfPointsIndices[i]);
    	Eigen::Vector4d bestplane = pf.getFittedPlane();

    	//discard outlying groups
    	if (pf.getProjectedInliers().size() < 5){
    			discarded.push_back(i);
    	}else{
    		//get with index 18
    		if (i == 18){
        		//if (pf.getProjectedInliers().size()  > 30){
    			maxid = i;
    			best18_projectedGroupsOfPoints = pf.getProjectedInliers();
    			best18_Indices = ids;
    			best18plane = bestplane;
    	   		validGroupsOfPoints.push_back(groupsOfPoints[i]);

//        		projectedGroupsOfPoints.push_back(pf.getProjectedInliers());
    			break;
    		}
    		fittedPlanes.push_back(bestplane);
    		projectedGroupsOfPoints.push_back(pf.getProjectedInliers());
   		//validGroupsOfPoints.push_back(groupsOfPoints[i]);
    	}
    }

    // fittedPlanes inline with clearedGroupsOfPoints and projectedGroupsOfPoints
    // discarded vector is valid only for groupsOfPoints vector

    ////////
    // VISUALIZE GROUPS
    writeGroupsToVRML(validGroupsOfPoints,"colored_groups.wrl", 0.95);
   // writeGroupsToVRML(projectedGroupsOfPoints,"colored_groups_proj.wrl", 0.95);

    ////////
    // VISUALIZE PLANES
    // The two functions must be called together, to show both points and planes on them
    Eigen::Vector3f color255(255.0,110.,110.);
    writeGroupsToVRML(validGroupsOfPoints,"fitted_planes.wrl", 0.95);
    writePlanesToVRML(allPoints,fittedPlanes,"fitted_planes.wrl", 0.95, true);

    vector<vector<Eigen::Vector3d> > maxGroupsOfPoints; maxGroupsOfPoints.push_back(groupsOfPoints[maxid]);
  //  writeGroupsToVRML(maxGroupsOfPoints,"fitted_planes.wrl", 0.95);

  vector<Eigen::Vector4d>  planes; planes.push_back(best18plane);
//    writePlanesToVRML(groupsOfPoints[maxid],planes,"fitted_planes.wrl", 0.99, true);



  //visualize group18 projection
	//projectGroup(inpM,best18_projectedGroupsOfPoints);
	//projectGroup(inpM,validGroupsOfPoints[0]);
    //return 1;


    // -----------------------------------------------------------------------
    // LATTICE FITTING
    // -----------------------------------------------------------------------

    cout << "=Lattice fitting======"<< endl;

    vector<LatticeStructure> lattices;

    /*max_projectedGroupsOfPoints contains the maximal lattice*/
    /*maxplane contains the corresponding plane*/

    int loadgroup = 1;
    if (!loadgroup){
    	writeToFile(best18_projectedGroupsOfPoints, "group18.csv");
    }
    else {
    	best18_projectedGroupsOfPoints.clear();
    	readFromFile(best18_projectedGroupsOfPoints,3,"group18.csv");
    }
   // for (int i=0;i < projectedGroupsOfPoints.size(); i++ ){
	cout << best18_projectedGroupsOfPoints.size() << endl;

    	LatticeDetector Ld;
    	//Ld.reconstructedPoints = projectedGroupsOfPoints[i];
    	Ld.reconstructedPoints = best18_projectedGroupsOfPoints;
    	//Ld.plane = fittedPlanes[i];
    	Ld.plane = best18plane;
    	Ld.inpManager = &inpM;


    	cout << "will compute basis: " << endl;
        vector<Eigen::Vector3d> candidateBasisVecs = Ld.calculateCandidateVectors(0);

    	writeToFile(best18_projectedGroupsOfPoints, "candidates18.csv");

    	cout << "candidate basis vecs: " << endl;

        for (size_t i=0; i< candidateBasisVecs.size();i++){
         	cout << candidateBasisVecs[i] << endl;
        	cout << "---- " << endl;

        }

        vector<Eigen::Vector3d> finalBasisVecs;
        int loadbasisvecs = 0;
        if (!loadbasisvecs){
        	cout << "calculating final basis" << endl;
			finalBasisVecs = Ld.getFinalBasisVectors(candidateBasisVecs);
	    	writeToFile(finalBasisVecs, "finalbasis18.csv");
        }else{
    	//If LOAD basis vectors:
        	readFromFile(finalBasisVecs,3,"finalbasis18.csv");
        }
		cout << "calculated final bvecs: " << endl;
		cout << finalBasisVecs[0] << endl;
		cout << "---- " << endl;
		cout << finalBasisVecs[1] << endl;
		cout << "---- " << endl;

		int width;
		int height;
		Vector3d lowerLeftCorner;

//-------------------------Boundary computation-----------------------
    	int loadbounds = 0;
    	if (!loadbounds){
    		Ld.calculateLatticeBoundary(finalBasisVecs[0], finalBasisVecs[1], lowerLeftCorner, width, height);
    		ofstream file4("boundaries18.csv");
    		file4 << lowerLeftCorner.format(CommaInitFmt)<<endl;
    		file4 << width << endl;
    		file4 << height << endl;
    		file4.close();
		}
    	else {
    		ifstream file4("boundaries18.csv");
    		string val;
    		getline(file4, val, ',');
    		lowerLeftCorner[0] = stod(val,NULL);
    		getline(file4, val, ',');
			lowerLeftCorner[1] = stod(val,NULL);
			getline(file4, val, '\n');
			lowerLeftCorner[2] = stod(val,NULL);
    		file4 >> width;
    		file4 >> height;
    		file4.close();
    	}
    	cout << "boundary computed.: " << endl;

    	cout << lowerLeftCorner << endl;
    	cout << "---- " << endl;
    	cout << "width: " << width << ". height: " << height << "." << endl;

    	//-------------------------Create Lattice structure-----------------------

        LatticeStructure L;
        //TODO This should maybe be changed back?
    	L.plane = best18plane;
//        L.plane = fittedPlanes[i];
        L.basisVectors = finalBasisVecs;
        L.lowerLeftCorner = lowerLeftCorner;
        L.width = width;
        L.height = height;
        lattices.push_back(L);

        cout << best18_Indices[0] << endl;
        //-------------------------Grid point indices-----------------------
        vector<pair<int, vector<int> > > latticeGridIndices = Ld.getOnGridIndices(best18_Indices, L);

        cout << "***   " << latticeGridIndices[0].first << "   ***   " << latticeGridIndices[0].second[0] << "   ***   " << latticeGridIndices[0].second[1] << endl;

        // ignore: writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, wrlName,true)

    //}

    //vector<LatticeStructure> consolidatedLattices = consolidateLattices(lattices);
    vector<LatticeStructure> consolidatedLattices; consolidatedLattices.push_back(lattices.at(0));

    //for the sake of testing max
   // vector<vector<Eigen::Vector3d> > maxGroupsOfPoints; maxGroupsOfPoints.push_back(groupsOfPoints[maxid]);
   // vector<Eigen::Vector4d>  planes; planes.push_back(best18plane);
    //max testing

    for (int i=0;i < consolidatedLattices.size(); i++ ){
    	LatticeStructure L = lattices[i];
    	cout << "---- " << endl;

        cout<< L.plane<<endl;
    	cout << "---- " << endl;
    	cout << "---- " << endl;

        cout<< L.basisVectors[0]<<endl;
    	cout << "---- " << endl;

        cout<< L.basisVectors[1]<<endl;
    	cout << "---- " << endl;
    	cout << "---- " << endl;

        cout<< L.lowerLeftCorner<<endl;
    	cout << "---- " << endl;

        cout<< L.width<<endl;
    	cout << "---- " << endl;

    	cout<< L.height<<endl;
    	cout << "---- " << endl;

    	writeGroupsToVRML(groupsOfPoints,"fitted_latts.wrl", 0.99);
   	    //writePlanesToVRML(concatenate(clearedGroupsOfPoints),planes,"fitted_latts.wrl", 0.95, true);

    	writeLatticeToVRML(L.plane,L.basisVectors,L.lowerLeftCorner,L.width,L.height, "fitted_latts.wrl",true);

    }

    return 0;
}