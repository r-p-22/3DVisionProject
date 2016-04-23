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
#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"

#include <iostream>
#include <fstream>
#include <string>     // std::string, std::stof
//#include <vector>
#include "my_v3d_vrmlio.h"
#include "3dtools.h"
#include "inputManager.h"

using namespace std;


template<typename T>
vector<T> concatenate(vector<vector<T> > V){
	vector<T> out(V.at(0));
	//out.push_back();
	for (int i = 1; i< V.size(); i++){
		out.insert(out.end(), V.at(i).begin(), V.at(i).end());
	}

	return out;
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

    cout << "will get groups" << endl;
    // compute groups
    groupsOfPoints = myRepPoints.getGroups();

    // print results
    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    // write results to file in grouping folder
    //myRepPoints.writeGroupsToFile();

    ////// Import STUFF


    inputManager inpM(argv);
    vector<Eigen::Vector3d> allPoints = inpM.getPoints();

    // -----------------------------------------------------------------------
    // PLANE FITTING
    // -----------------------------------------------------------------------

    //groupsOfPoints contains the groups

    vector<vector<Eigen::Vector3d> > projectedGroupsOfPoints;
    vector<vector<Eigen::Vector3d> > clearedGroupsOfPoints;

    std::vector<Eigen::Vector4d> fittedPlanes;


    PlaneFitter pf;

    vector<int> discarded;

    int maxid = 0;
    int maxx = -1;
    vector<Eigen::Vector3d> best18_projectedGroupsOfPoints;
    Eigen::Vector4d best18plane;
    for (int i = 0; i < groupsOfPoints.size(); i++){
    	pf.ransacFit(groupsOfPoints[i]);
    	Eigen::Vector4d bestplane = pf.getFittedPlane();

    	//discard outlying groups
    	if (pf.getProjectedInliers().size() < 10){
    			discarded.push_back(i);
    	}else{
    		//get with index 18
    		if (i == 18){
        		//if (pf.getProjectedInliers().size()  > 30){
    			maxid = i;
    			best18_projectedGroupsOfPoints = pf.getProjectedInliers();
    			best18plane = bestplane;
    		}
    		fittedPlanes.push_back(bestplane);
    		projectedGroupsOfPoints.push_back(pf.getProjectedInliers());
    		clearedGroupsOfPoints.push_back(groupsOfPoints[i]);
    	}
    }

    // fittedPlanes inline with clearedGroupsOfPoints and projectedGroupsOfPoints
    // discarded vector is valid only for groupsOfPoints vector

    ////////
    // VISUALIZE GROUPS
    writeGroupsToVRML(clearedGroupsOfPoints,"colored_groups.wrl", 0.95);

    ////////
    // VISUALIZE PLANES
    // The two functions must be called together, to show both points and planes on them
    Eigen::Vector3f color255(255.0,110.,110.);
    writeGroupsToVRML(clearedGroupsOfPoints,"fitted_planes.wrl", 0.95);
    writePlanesToVRML(allPoints,fittedPlanes,"fitted_planes.wrl", 0.95, true);

    vector<vector<Eigen::Vector3d> > maxGroupsOfPoints; maxGroupsOfPoints.push_back(groupsOfPoints[maxid]);
  //  writeGroupsToVRML(maxGroupsOfPoints,"fitted_planes.wrl", 0.95);

  vector<Eigen::Vector4d>  planes; planes.push_back(best18plane);
//    writePlanesToVRML(groupsOfPoints[maxid],planes,"fitted_planes.wrl", 0.99, true);



    //return 1;


    // -----------------------------------------------------------------------
    // LATTICE FITTING
    // -----------------------------------------------------------------------

    cout << "=Lattice fitting======"<<endl;

    vector<LatticeStructure> lattices;

    /*max_projectedGroupsOfPoints contains the maximal lattice*/
    /*maxplane contains the corresponding plane*/
	IOFormat CommaInitFmt(StreamPrecision, DontAlignCols, ", ", ", ", "", "", "", "");

    int loadgroup = 1;
    if (!loadgroup){
		ofstream file("group18.csv");
		for (int i=0; i<best18_projectedGroupsOfPoints.size(); i++){
			file << best18_projectedGroupsOfPoints[i].format(CommaInitFmt)<<endl;
		}

		file.close();
    }
    else {
    	best18_projectedGroupsOfPoints.clear();
    	ifstream ingroup("group18good.csv");
    	Eigen::Vector3d groupvec;
    	string val;
		getline(ingroup, val, ',');
    	while (true){
    		groupvec[0] = stod(val,NULL);
    		getline(ingroup, val, ',');
    		groupvec[1] = stod(val,NULL);
    		getline(ingroup, val, '\n');
    		groupvec[2] = stod(val,NULL);
    		best18_projectedGroupsOfPoints.push_back(groupvec);
    		cout << groupvec << endl;
    		getline(ingroup, val, ',');
    		if (ingroup.eof())
    			break;
    	}
    	ingroup.close();
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

        ofstream file2("candidates18.csv");
           for (int i=0; i<candidateBasisVecs.size(); i++){
               file2 << candidateBasisVecs[i].format(CommaInitFmt)<<endl;
           }
           file2.close();

    	cout << "candidate basis vecs: " << endl;

        for (size_t i=0; i< candidateBasisVecs.size();i++){
         	cout << candidateBasisVecs[i] << endl;
        	cout << "---- " << endl;

        }

        vector<Eigen::Vector3d> finalBasisVecs;
        int loadbasisvecs = 1;
        if (!loadbasisvecs){
        	cout << "calculating final basis" << endl;
			finalBasisVecs = Ld.getFinalBasisVectors(candidateBasisVecs);
			cout << "done!. written to file" << endl;
				ofstream file3("finalbasis18.csv");
				for (int i=0; i<finalBasisVecs.size(); i++){
				   file3 << finalBasisVecs[i].format(CommaInitFmt)<<endl;
				}
				file3.close();

        }else{
    	//If LOAD basis vectors:

        	Eigen::Vector3d a; a << 0,0,0;
			finalBasisVecs.push_back(a);
			finalBasisVecs.push_back(a);
			ifstream inbasisvecs("finalbasis18good.csv");
			string val;
			getline(inbasisvecs, val, ',');
			finalBasisVecs[0][0] = stod(val,NULL);
			getline(inbasisvecs, val, ',');
			finalBasisVecs[0][1] = stod(val,NULL);
			getline(inbasisvecs, val, '\n');
			finalBasisVecs[0][2] = stod(val,NULL);
			getline(inbasisvecs, val, ',');
			finalBasisVecs[1][0] = stod(val,NULL);
			getline(inbasisvecs, val, ',');
			finalBasisVecs[1][1] = stod(val,NULL);
			getline(inbasisvecs, val, '\n');
			finalBasisVecs[1][2] = stod(val,NULL);

			inbasisvecs.close();

        }
		cout << "calculated final bvecs: " << endl;
		cout << finalBasisVecs[0] << endl;
		cout << "---- " << endl;
		cout << finalBasisVecs[1] << endl;
		cout << "---- " << endl;


		vector<Vector3d> latticeBoundaries = Ld.calculateLatticeBoundary(finalBasisVecs[0], finalBasisVecs[1]);
    	cout << "boundary computed.: " << endl;

    	cout << latticeBoundaries[0] << endl;
    	cout << "---- " << endl;

    	cout << latticeBoundaries[1] << endl;

    	ofstream file4("boundaries18.csv");
		   for (int i=0; i<latticeBoundaries.size(); i++){
			   file4 << latticeBoundaries[i].format(CommaInitFmt)<<endl;
		   }
		   file4.close();

        LatticeStructure L;
    	L.plane = best18plane;
//        L.plane = fittedPlanes[i];
        L.basisVectors = finalBasisVecs;

        L.boundary.push_back(latticeBoundaries[0]);
        L.boundary.push_back(latticeBoundaries[1]);
        lattices.push_back(L);

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

        cout<< L.boundary[0]<<endl;
    	cout << "---- " << endl;

        cout<< L.boundary[1]<<endl;
    	cout << "---- " << endl;

    	writeGroupsToVRML(maxGroupsOfPoints,"fitted_latts.wrl", 0.99);
   	    //writePlanesToVRML(concatenate(clearedGroupsOfPoints),planes,"fitted_latts.wrl", 0.95, true);
      	writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, "fitted_latts.wrl",true);

    }

    return 0;
}
