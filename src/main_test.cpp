//============================================================================
// Name        : main.cpp
// Author      : group 15
// Version     : 6. April 2016
// Copyright   :
// Description : main function to run 3D vision Project
//============================================================================

#include <opencv2/highgui/highgui.hpp>

#include "detectRepPoints.h"

#include <Eigen/Dense>
#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"

#include <iostream>
#include <fstream>
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
    vector<TriangulatedPoint> pointModel = inpM.getPointModel();
    vector<Eigen::Matrix<double,3,4>> camPoses = inpM.getCamPoses();
    Eigen::Matrix<double,3,3> camK = inpM.getK();

    //vector<vector<Eigen::Vector3d> > groupsOfPoints;
    //groupsOfPoints.push_back(allPoints);

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
    vector<Eigen::Vector3d> max_projectedGroupsOfPoints;
    Eigen::Vector4d maxplane;
    for (int i = 0; i < groupsOfPoints.size(); i++){
    	pf.ransacFit(groupsOfPoints[i]);
    	Eigen::Vector4d bestplane = pf.getFittedPlane();

    	//discard outlying groups
    	if (pf.getProjectedInliers().size() < 10){
    			discarded.push_back(i);
    	}else{
    		//find maximum
    		if (pf.getProjectedInliers().size() > maxid){
    			maxid = i;
    			max_projectedGroupsOfPoints = pf.getProjectedInliers();
    			maxplane = bestplane;
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
    writePlanesToVRML(concatenate(clearedGroupsOfPoints),fittedPlanes,"fitted_planes.wrl", 0.95, true);

    //return 1;


    // -----------------------------------------------------------------------
    // LATTICE FITTING
    // -----------------------------------------------------------------------

    cout << "=Lattice fitting======"<<endl;

    vector<LatticeStructure> lattices;

    /*max_projectedGroupsOfPoints contains the maximal lattice*/
    /*maxplane contains the corresponding plane*/

   // for (int i=0;i < projectedGroupsOfPoints.size(); i++ ){

    	LatticeDetector Ld;
    	//Ld.reconstructedPoints = projectedGroupsOfPoints[i];
    	Ld.reconstructedPoints = max_projectedGroupsOfPoints;
    	//Ld.plane = fittedPlanes[i];
    	Ld.plane = maxplane;
    	Ld.inpManager = &inpM;


    	cout << "will compute basis: " << endl;
        vector<Eigen::Vector3d> candidateBasisVecs = Ld.calculateCandidateVectors(0);

    	cout << "candidate basis vecs: " << endl;

        for (size_t i=0; i< candidateBasisVecs.size();i++){
         	cout << candidateBasisVecs[i] << endl;
        	cout << "---- " << endl;

        }
        vector<Eigen::Vector3d> finalBasisVecs = Ld.getFinalBasisVectors(candidateBasisVecs);

        vector<Vector3d> latticeBoundaries = Ld.calculateLatticeBoundary(finalBasisVecs[0], finalBasisVecs[1]);
        LatticeStructure L;
    	L.plane = maxplane;
//        L.plane = fittedPlanes[i];
        L.basisVectors = finalBasisVecs;

        //testing: comment out
        //Vector3d l1; l1 << 1.73259, -1.61935, -4.49751;
        //Vector3d l2; l2 << 10.7133, 1.22839, -4.12382;
        L.boundary.push_back(latticeBoundaries[0]);
        L.boundary.push_back(latticeBoundaries[1]);
        lattices.push_back(L);

        // ignore: writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, wrlName,true)

    //}

    //vector<LatticeStructure> consolidatedLattices = consolidateLattices(lattices);
    vector<LatticeStructure> consolidatedLattices; consolidatedLattices.push_back(lattices.at(0));
    cout << consolidatedLattices.size() << endl;

    //for the sake of testing max
    vector<vector<Eigen::Vector3d> > maxGroupsOfPoints; maxGroupsOfPoints.push_back(groupsOfPoints[maxid]);
    vector<Eigen::Vector4d>  planes; planes.push_back(maxplane);
    //max testing

    for (int i=0;i < consolidatedLattices.size(); i++ ){
    	LatticeStructure L = lattices[i];
        cout<< L.plane<<endl;
        cout<< L.basisVectors[0]<<endl;
        cout<< L.basisVectors[1]<<endl;
        cout<< L.boundary[0]<<endl;
        cout<< L.boundary[1]<<endl;

    	writeGroupsToVRML(maxGroupsOfPoints,"fitted_latts.wrl", 0.95);
   	    writePlanesToVRML(concatenate(clearedGroupsOfPoints),planes,"fitted_latts.wrl", 0.95, true);
      	writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, "fitted_latts.wrl",true);

    }

    return 0;
}
