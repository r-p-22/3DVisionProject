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

#include "my_v3d_vrmlio.h"
#include "3dtools.h"
#include "inputManager.h"
using namespace std;

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
    myRepPoints.writeGroupsToFile();


    ////////
    // VISUALIZE GROUPS

    writeGroupsToVRML(groupsOfPoints,"colored_groups.wrl", 0.95);



    //
    //
    // TESTING: ON
    //
    //
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

    std::vector<Eigen::Vector4d> fittedPlanes;


    PlaneFitter pf;
    vector<vector<Eigen::Vector3d> >::iterator it;

    for (it = groupsOfPoints.begin(); it != groupsOfPoints.end(); it++){

    	pf.ransacFit(*it);
    	Eigen::Vector4d bestplane = pf.getFittedPlane();

    	cout << "best:==" << endl;
    	cout << bestplane << endl;
    	cout << "==" << endl;

    	fittedPlanes.push_back(bestplane);
    	projectedGroupsOfPoints.push_back(pf.getProjectedInliers());

    }

    ////////
    // VISUALIZE PLANES
    // The two functions must be called together, to show both points and planes on them
    Eigen::Vector3f color255(255.0,110.,110.);
    writeGroupsToVRML(groupsOfPoints,"fitted_planes.wrl", 0.95);
    writePlanesToVRML(allPoints,fittedPlanes,"fitted_planes.wrl", 0.95, true);


    /* Testing stuff ====>*/
    vector<Eigen::Vector3d> inls = pf.getProjectedInliers();
    cout<< "diff1:==="<<endl;
    cout<< inls[2] - inls[10]<<endl;
    cout<< "diff2:==="<<endl;
    cout<< inls[52] - inls[44]<<endl;

    vector<Eigen::Vector3d> finalBasisVecs;
    finalBasisVecs.push_back(inls[2] - inls[10]);
    finalBasisVecs.push_back(inls[52] - inls[60]);
    /* <=====*/

    // -----------------------------------------------------------------------
    // LATTICE FITTING
    // -----------------------------------------------------------------------


    vector<LatticeStructure> lattices;

    for (int i=0;i < projectedGroupsOfPoints.size(); i++ ){
    	LatticeDetector Ld;
    	Ld.reconstructedPoints = projectedGroupsOfPoints[i];
    	Ld.plane = fittedPlanes[i];
    	Ld.inpManager = &inpM;

    	/*
    	cout << "basis vectors: " << endl;
        vector<Eigen::Vector3d> naiveBasisVecs = Ld.calculateCandidateVectors(1);
        for (size_t i=0; i< naiveBasisVecs.size();i++){
         	cout << naiveBasisVecs[i] << endl;
        }
        vector<Eigen::Vector3d> finalBasisVecs = Ld.getFinalBasisVectors(naiveBasisVecs);
*/
        LatticeStructure L;
        L.plane = fittedPlanes[i];
        L.basisVectors = finalBasisVecs;

        Vector3d l1; l1 << 1.73259, -1.61935, -4.49751;
        Vector3d l2; l2 << 10.7133, 1.22839, -4.12382;
        L.boundary.push_back(l1);
        L.boundary.push_back(l2);
        lattices.push_back(L);

        // ignore: writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, wrlName,true)

    }
    writeGroupsToVRML(groupsOfPoints,"fitted_latts.wrl", 0.95);

    //vector<LatticeStructure> consolidatedLattices = consolidateLattices(lattices);
    vector<LatticeStructure> consolidatedLattices = lattices;

    for (int i=0;i < consolidatedLattices.size(); i++ ){
    	LatticeStructure L = lattices[i];
        cout<< L.plane<<endl;
        cout<< L.basisVectors[0]<<endl;
        cout<< L.basisVectors[1]<<endl;
        cout<< L.boundary[0]<<endl;
        cout<< L.boundary[1]<<endl;

    	writeLatticeToVRML(L.plane,L.basisVectors,L.boundary, "fitted_latts.wrl",true);
    }

    return 0;
}
