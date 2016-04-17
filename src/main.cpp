//============================================================================
// Name        : main.cpp
// Author      : group 15
// Version     : 6. April 2016
// Copyright   :
// Description : main function to run 3D vision Project
//============================================================================

#include "main.h"
#include "detectRepPoints.h"
#include "latticeDetectorTester.h"
#include <opencv2/highgui/highgui.hpp>


#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
    // check argc
    if(argc != 3)
    {
        cout << "Usage: ./3DVisionProject data/images.txt data/points.txt" << endl;
        return -1;
    }

    // -----------------------------------------------------------------------
    // REPETITIVE POINTS
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Computing groups of repetitive points." << endl << endl;
    // new class, choose wheater to recompute of read data from file using second argment:
    // 0: fast:     read groups from file
    // 1: medium:   take sift features from file, recompute groups
    // 2: slow:     recompute sift features and groups

    detectRepPoints myRepPoints(argv,1);
    vector<vector<Eigen::Vector3d> > groupsOfPoints;

    // compute groups
    groupsOfPoints = myRepPoints.getGroups();

    // print results
    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    // -----------------------------------------------------------------------
    // PLANE FITTING
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Running plane Fit function" << endl << endl;
    //planeFit();

    LatticeDetectorTester myLatticeDetectorTester = LatticeDetectorTester();
    myLatticeDetectorTester.test();

    return 0;
}
