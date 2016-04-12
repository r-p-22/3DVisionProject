//============================================================================
// Name        : main.cpp
// Author      : group 15
// Version     : 6. April 2016
// Copyright   :
// Description : main function to run 3D vision Project
//============================================================================

#include "main.h"
#include "detectRepPoints.h"
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
    detectRepPoints myRepPoints(argv,0);                      // new class, 1: compute from images, 0: take sift features from file
    vector<vector<Eigen::Vector3d> > groupsOfPoints;

    // compute groups
    groupsOfPoints = myRepPoints.getGroups();

    // print results
    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    // write results to file in grouping folder
    myRepPoints.writeGroupsToFile();


    // -----------------------------------------------------------------------
    // PLANE FITTING
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Running plane Fit function" << endl << endl;
    planeFit();

    return 0;
}
