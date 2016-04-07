//============================================================================
// Name        : main.cpp
// Author      : group 15
// Version     : 6. April 2016
// Copyright   :
// Description : main function to run 3D vision Project
//============================================================================

#include "main.h"
#include "detectRepPoints.h"

#include <iostream>
using namespace std;

int main(int argc, char** argv)
{
    // check argc
    if(argc != 3)
    {
        cout << "Usage: ./3DVisionProject data/cams.txt data/points.txt" << endl;
        return -1;
    }


    // new repetitive points class
    cout << "-----------------------------" << endl;
    cout << "Computing repetitive Points." << endl << endl;
    detectRepPoints myRepPoints(argv);
    vector<vector<int> > groupToPoints;
    vector<int> pointToGroup;
    myRepPoints.getRepetitivePoints(pointToGroup,groupToPoints);

    cout << "Statistics and Group members:" << endl;
    myRepPoints.printGroupMembers();

    cout << "Point Coordinates:" << endl;
    for(int i = 0; i<myRepPoints.n_points; i++)
    {
        cout << "point "<< i << " has coordinates" << myRepPoints.get3dFromPointIdx(i).transpose() << endl;
    }

    // run planeFit
    cout << "-----------------------------" << endl;
    cout << "Running plane Fit function" << endl << endl;
    planeFit();

    return 0;
}
