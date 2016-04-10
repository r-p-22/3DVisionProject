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
    if(argc != 5)
    {
        cout << "Usage: ./3DVisionProject data/images.txt data/points99.txt data/outputPoints.txt data/outSiftFeaturesVector.txt" << endl;
        return -1;
    }

    // -----------------------------------------------------------------------
    // REPETITIVE POINTS
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Computing groups of repetitive points." << endl << endl;
    detectRepPoints myRepPoints(argv);                      // new class
    vector<vector<Eigen::Vector3d> > groupsOfPoints;

    groupsOfPoints = myRepPoints.getGroups();               // compute groups

    // cout << "Statistics and Group members:" << endl;
    // myRepPoints.printGroupMembers();                     // print results

    // workings of groupToPoints (indexing container) and groupOfPoints (vector of 3d points)
    cout << "Groups with member points and their coordinates" << endl;

    int indexForGroupOfPoints = 0;          // needed because groupOfPoints.size() < groupToPoints.size()  (left out empty groups)
    for(int i = 0; i<myRepPoints.groupToPoints.size(); i++)
    {
        if(myRepPoints.groupToPoints.at(i).size()==0)
            cout << "-----" << endl << "Group " << i << " (class internal index) is empty due to recycling and is not part of the groupOfPoints vector." << endl;
        else
        {
            cout << "-----" << endl;
            cout << "Group " << i << " members with their coordinates:" << endl;
            for(int j = 0; j<groupsOfPoints.at(indexForGroupOfPoints).size(); j++)
            {
                cout << "P" << myRepPoints.groupToPoints.at(i).at(j) << ": " << groupsOfPoints.at(indexForGroupOfPoints).at(j).transpose() << endl;
            }
            indexForGroupOfPoints++;
        }
    }

    myRepPoints.writeGroupsToFile(argv[3]);


    // -----------------------------------------------------------------------
    // PLANE FITTING
    // -----------------------------------------------------------------------

    cout << "-----------------------------" << endl;
    cout << "Running plane Fit function" << endl << endl;
    planeFit();

    return 0;
}
