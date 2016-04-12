//============================================================================
// Name        : detectRepPoints.h
// Author      : Ryen Elith
// Version     : 6. April 2016
// Copyright   :
// Description : Headerfile for Class to find repetitive 3D points
//============================================================================
#ifndef DETECTREPPOINTS_H
#define DETECTREPPOINTS_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/nonfree/features2d.hpp> //Thanks to Alessandro
#include <string>

using namespace std;

class detectRepPoints{

private:

        // output filenames
        string file1, file2;
        const char *outputPoints, *outputSiftFeatures;

        // generic data info
        int siftFeatureDim;                                                // dimension: 128 for sift

        // type def for for loops
        typedef unsigned long forLooptype;

        // copy of main argument vector
        char** classArgv;

        // point visibility stuff

            // vector holding image names
            vector<string> imageNames;

            // class own container with images each point is seen in
            Eigen::MatrixXf pointsInImage;

            // points to compare against eachother (2 points that can be seen in same image - any image)
            Eigen::MatrixXf pointsToTest;

            // function to get point visibilities in images (reads point file and fills pointsToSift (partly) and pointsInImage)
            int get3DPointVisibility();

            // builds binary matrix pointsToTest with comparisons to execute
            int getPointsToTest();


        // Sift stuff

            // choice wheather to calculate sift from images or take from fiel
            bool computeFromImages;

            // struct to store 3D point information
            struct siftFeatures {
                int pointIndex;                     // point index -> trivial, equal to position of struct in pointsToSift vector
                Eigen::Vector3d pos;                // point position in 3D
                vector<int> imIndex;                // image indexes where point is seen
                vector<Eigen::Vector2f> siftPos;    // position of sift features in corresponding image
            };

            // container storing: lookup point -> sift descriptors
            vector<vector<Eigen::MatrixXf> > siftFeatureVector;

            // container storing: 3d point -> point's information (siftFeature struct)
            vector<struct siftFeatures> pointsToSift;

            // function to write siftFeatures results to file data/outputSiftFeatures.txt
            int writeSiftFeaturesToFile();

            // function to calculate angle between two descriptors
            double angleOfTwoSift(Eigen::MatrixXf sift1, Eigen::MatrixXf sift2);

            // calculate median of vector
            template<typename T1> T1 median(vector<T1> &v);

            // function to computes all sift descriptors and fill siftFeatureVector and complete pointsToSift structs (with sift indexes)
            int get3DPointSiftRepresentations();


            // function to compute siftDescriptor of one image using openCV
            int computeSiftDescriptor(int imageIndex,Eigen::Vector2f pos,Eigen::MatrixXf &outSingleFeatureVector);

        // grouping stuff

            // group organisation variables
            int highestGroupIdx;                                        // highest group index in use
            vector<int> nextGroupIdx;                                   // next group index to assign for new groups (holds recycled groups)
            double tol_angle;                                           // decision criteria angle for repetitive points
            int countComparisons;                                       // count of comparisons executed to find groups

            // function to compare two 3D points based on their sift descriptors
            int compare3DPoints(int pointIdx1, int pointIdx2);

            // function to add new group for point index
            int addGroup(int pointIdx);

            // update group hierarchy given two compared points: result of comparison: 0->no match, 1->match
            void updateGroups(int pointIdx1, int pointIdx2, int comparisonResult);

            // function to ouput content of group given a group index
            void coutGroupContent(int groupIdx);

            // main function to get repetitive points indexes
            int getRepetitivePoints();

            // getter method to get 3d location of given point index
            Eigen::Vector3d get3dFromPointIdx(int pointIndex);

            // vector with grouped 3d points (no recycled groups included)
            vector<vector<Eigen::Vector3d> > groupsOfPoints;


public:
        // constructor
        detectRepPoints(char** classArgv, bool computeFromImagesArg);
        // destructor
        ~detectRepPoints();

        // number of points and images
        int n_img;
        forLooptype n_points;

        //matrix pointToGroup with 1 to 1 relation: point index -> group index
        vector<int> pointToGroup;

        // 2d vector with 1 to n relation: group index -> point index
        vector<vector<int> > groupToPoints;

        // function to print grouping results (as indexes) with some statistics
        int printGroupMembers();

        // main function function to use to get groups consisting of 3d points
        vector<vector<Eigen::Vector3d> > getGroups();

        // function to write result to a text file data/outputPoints.txt
        int writeGroupsToFile();


};


#endif // detectRepPoints_H
