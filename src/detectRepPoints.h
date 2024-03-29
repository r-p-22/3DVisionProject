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
        string file1, file2, file3;
        const char *outputPoints, *outputSiftFeatures, *outputPointsToTest;

        // generic data info
        int siftFeatureDim;                                                // dimension: 128 for sift

        // type def for for loops
        typedef unsigned long forLooptype;

        // copy of main argument vector
        char** classArgv;

        // reading files stuff
        bool readSiftFeatures;
        bool readGroups;

        // point visibility stuff

            // vector holding image names
            vector<string> imageNames;

            // class own container with images each point is seen in
            vector<vector<bool> > pointsInImage;

            // points to compare against eachother (2 points that can be seen in same image - any image)
            vector<vector<bool> > pointsToTest;

            // struct to store 3D point information
            struct siftFeatures {
                int pointIndex;                     // point index -> trivial, equal to position of struct in pointsToSift vector
                Eigen::Vector3d pos;                // point position in 3D
                vector<int> imIndex;                // image indexes where point is seen
                vector<Eigen::Vector2f> siftPos;    // position of sift features in corresponding image
            };

            // function to read image names and get number of images
            int getNumberOfImages();

            // function to get point visibilities in images (reads point file and fills pointsToSift and pointsInImage)
            int get3DPointVisibility();

            // builds binary matrix pointsToTest with comparisons to execute
            int getPointsToTest();

        // Sift stuff

            // choice wheather to calculate sift from images or take from fiel
            int computeOrRead;

            // container storing: lookup point -> sift descriptors
            vector<vector<Eigen::VectorXf> > siftFeatureVector;

            // container storing: 3d point -> point's information (siftFeature struct)
            vector<struct siftFeatures> pointsToSift;

            // function to write siftFeatures results to file data/outputSiftFeatures.txt
            int writeSiftFeaturesToFile();

            // function to calculate angle between two descriptors
            double angleOfTwoSift(Eigen::VectorXf sift1, Eigen::VectorXf sift2);

            // calculate median of vector
            template<typename T1> T1 median(vector<T1> &v);

            // function to computes all sift descriptors and fill siftFeatureVector
            int get3DPointSiftRepresentations();

        // grouping stuff

            // group organisation variables
            int highestGroupIdx;                                        // highest group index in use
            vector<int> nextGroupIdx;                                   // next group index to assign for new groups (holds recycled groups)
            double tol_angle;                                           // decision criteria angle for repetitive points
            int countComparisons;                                       // count of comparisons executed to find groups
            int comparisonsToDo;                                        // number of comparisons to execute
            int minGroupSize;                                           // minimum number of points needed to form a group
            int maxGroupSize;                                           // upper bound size of groups where merging still possible
            double validGroupPCA;                                       // min ratio between largest two eigenvalues for valid group
            double validGroupPCARatio;                                  // min ratio between largest two eigenvalues for valid group
            double validGroupPCAEvSize;                                 // min size of largest eigenvalue of group for valid group
            bool PCAfilter;                                             // toggle PCA filtering for groups on/off

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

            // method to read groupOfPoints from file
            int readGroupsFromFile();

            // function to write result to a text file data/grouping/outputPoints.txt
            int writeGroupsToFile();

            // functions to read and write pointToTest from/to file
            int writePointsToTestToFile();
            int readPointsToTestFromFile(ifstream &is);

            // vector with grouped 3d points (no recycled groups included) and indices
            vector<vector<Eigen::Vector3d> > groupsOfPoints;
            vector<vector<int> > groupsOfPointsIndices;

            //matrix pointToGroup with 1 to 1 relation: point index -> group index
            vector<int> pointToGroup;

            // 2d vector with 1 to n relation: group index -> point index
            vector<vector<int> > groupToPoints;

            // bitwise compare: return true if two binary vectors have value true in same position
            bool bitwiseCompare(vector<bool> vec1,vector<bool> vec2);

            // PCA of group points to see if usefull for fitting lattice
            bool analyseGroupWithPCA(int externalGroupIdx);


public:
        // constructor
        detectRepPoints(char** classArgv, int computeOrReadArg);
        // destructor
        ~detectRepPoints();

        // number of points and images
        int n_img;
        forLooptype n_points;

        // function to compute siftDescriptor of one image using openCV
        int computeSiftDescriptor(int imageIndex, Eigen::Vector2f pos, Eigen::VectorXf &outSingleFeatureVector);

        // main function function to use to get groups consisting of 3d points
        vector<vector<Eigen::Vector3d> > getGroups();

        // main function function to use to get group indices consisting of 3d points
        vector<vector<int> > getGroupIndices();

        // function to print grouping results (as indexes) with some statistics
        int printGroupMembers();

        // function to visualise grouping results
        int visualiseGroups();

};


#endif // detectRepPoints_H
