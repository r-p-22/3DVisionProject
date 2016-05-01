//============================================================================
// Name        : detectRepPoints.cpp
// Author      : Ryen Elith
// Version     : 6. April 2016
// Copyright   :
// Description : detectRepPoints Class to find repetitive 3D points
//============================================================================

#include "detectRepPoints.h"

// constructor
detectRepPoints::detectRepPoints(char** argv, int computeOrReadArg)
{
    // grouping parameters
    tol_angle = 0.25;               // sift comparison tolerance
    minGroupSize = 8;               // minimum number of members required to form a group
    maxGroupSize = 100;              // if groups are too big, they won't be merged together.

    PCAfilter = false;              // use PCA to filter out less suited groups
    validGroupPCARatio = 0.04;      // min ratio between largest two eigenvalues for valid group
    validGroupPCAEvSize = 1;        // min size of largest eigenvalue of group for valid group

    // save arguments vector locally in class
    classArgv = argv;
    computeOrRead = computeOrReadArg;   // 0: all from file (fastest), 1: recompute grouping, 2: recompute sift descriptors and grouping
    if(computeOrRead == 0)
    {
        readGroups = true;
        readSiftFeatures = true;
    }
    if(computeOrRead == 1)
    {
        readGroups = false;
        readSiftFeatures = true;
    }
    if(computeOrRead == 2)
    {
        readGroups = false;
        readSiftFeatures = false;
    }

    // class internal file names
    file1 = "data/grouping/outputPoints.txt";
    file2 = "data/grouping/outSiftFeaturesVector.txt";
    file3 = "data/grouping/outPointsToTest.txt";
    outputPoints = file1.c_str();
    outputSiftFeatures = file2.c_str();
    outputPointsToTest = file3.c_str();

    // sift feature dimensions
    siftFeatureDim = 128;

    // group organisation variables
    highestGroupIdx = 0;
    countComparisons = 0;
    comparisonsToDo = 0;
    nextGroupIdx.push_back(highestGroupIdx);

    // number of images and image names
    getNumberOfImages();

    // point visibility stuff
    if(!readGroups)    // groups not read from file
    {
        // get neccessary things for grouping
        cout << "getting visibility" << endl;
        get3DPointVisibility();              // reads point file and fills pointsToSift and pointsInImage
        cout << "getting points to test" << endl;
        getPointsToTest();                   // function to fill pointsToTest neede for grouping
        cout << "getting sift representation" << endl;
        pointToGroup = vector<int>(n_points,-1);
        get3DPointSiftRepresentations();     // fills siftFeatureVector
    }
}

// Destructor
detectRepPoints::~detectRepPoints()
{

}

// reads point file and fills pointsToSift (partly) and pointsInImage
int detectRepPoints::get3DPointVisibility()
{
    // toggle console output off
    streambuf *old = cout.rdbuf(0);

    // read 3D point positions
    ifstream is(classArgv[2]);
    if(!is.good())
    {
        cout << "Error finding file" << endl;
        return -1;
    }
    is >> n_points;
    cout << "Total number of points: " << n_points << endl;

    // initialise containers
    pointsToSift = vector<struct siftFeatures>(n_points);       // structure with 3d point information
    for(int i = 0; i<n_points;i++)                                 // binary matrix showing which points are seen in which image
    {
        vector<bool> filler(n_img,0);
        pointsInImage.push_back(filler);
    }

    // fill containers
    for(forLooptype i = 0; i<n_points; i++)
    {
        cout << "Point " << i << " is represented by " << endl;

        // index is trivial
        pointsToSift[i].pointIndex = i;

        // get 3D position
        Eigen::Vector3d pos3D;
        float pos1,pos2,pos3;
        is >> pos1 >> pos2 >> pos3;
        pos3D << pos1,pos2,pos3;

        pointsToSift[i].pos = pos3D;

        // get views
        int nMeasurements = 0;
        is >> nMeasurements;        // point i has nMeasurement sift descriptors

        for (int k = 0; k < nMeasurements; ++k)
        {
            // read remaining 3d point information
            int imIdx, imId;            // image index (view) and image id (not needed)
            Eigen::Vector2f imPos;      // position of projected 3d point in image
            float imPos1,imPos2;
            is >> imIdx >> imId >> imPos1 >> imPos2;
            imPos << imPos1,imPos2;

            // update the pointsToSift vector for given point i
            pointsToSift.at(i).imIndex.push_back(imIdx);
            pointsToSift.at(i).siftPos.push_back(imPos) ;

            // update points in image matrix
            pointsInImage.at(i).at(imIdx) = 1;

            cout << "Sift descriptor " << k << " in image " <<  pointsToSift.at(i).imIndex.back()
                 << " (" << pointsToSift.at(i).siftPos.back()(0) << ","
                 <<  pointsToSift.at(i).siftPos.back()(1) << ")" << endl;
        }
    }

    is.close();

    // toggle cout on
    cout.rdbuf(old);

    return 0;
}

int detectRepPoints::getPointsToTest()
{
    // if file exists, read from file
    ifstream is(outputPointsToTest);
    if(is.good())
    {
        readPointsToTestFromFile(is);
    }
    else
    {
        // binary matrix with comparisons to execute (upper triangle interesting)
        for(int i = 0; i<n_points;i++)
        {
            vector<bool> filler(n_points,0);
            pointsToTest.push_back(filler);
        }

        // check column i with columns i+1 to end
        for (forLooptype i = 0; i<n_points; i++)
        {
            for (forLooptype j = i+1; j<n_points; j++)
            {
                // comparison needed if the two points are seen together in at least one image
                if(bitwiseCompare(pointsInImage[i],pointsInImage[j]))
                {
                    pointsToTest.at(i).at(j) = 1;
                    comparisonsToDo++;
                }
                else
                    pointsToTest.at(i).at(j) = 0;
            }
            // progress of finding points to test
            if(i % ((int)(n_points/(double)100)+1) == 0)
                cout << "progress pointsToTest: " << floor(i/double(n_points)*100) << " %" <<  endl;
        }

        // write results to file
        writePointsToTestToFile();
    }

    return 0;

}

int detectRepPoints::get3DPointSiftRepresentations()
{
    // open file incase reading sift features
    ifstream is(outputSiftFeatures);
    if(!is.good())
    {
        cout << "Problems opening " << outputSiftFeatures << endl;
        return -1;
    }
    int n_img_not_needed;
    is >> n_img_not_needed;

    // get sift features for each point (either compute or read from file)
    for (forLooptype i = 0;i<n_points; i++)
    {
        if(readSiftFeatures)
            cout << "Read sift features for point " << i << endl;
        else
        {
            if(i % ((int)(n_points/(double)100)+1) == 0)
                cout << "Progress getting sift representations: " << floor(i/double(n_points)*100) << " %" <<  endl;
        }

        // toggle console output off
        streambuf *old = cout.rdbuf(0);

        cout << "#############################################" << endl;
        cout << "Point " << i << " has sift features: " << endl;

        // vector for currrent point to add to siftFeatureVector
        vector<Eigen::VectorXf > FeaturesOfOnePoint;

        // get number of views
        int currentViews;
        if(readSiftFeatures)
            is >> currentViews;
        else
            currentViews = pointsToSift[i].imIndex.size();


        for (forLooptype j = 0; j<currentViews; j++) // for each sift feature
        {

            cout << j << ": ";

            Eigen::VectorXf singleFeatureVector(siftFeatureDim);  // filled by function below,assuming descriptor has 128 elements

            if(readSiftFeatures)   // read sift features from file
            {
                for(int h = 0; h<siftFeatureDim; h++)
                {
                    string num;
                    is >> num;
                    int siftElement = atof(num.c_str());
                    singleFeatureVector(h) = siftElement;
                }
            }
            else // recompute
            {
                // collect relevant information for current point
                int image = pointsToSift[i].imIndex.at(j);
                Eigen::Vector2f pos;
                pos << pointsToSift.at(i).siftPos.at(j)(0), pointsToSift.at(i).siftPos.at(j)(1);

                // compute sift vector
                computeSiftDescriptor(image,pos,singleFeatureVector);
            }

            // add new feature 1D-matrix to features of current point
            FeaturesOfOnePoint.push_back(singleFeatureVector);
            cout << FeaturesOfOnePoint.back().transpose() << endl;
        }
        // add all features of this image to siftFeature vector
        siftFeatureVector.push_back(FeaturesOfOnePoint);

        // toggle cout on
        cout.rdbuf(old);
    }
    // close file reading
    is.close();

    // if new featureVector build from image-> store results in txt file
    if(!readSiftFeatures)
    {
        writeSiftFeaturesToFile();
        cout << "Saved point feature vectors to " << outputSiftFeatures << endl;
    }

    return 0;
}

// function to compute siftDescriptor using openCV
int detectRepPoints::computeSiftDescriptor(int imageIndex, Eigen::Vector2f pos, Eigen::VectorXf &outSingleFeatureVector)
{
    cv::Mat B;
    float x = pos(0);
    float y = pos(1);
    bool showRoi = false;
    bool showKeypoints = false;

    const cv::Mat input = cv::imread("data/"+imageNames[imageIndex], 0); //Load as grayscale

    if(! input.data )                              // Check for invalid input
        {
            cout <<  "Could not open or find the image" << std::endl ;
            return -1;
        }
    assert(x>0 && y>0 && x< input.cols && y < input.rows);

    // create mask to search for nearby keypoints to use their size for our keypoint
    cv::Mat mask;
    input.copyTo(mask);
    mask = cv::Scalar(0);
    float maskHalfDim = 30;
    int rectTLx,rectTLy,rectBRx,rectBRy; // for solving boundary issues

    if(x-maskHalfDim<0)
        rectTLx = 0;
    else
        rectTLx = x-maskHalfDim;

    if(y-maskHalfDim<0)
        rectTLy = 0;
    else
        rectTLy = y-maskHalfDim;

    if(x+maskHalfDim>mask.cols)
        rectBRx = mask.cols;
    else
        rectBRx = x+maskHalfDim;

    if(y+maskHalfDim>mask.rows)
        rectBRy = mask.rows;
    else
        rectBRy = y+maskHalfDim;

    cv::Mat roi(mask, cv::Rect(rectTLx,rectTLy,rectBRx-rectTLx,rectBRy-rectTLy));
    roi = cv::Scalar(255);

    // show roi
    if(showRoi)
    {
        cv::Scalar color = cv::Scalar( 100, 100, 100 );
        cv::Point2f center(x, y);
        cv::circle(mask, center, 3,color, 2, 8, 0 );

        cv::namedWindow("Region of interest", CV_WINDOW_NORMAL);
        cv::imshow("Region of interest",mask);
        // cv::waitKey(0);
    }

    // detect sift keypoints in roi
    cv::SiftFeatureDetector detector;
    std::vector<cv::KeyPoint> keypoints;
    detector.detect(input, keypoints, mask);

    // calculate median of keypoint sizes of roi
    float KPsize = 0;
    int n = keypoints.size();
    if(n!=0)
    {
        cout << "Found " << n << " keypoint(s)." << endl;
        vector<float> keypointSizes;
        for(int i =0;i<n;i++)
        {
            keypointSizes.push_back(keypoints[i].size);
            //cout << keypoints[i].size << endl;
        }
        KPsize = median(keypointSizes);
    }
    else
    {
        KPsize = 5;
        cout << "No keypoints detected for roi! " ;
    }

    cout << "Using keypoint size = " << KPsize << endl;

    // Add results to image and save.
    if(showKeypoints)
    {
        cv::Mat output;
        cv::drawKeypoints(input, keypoints, output);
        cv::Scalar boxColor = cv::Scalar(0,255,0);
        cv::rectangle(output,cv::Rect(rectTLx,rectTLy,rectBRx-rectTLx,rectBRy-rectTLy),boxColor,3);
        cv::namedWindow("ROI with keypoints", CV_WINDOW_NORMAL);
        cv::imshow("ROI with keypoints",output);
        cv::waitKey(65);
    }

    // compute descriptor of desired keypoint location using calculated size
    cv::SIFT siftDetector;
    cv::KeyPoint myKeypoint(x, y, KPsize);
    vector<cv::KeyPoint> myKeyPoints;
    myKeyPoints.push_back(myKeypoint);
    cv::Mat descriptor;

    siftDetector(input,mask,myKeyPoints,descriptor, true);

    cout << "Descriptor of point with coordinates (" << x << "," << y << ")" << endl;
    //cout << descriptor << endl;

    // convert descriptor to Eigen::VectorXf (column vector)
    Eigen::VectorXf outDescriptor(descriptor.cols,1);
    for(int i = 0; i<descriptor.cols;i++)
    {
        outDescriptor(i) = descriptor.at<float>(0,i);
    }

    // pass result to specified Eigen::VectorXf
    outSingleFeatureVector = outDescriptor;

    return 0;
}

// function to calculate angle between two descriptors
double detectRepPoints::angleOfTwoSift(Eigen::VectorXf sift1, Eigen::VectorXf sift2)
{
    // check dimensions
    if(sift1.rows()!=sift2.rows())
    {
        cout << "Sift descriptors have wrong dimentions to calculate angle!" << endl;
        return -1;
    }

    // calculate dot product
    double dotProduct = 0.0;
    for (int i = 0; i<sift1.rows();i++)
    {
        dotProduct += sift1(i,0)*sift2(i,0);
    }

    // return the angle (in radians)
    return acos(dotProduct/double(sift1.norm()*sift2.norm()));
}

// calculate median of vector
template<typename T1> T1 detectRepPoints::median(vector<T1> &v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

// compare two 3D points based on their sift descriptors
int detectRepPoints::compare3DPoints(int pointIdx1, int pointIdx2)
{
    // decision result
    int sameGroup = 0;                   // 0: not same group, 1: repetitive -> put in same group

    // get number of sift features for each point
    int n_sift_point1 = pointsToSift[pointIdx1].imIndex.size();
    int n_sift_point2 = pointsToSift[pointIdx2].imIndex.size();

    // compare each sift descriptor of point1 with each sift descriptor of point 2
    for(int i = 0; i < n_sift_point1; i++)  // each sift descriptor of point 1
    {
        //get the descriptor
        Eigen::VectorXf sift1 = siftFeatureVector.at(pointIdx1).at(i);

        for(int j = 0; j < n_sift_point2; j++) // each sift descriptor of point 2
        {

            // get the descriptor
            Eigen::VectorXf sift2 = siftFeatureVector.at(pointIdx2).at(j);

            // if angle between any two sift descriptors is small enough -> classify as repetitive
            if(angleOfTwoSift(sift1,sift2) < tol_angle)
                sameGroup = 1;
        }
    }
    return sameGroup;
}

// output group content to console
void detectRepPoints::coutGroupContent(int groupIdx)
{
    cout << "group "<< groupIdx << " contains points " ;

    for(forLooptype i=0; i<groupToPoints.at(groupIdx).size();i++)
    {
        cout << groupToPoints.at(groupIdx).at(i) << " ";
    }
    cout << endl;
}

// function to add new group for point index
int detectRepPoints::addGroup(int pointIdx)
{
    // add to pointToGroup vector
    int newGroupIdx = nextGroupIdx.back();
    nextGroupIdx.pop_back();
    pointToGroup[pointIdx] = newGroupIdx;

    if(newGroupIdx == highestGroupIdx)      // no recycle group index to use
    {
        cout << "No group recycling ->";

        // add new group
        vector<int> newGroup;
        newGroup.push_back(pointIdx);
        groupToPoints.push_back(newGroup);

        // update recycle vector
        highestGroupIdx++;
        nextGroupIdx.push_back(highestGroupIdx);
    }
    else // recycle group index
    {
        cout << "Recycling group " << newGroupIdx << "->";
        groupToPoints.at(newGroupIdx).push_back(pointIdx);
    }
    coutGroupContent(newGroupIdx);

    return 0;
}

// update group hierarchy given two compared points
void detectRepPoints::updateGroups(int pointIdx1, int pointIdx2,        // compared points
                 int comparisonResult                                   // result of comparison: 0->no match, 1->match
                 )
{

    // case no repetitive points
    if(comparisonResult == 0)
    {
        if(pointToGroup[pointIdx1]==-1) // point1 has no group assigned
        {
            addGroup(pointIdx1);
        }

        if(pointToGroup[pointIdx2]==-1) // point2 has no group assigned
        {
            addGroup(pointIdx2);
        }
    }

    // case repetitive point
    if(comparisonResult == 1)
    {
        // 1) both points have not been assigned to a group -> add both to new group
        if(pointToGroup[pointIdx1]==-1 && pointToGroup[pointIdx2]==-1)
        {
            cout << "grouping case 1" << endl;
            addGroup(pointIdx1);                                                    // new group for point 1
            pointToGroup[pointIdx2] = pointToGroup[pointIdx1];                      // use same group for point 2
            groupToPoints.at(pointToGroup[pointIdx2]).push_back(pointIdx2);
            cout << "New common group generated for point " << pointIdx1 << " and " << pointIdx2 << ": ";
            coutGroupContent(pointToGroup[pointIdx1]);
        }

        // 2) one of the points has been assigned
        else if(pointToGroup[pointIdx1]!=-1 && pointToGroup[pointIdx2]==-1)  // point 1 is assigned to a group
        {
            cout << "grouping case 2" << endl;
            // add point 2 to group of point 1 if group isn't already to big
            if(groupToPoints.at(pointToGroup[pointIdx1]).size() < maxGroupSize)
            {
                pointToGroup[pointIdx2] = pointToGroup[pointIdx1];                                  // use same group for point 2
                groupToPoints.at(pointToGroup[pointIdx2]).push_back(pointIdx2);
                cout << "Point " << pointIdx2 << " added to group of point " << pointIdx1 << ": ";
                coutGroupContent(pointToGroup[pointIdx1]);
            }
            else
                addGroup(pointIdx2);        // otherwise make new group for point 2


        }
        else if(pointToGroup[pointIdx1]==-1 && pointToGroup[pointIdx2]!=-1)  // point 2 is assigned to a group
        {
            // add point 1 to group of point 2 if group isn't already to big
            if(groupToPoints.at(pointToGroup[pointIdx2]).size() < maxGroupSize)
            {
                pointToGroup[pointIdx1] = pointToGroup[pointIdx2];                                  // use same group for point 2
                groupToPoints.at(pointToGroup[pointIdx1]).push_back(pointIdx1);
                cout << "Point " << pointIdx1 << " added to group of " << pointIdx2 << ": ";
                coutGroupContent(pointToGroup[pointIdx2]);
            }
            else
                addGroup(pointIdx1);    // otherwise make new group for point 1
        }

        // 3) both points have been assinged -> merge the groups
        else if(pointToGroup[pointIdx1]!=-1 && pointToGroup[pointIdx2]!=-1)
        {
            // check groups aren't too large, else just don't merge, do nothing
            if(groupToPoints.at(pointToGroup[pointIdx1]).size() < maxGroupSize &&
                    groupToPoints.at(pointToGroup[pointIdx2]).size() < maxGroupSize)
            {
                cout << "grouping case 3: "
                     << "(point " << pointIdx1 << "->group " << pointToGroup[pointIdx1]
                     << ", point " << pointIdx2 << "->group " << pointToGroup[pointIdx2] <<")"
                     << endl;

                // take smaller group index and recycle larger one
                int mergedGroupIdx, recycleGroupIdx;
                if(pointToGroup[pointIdx1]<pointToGroup[pointIdx2])
                {
                    mergedGroupIdx = pointToGroup[pointIdx1];
                    recycleGroupIdx = pointToGroup[pointIdx2];
                }
                else
                {
                    mergedGroupIdx = pointToGroup[pointIdx2];
                    recycleGroupIdx = pointToGroup[pointIdx1];
                }

                // merge groups
                cout << "Merging groups: " << endl;
                coutGroupContent(mergedGroupIdx);
                cout << "with " << endl;
                coutGroupContent(recycleGroupIdx);

                // update point to group vector
                for(forLooptype i = 0; i<pointToGroup.size(); i++)
                {

                    if(pointToGroup[i] == recycleGroupIdx)
                        pointToGroup[i] = mergedGroupIdx;
                }

                // merge groupToPoints vectors
                groupToPoints.at(mergedGroupIdx).
                        insert(groupToPoints.at(mergedGroupIdx).end(),      // append to end of existing group
                               groupToPoints.at(recycleGroupIdx).begin(),  // all elements of the group to be recyled
                               groupToPoints.at(recycleGroupIdx).end());
                groupToPoints.at(recycleGroupIdx).clear();                  // clear elements from group for recycling
                nextGroupIdx.push_back(recycleGroupIdx);                    // use recycled group index for next new group

                cout << "Leading to group" << endl;
                coutGroupContent(mergedGroupIdx);

            }
        } // end of case merging
        else
            cout << "This case should never occur." << endl;
    } // end case repetitive points

}

// main function to get repetitive points
int detectRepPoints::getRepetitivePoints()
{
    int matchResult;
    for (forLooptype i = 0; i<n_points; i++)
    {
        for (forLooptype j = i+1; j<n_points; j++)
        {
            // show progress
            if((countComparisons+1) % ((int)(comparisonsToDo/(double)100)+1) == 0)
                cout << "Progress of grouping: " << floor((countComparisons+1)/double(comparisonsToDo)*100) << " %" <<  endl;

            // toggle console output off
            streambuf *old = cout.rdbuf(0);

            // compare points if needed and {not already in same group, but assiged}
            if(pointsToTest.at(i).at(j)==1 && pointToGroup[i]==pointToGroup[j] && pointToGroup[i]!=-1)
            {
                cout << "Points " << i << " and " << j << " are already in same group." << endl;
                comparisonsToDo--; // one less to do
            }

            if(pointsToTest.at(i).at(j)==1 && (pointToGroup[i]!=pointToGroup[j] || pointToGroup[i]==-1 || pointToGroup[j]==-1))
            {
                cout << "Comparing point " << i << " with point " << j << " ";
                countComparisons++;
                if(compare3DPoints(i, j)==1)
                {
                    cout << "match!" << endl;
                    matchResult = 1;
                }
                else
                {
                    cout << "No match." << endl;
                    matchResult = 0;
                }

                cout << "Updating groups! ";
                updateGroups(i,j,matchResult);
            }

            // toggle cout on
            cout.rdbuf(old);
        }
    } // end of comparisons

    return 0;
}

// function to print group results
int detectRepPoints::printGroupMembers()
{
    if(readGroups)
    {
        // read n_points from file
        ifstream is(classArgv[2]);
        if(!is.good())
        {
            cout << "Error finding points.txt file" << endl;
            return -1;
        }
        is >> n_points;
    }


    // output statistics
    cout << "Statistics: " << endl;
    cout << "Total number of points (" << n_points
         << "), exhaustive would be (" << 0.5*(n_points*n_points-n_points) << ") comparisons." << endl;

    if(computeOrRead == 0)
    {
        for(int i = 0; i<groupsOfPoints.size(); i++)
        {
            cout << "group " << i << " contains " << groupsOfPoints[i].size() << " points." << endl;
        }
        cout << "No further statistics being displayed, because data read from file" << endl;
    }
    else
    {
        cout << "Total number of comparisons executed: " << countComparisons << endl;

        // output groups
        for(forLooptype i = 0; i<groupToPoints.size();i++)
        {
            coutGroupContent(i);
        }

        // output point coordinates
        for(forLooptype i = 0; i<n_points; i++)
        {
            cout << "point "<< i << " has coordinates" << get3dFromPointIdx(i).transpose() << endl;
        }
    }

    return 0;
}

// getter method to get 3d location of given point index
Eigen::Vector3d detectRepPoints::get3dFromPointIdx(int pointIndex)
{
    return pointsToSift[pointIndex].pos;
}

// function to print group results
vector<vector<Eigen::Vector3d> > detectRepPoints::getGroups()
{
    if(computeOrRead == 0)
    {
        readGroupsFromFile();
    }
    else
    {
        // calculate group members indexes
        getRepetitivePoints();

        // build a vector where each element contains a vector of 3d points that belong to that group
        for(forLooptype i = 0; i<groupToPoints.size(); i++)
        {
            // build vector with points of group
            vector<Eigen::Vector3d> currentGroupPoints;
            for(forLooptype j = 0; j<groupToPoints.at(i).size();j++)
            {
                currentGroupPoints.push_back(get3dFromPointIdx(groupToPoints.at(i).at(j)));
            }

            // add current group to groupsOfPoints enough points contained
            if(currentGroupPoints.size() >= minGroupSize)                        // indexes of groupsToPoints != group index (empty groups left out)
            {
                groupsOfPoints.push_back(currentGroupPoints);

                // check PCA constraints -> if not satisfied, remove just added group.
                if(PCAfilter)
                {
                    if(!analyseGroupWithPCA(groupsOfPoints.size()-1))
                        groupsOfPoints.pop_back();
                }
            }
        }

        // retrieve indexing
        computeGroupIndices();

        // save results
        writeGroupsToFile();
    }

    // copy groupOfPoint to vector that is returned
    vector<vector<Eigen::Vector3d> > outGroupsOfPoints;
    outGroupsOfPoints = groupsOfPoints;

    return outGroupsOfPoints;
}

// function to write result to a text file data/outputPoints.txt
int detectRepPoints::writeGroupsToFile()
{
    ofstream of(outputPoints);
    if(!of.good())
    {
        cout << "couldn't open outputfile for Points grouping." << endl;
    }

    // output grouped points to text file
    of << groupsOfPoints.size() << endl;
    for(forLooptype i = 0; i<groupsOfPoints.size();i++)
    {
        of << groupsOfPoints[i].size() << endl;
        for(forLooptype j = 0; j<groupsOfPoints[i].size();j++)
        {
            of << groupsOfPointsIndices.at(i).at(j) << " ";
            of << groupsOfPoints.at(i).at(j).transpose() << endl;
        }
    }
    cout << "Saved groups of points to: " << outputPoints << endl;
    of.close();


    return 0;
}

// read groups from file
int detectRepPoints::readGroupsFromFile()
{
    ifstream is(outputPoints);
    if(!is.good())
    {
        cout << "couldn't open outputfile for Points grouping." << endl;
        return -1;
    }

    // read grouped points from text file
    int numberOfGroups;             // number of groups
    is >> numberOfGroups;
    for(int i = 0; i<numberOfGroups;i++)
    {
        int numberOfPointsInGroup;  // number of points in current group
        is >> numberOfPointsInGroup;
        vector<Eigen::Vector3d> currentGroupPoints;
        vector<int> currentGroupPointsIndices;
        for(int j = 0; j<numberOfPointsInGroup;j++)
        {
            // get index
            int index;
            is >> index;

            // get 3D position
            Eigen::Vector3d pos3D;
            float pos1,pos2,pos3;
            is >> pos1 >> pos2 >> pos3;
            pos3D << pos1,pos2,pos3;

            currentGroupPoints.push_back(pos3D);
            currentGroupPointsIndices.push_back(index);
        }

        // add current group to groupsOfPoints
        groupsOfPoints.push_back(currentGroupPoints);
        groupsOfPointsIndices.push_back(currentGroupPointsIndices);

        // check PCA constraints -> if not satisfied, remove just added group.

        if(PCAfilter)
        {

            if(!analyseGroupWithPCA(groupsOfPoints.size()-1))
            {
                groupsOfPoints.pop_back();
                groupsOfPointsIndices.pop_back();
            }
        }

        cout << "Read group " << i << endl;
    }
    is.close();

    return 0;
}

// function to write siftFeatures results to file data/outputSiftFeatures.txt
int detectRepPoints::writeSiftFeaturesToFile()
{
    ofstream of(outputSiftFeatures);
    if(!of.good())
    {
        cout << "couldn't open outputfile for siftFeature vector." << endl;
    }

    // output grouped points to text file
    of << n_img << endl;        // output number of images
    for(forLooptype i = 0; i<n_points;i++)
    {
        of << siftFeatureVector[i].size() << endl;  // new point
        for(forLooptype j = 0; j<siftFeatureVector[i].size();j++)
        {
            of << siftFeatureVector.at(i).at(j).transpose() << endl;
        }
    }
    return 0;
}

// function to read image names and get number of images
int detectRepPoints::getNumberOfImages()
{
    n_img = 0;

    // check images.txt is present
    ifstream is(classArgv[1]);
    if(!is.good())
    {
        cout << "Error finding images-file: Don't forget to load images folder (with all images) and image.txt file (with all indexes) to data/" << endl;
    }

    // read image names
    while(!is.eof())
    {
        string nextImage;
        is >> nextImage;
        imageNames.push_back(nextImage);
        // cout << nextImage << " referenced." << endl;
        n_img++;
    }
    is.close();

    cout << "Total number of images: " << n_img << endl;        // set by reading images file or siftFeatures file

    return 0;
}

// method to display groups when read from file
int detectRepPoints::visualiseGroups()
{
    if(readGroups)
        get3DPointVisibility();

    for(int i = 0; i<groupsOfPoints.size(); i++)
    {
        // get name of first image where first point of group is visible
        cout << "External group: " << i << endl;

        int imageIndex = pointsToSift.at(groupsOfPointsIndices.at(i).at(0)).imIndex.at(0);
        cout << "Image used for visualisation: " << "data/"+imageNames[imageIndex] << endl;

        // open image for visualisation
        cv::Mat input = cv::imread("data/"+imageNames[imageIndex]); //Load as grayscale
        if(! input.data )                              // Check for invalid input
        {
            cout <<  "Could not open or find the image" << std::endl ;
            return -1;
        }

        // draw circles at points of group
        for(int k = 0; k<groupsOfPoints.at(i).size(); k++)
        {

            // point index
            int pointIdx = groupsOfPointsIndices.at(i).at(k);

            // get position to draw circle
            //check if point is visible in chosen image
            for(int j = 0; j<pointsToSift.at(pointIdx).imIndex.size(); j++)
            {
                if(pointsToSift.at(pointIdx).imIndex.at(j) == imageIndex)
                {
                    double pos_x, pos_y;
                    pos_x = pointsToSift.at(pointIdx).siftPos.at(j)(0);
                    pos_y = pointsToSift.at(pointIdx).siftPos.at(j)(1);
                    cv::Point2f center(pos_x, pos_y);

                    // draw circle
                    cv::circle(input,center,15,cv::Scalar(0,255,0),4);
                }
            }
        }

        // show window
        cv::namedWindow("Visualisation of external group",cv::WINDOW_NORMAL);
        cv::imshow("Visualisation of external group",input);
        cout << "Press key to continue... "<< endl;
        cv::waitKey();
    }

    return 0;
}

// bitwise compare: return true if two binary vectors have value true in same position
bool detectRepPoints::bitwiseCompare(vector<bool> vec1,vector<bool> vec2)
{
    // check vectors have same dimensions
    if(vec1.size() != vec2.size())
    {
        cout<< "Bitwise comparison only possible for equally sized vectors." << endl;
        return false;
    }

    // check if vectors have at least one possition where both are true
    for(int i = 0; i<vec1.size(); i++)
    {
        if(vec1[i] == true && vec2[i] == true)
            return true;
    }

    // if nothing found
    return false;
}

// PCA of group points to see if usefull for fitting lattice
// source: http://docs.opencv.org/3.1.0/d1/dee/tutorial_introduction_to_pca.html#gsc.tab=0
bool detectRepPoints::analyseGroupWithPCA(int externalGroupIdx)
{
    //Construct a buffer used by the pca analysis
    int sz = groupsOfPoints[externalGroupIdx].size();
    cv::Mat data_pts = cv::Mat(sz, 3, CV_64FC1);
    for (int i = 0; i < data_pts.rows; ++i)
    {
        data_pts.at<double>(i, 0) = groupsOfPoints[externalGroupIdx].at(i)(0);
        data_pts.at<double>(i, 1) = groupsOfPoints[externalGroupIdx].at(i)(1);
        data_pts.at<double>(i, 2) = groupsOfPoints[externalGroupIdx].at(i)(2);
    }
    //Perform PCA analysis
    cv::PCA pca_analysis(data_pts, cv::Mat(), CV_PCA_DATA_AS_ROW);

    //Store the center of the object
    //Point cntr = Point(static_cast<int>(pca_analysis.mean.at<double>(0, 0)),
    //                  static_cast<int>(pca_analysis.mean.at<double>(0, 1)));

    //Store the eigenvalues and eigenvectors
    vector<Eigen::VectorXf > eigen_vecs;
    vector<double> eigen_val;

    for (int i = 0; i < 3; ++i)
    {
        // eigenvector
        Eigen::VectorXf eigen_vector(3);
        eigen_vector << pca_analysis.eigenvectors.at<double>(i, 0),
                        pca_analysis.eigenvectors.at<double>(i, 1),
                        pca_analysis.eigenvectors.at<double>(i, 2);
        eigen_vecs.push_back(eigen_vector);

        // eigenvalues
        eigen_val.push_back(pca_analysis.eigenvalues.at<double>(0, i));
        //cout << "group "  << externalGroupIdx <<  ": eigen_val" << i << ": " << eigen_val[i]
        //     << " with eigen_vec: " << eigen_vecs.at(i)(0) << " " << eigen_vecs.at(i)(1) << " " << eigen_vecs.at(i)(2) << endl;
    }

    // find two largest eigenvalues
    sort (eigen_val.begin(), eigen_val.end());
    if(fabs(eigen_val[1]) > validGroupPCARatio* fabs(eigen_val[2])  // ratio constraint
            && fabs(eigen_val[2])>validGroupPCAEvSize)              // largest ev constraint
    {
        // cout << "group " << externalGroupIdx << " is good for eigenval" << endl;

        return true;
    }

    return false;
}

// function to write result to a text file data/groupingoutputPointsToTest.txt
int detectRepPoints::writePointsToTestToFile()
{
    // file to ouput results to file
    ofstream os(outputPointsToTest);
    if(!os.good())
    {
        cout << "Couldn't open file to save pointsToTest." << endl;
        return -1;
    }

    for(int i=0; i<n_points; i++)
    {
        for(int j=0; j<n_points; j++)
        {
            os << pointsToTest.at(i).at(j) << " ";
        }

        // new line
        os << endl;
    }
    os.close();
    cout << "Saved pointsToTest to file " << outputPointsToTest << endl;

    return 0;
}
// function to read pointsToTest from file
int detectRepPoints::readPointsToTestFromFile(ifstream &is)
{
    // file to read results from
    if(!is.good())
    {
        cout << "Couldn't open file to read pointsToTest." << "." << endl;
    }

    vector<bool> rowRead;
    bool dummy;

    for(int i=0; i<n_points; i++)
    {
        for(int j=0; j<n_points; j++)
        {
            is >> dummy;
            rowRead.push_back(dummy);
        }

        // new line
        pointsToTest.push_back(rowRead);
    }
    is.close();

    cout << "Read pointsToTest from file " << outputPointsToTest << "." << endl;

    return 0;
}

// main function function to use to get group indices consisting of 3d points
int detectRepPoints::computeGroupIndices()
{
    // check get groups has been called
    if(groupsOfPoints.size() == 0)
    {
        cout << "Trying to display groups although, groupsOfPoints is empty. Run getGroups() first or adjust grouping parameters (minGroupSize)." << endl;
        return -1;
    }

    // visualise one group after the other
    for(int i = 0; i<groupsOfPoints.size(); i++)
    {
        cout << "Retrieving indexes for points of group " << i << endl;

        // container to save external group member to point indexes
        vector<int> externalGroupPointToIndex;

        // for each point j in current group find visibility index
        for(int j = 0; j<groupsOfPoints[i].size(); j++)
        {
            // search through pointsToSift for right index
            int h = 0;
            while((pointsToSift[h].pos-groupsOfPoints.at(i).at(j)).norm() > 0.0001)
            {
                h++;
                if(h > n_points)
                    cout << "Didn't find index of point in grouping results." << endl;
            }


            // save index h as point j's index in the points.txt file
            externalGroupPointToIndex.push_back(h);
        }

        groupsOfPointsIndices.push_back(externalGroupPointToIndex);
    }

    return 0;
}

// main function function to use to get group indices consisting of 3d points
vector<vector<int> > detectRepPoints::getGroupIndices()
{
    return groupsOfPointsIndices;
}

