//============================================================================
// Name        : detectRepPoints.cpp
// Author      : Ryen Elith
// Version     : 6. April 2016
// Copyright   :
// Description : detectRepPoints Class to find repetitive 3D points
//============================================================================

#include "detectRepPoints.h"

// constructor
detectRepPoints::detectRepPoints(char** argv, bool computeFromImagesArg)
{
    // save arguments vector locally in class
    classArgv = argv;
    computeFromImages = computeFromImagesArg;

    // number of images
    n_img=0;

    // sift feature dimensions
    siftFeatureDim = 128;

    ifstream is(classArgv[1]);
    if(!is.good())
    {
        cout << "Error finding images-file: Don't forget to load images folder (with all images) and image.txt file (with all indexes) to data/" << endl;
    }
    // toggle console output off
    streambuf *old = cout.rdbuf(0);

    if(computeFromImages)
    {
        while(!is.eof())
        {
            string nextImage;
            is >> nextImage;
            imageNames.push_back(nextImage);
            cout << nextImage << " referenced." << endl;
            n_img++;
        }
        is.close();
    }
    else
    {
        // use n_img from sift features file
        ifstream is(outputSiftFeatures);
        if(!is.good())
        {
            cout << "Error finding file:" << outputSiftFeatures << endl;
        }
        string num;
        is >> num;
        n_img = atof(num.c_str());
        is.close();
    }
    cout << "Total number of images: " << n_img << endl;        // set by reading images file or siftFeatures file

    // point visibility stuff
    get3DPointVisibility();                 // reads point file and fills pointsToSift (partly) and pointsInImage
    getPointsToTest();                      // function to fill pointsToTest


    // group organisation variables
    tol_angle = 0.25;
    highestGroupIdx = 0;
    countComparisons = 0;
    nextGroupIdx.push_back(highestGroupIdx);
    pointToGroup = vector<int>(n_points,-1);

    // toggle cout on
    cout.rdbuf(old);

    // Sift stuff
    get3DPointSiftRepresentations();                        // fills siftFeatureVector and pointsToSift (complete)
}

// Destructor
detectRepPoints::~detectRepPoints()
{

}

// reads point file and fills pointsToSift (partly) and pointsInImage
int detectRepPoints::get3DPointVisibility()
{

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
    pointsInImage = Eigen::MatrixXf::Zero(n_img,n_points);      // binary matrix showing which points are seen in which image

    // fill containers
    for(int i = 0; i<n_points; i++)
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
            pointsInImage(imIdx,i) = 1;

            cout << "Sift descriptor " << k << " in image " <<  pointsToSift.at(i).imIndex.back()
                 << " (" << pointsToSift.at(i).siftPos.back()(0) << ","
                 <<  pointsToSift.at(i).siftPos.back()(1) << ")" << endl;
        }
    }

    is.close();
    return 0;
}

int detectRepPoints::getPointsToTest()
{
    // binary matrix with comparisons to execute (upper triangle interesting)
    pointsToTest = Eigen::MatrixXf::Zero(n_points,n_points);

    // check column i with columns i+1 to end
    for (int i = 0; i<n_points; i++)
    {
        for (int j = i+1; j<n_points; j++)
        {
            // calculate sum of compared point columns in pointsInImage
            // extract column vector of base
            Eigen::MatrixXf sum = pointsInImage.block(0,i,n_img,1)
                    + pointsInImage.block(0,j,n_img,1);

            // comparison needed if the two points are seen together in at least one image
            if(sum.maxCoeff()==2)
                pointsToTest(i,j) = 1;
            else
                pointsToTest(i,j) = 0;
        }
    }

    // show points to be compared
    cout << "PointsToTest" << endl << pointsToTest << endl;

    return 0;
}

int detectRepPoints::get3DPointSiftRepresentations()
{
    // open file reading
    ifstream is(outputSiftFeatures);
    if(!is.good())
    {
        cout << "Problems opening " << outputSiftFeatures << endl;
        return -1;
    }
    int n_img_not_needed;
    is >> n_img_not_needed;

    // get sift features for each point
    for (int i = 0;i<n_points; i++)
    {
        cout << "Progress getting sift representations: "<< floor(100*(double)i/(double)n_points) << " %" << endl;

        // toggle console output off
        streambuf *old = cout.rdbuf(0);

        cout << "#############################################" << endl;
        cout << "Point " << i << " has sift features: " << endl;

        if(!computeFromImages)
        {
            string newPoint;  // read - symbol for checking file format
            is >> newPoint;
            cout << "new Point being read (read " << newPoint << " in " << outputSiftFeatures << endl;
        }

        // vector for currrent point to add to siftFeatureVector
        vector<Eigen::MatrixXf > FeaturesOfOnePoint;
        for (int j = 0; j<pointsToSift[i].imIndex.size(); j++) // for each sift feature
        {
            cout << j << ": ";

            Eigen::MatrixXf singleFeatureVector(siftFeatureDim,1);  // filled by function below,assuming descriptor has 128 elements

            if(computeFromImages)
            {
                // collect relevant information for current point
                int image = pointsToSift[i].imIndex.at(j);
                Eigen::Vector2f pos;
                pos << pointsToSift.at(i).siftPos.at(j)(0), pointsToSift.at(i).siftPos.at(j)(1);

                // compute sift vector
                computeSiftDescriptor(image,pos,singleFeatureVector);
            }
            else // read sift features from file
            {
                for(int h = 0; h<siftFeatureDim; h++)
                {
                    string num;
                    is >> num;
                    int siftElement = atof(num.c_str());
                    singleFeatureVector(h,0) = siftElement;
                }
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
    if(computeFromImages)
    {
        writeSiftFeaturesToFile();
        cout << "Saved point feature vectors to " << outputSiftFeatures << endl;
    }

    return 0;
}

// function to compute siftDescriptor using openCV
int detectRepPoints::computeSiftDescriptor(int imageIndex, Eigen::Vector2f pos, Eigen::MatrixXf &outSingleFeatureVector)
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

    // convert descriptor to Eigen::MatrixXf (column vector)
    Eigen::MatrixXf outDescriptor(descriptor.cols,1);
    for(int i = 0; i<descriptor.cols;i++)
    {
        outDescriptor(i,0) = descriptor.at<float>(0,i);
    }

    // pass result to specified Eigen::MatrixXf
    outSingleFeatureVector = outDescriptor;

    return 0;
}

// function to calculate angle between two descriptors
double detectRepPoints::angleOfTwoSift(Eigen::MatrixXf sift1, Eigen::MatrixXf sift2)
{
    // check dimensions
    if(sift1.cols()!=1 || sift2.cols()!=1 || sift1.rows()!=sift2.rows())
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
    return acos(dotProduct/(sift1.norm()*sift2.norm()));
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
        Eigen::MatrixXf sift1 = siftFeatureVector.at(pointIdx1).at(i);

        for(int j = 0; j < n_sift_point2; j++) // each sift descriptor of point 2
        {

            // get the descriptor
            Eigen::MatrixXf sift2 = siftFeatureVector.at(pointIdx2).at(j);

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
    for(int i=0; i<groupToPoints.at(groupIdx).size();i++)
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
            // add point 2 to group of point 1
            pointToGroup[pointIdx2] = pointToGroup[pointIdx1];                                  // use same group for point 2
            groupToPoints.at(pointToGroup[pointIdx2]).push_back(pointIdx2);
            cout << "Point " << pointIdx2 << " added to group of point " << pointIdx1 << ": ";
            coutGroupContent(pointToGroup[pointIdx1]);

        }
        else if(pointToGroup[pointIdx1]==-1 && pointToGroup[pointIdx2]!=-1)  // point 2 is assigned to a group
        {
            // add point 1 to group of point 2
            pointToGroup[pointIdx1] = pointToGroup[pointIdx2];                                  // use same group for point 2
            groupToPoints.at(pointToGroup[pointIdx1]).push_back(pointIdx1);
            cout << "Point " << pointIdx1 << " added to group of " << pointIdx2 << ": ";
            coutGroupContent(pointToGroup[pointIdx2]);
        }

        // 3) both points have been assinged -> merge the groups
        else if(pointToGroup[pointIdx1]!=-1 && pointToGroup[pointIdx2]!=-1)
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
                for(int i = 0; i<pointToGroup.size(); i++)
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

        } // end of case merging
        else
            cout << "This case should never occur." << endl;
    } // end case repetitive points

}

// main function to get repetitive points
int detectRepPoints::getRepetitivePoints()
{
    int matchResult;
    for (int i = 0; i<n_points; i++)
    {
        for (int j = i+1; j<n_points; j++)
        {
            // output progress of grouping
            cout << "Progress of grouping (more than): "<< floor(100*(double)(countComparisons)/double(n_points*n_points/2)) << " %"<< endl;

            // toggle console output off
            streambuf *old = cout.rdbuf(0);

            // compare points if needed and {not already in same group, but assiged}
            if(pointsToTest(i,j)==1 && pointToGroup[i]==pointToGroup[j] && pointToGroup[i]!=-1)
                cout << "Points " << i << " and " << j << " are already in same group." << endl;

            if(pointsToTest(i,j)==1 && (pointToGroup[i]!=pointToGroup[j] || pointToGroup[i]==-1 || pointToGroup[j]==-1))
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
    // output statistics
    cout << "Statistics: " << endl;
    cout << "Total number of points (" << n_points
         << "), exhaustive would be (" << 0.5*(n_points*n_points-n_points) << ") comparisons." << endl;
    cout << "Total number of comparisons executed: " << countComparisons << endl;

    // output groups
    for(int i = 0; i<groupToPoints.size();i++)
    {
        coutGroupContent(i);
    }

    // output point coordinates
    for(int i = 0; i<n_points; i++)
    {
        cout << "point "<< i << " has coordinates" << get3dFromPointIdx(i).transpose() << endl;
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
    // calculate group members indexes
    getRepetitivePoints();

    // build a vector where each element contains a vector of 3d points that belong to that group
    for(int i = 0; i<groupToPoints.size(); i++)
    {
        // build vector with points of current group
        vector<Eigen::Vector3d> currentGroupPoints;
        for(int j = 0; j<groupToPoints.at(i).size();j++)
        {
            currentGroupPoints.push_back(get3dFromPointIdx(groupToPoints.at(i).at(j)));
        }

        // add current group to groupsOfPoints if not empty
        if(currentGroupPoints.size()!=0)                        // indexes of groupsToPoints != group index (empty groups left out)
            groupsOfPoints.push_back(currentGroupPoints);
    }
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
    for(int i = 0; i<groupsOfPoints.size();i++)
    {
        of << "-" << endl;
        for(int j = 0; j<groupsOfPoints[i].size();j++)
        {
            of << groupsOfPoints.at(i).at(j).transpose() << endl;
        }
    }
    of.close();

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
    for(int i = 0; i<n_points;i++)
    {
        of << "-" << endl;  // new point
        for(int j = 0; j<siftFeatureVector[i].size();j++)
        {
            of << siftFeatureVector.at(i).at(j).transpose() << endl;
        }
    }
}
