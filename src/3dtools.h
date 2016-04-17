#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>

#include <opencv2/core/core.hpp>
#include <opencv2/nonfree/features2d.hpp>
#include <opencv2/highgui/highgui.hpp>

#include "CImg.h"
#include "camera.h"


using namespace std;
using namespace cimg_library;

class snake
{
private:
    vector< Eigen::Vector2f > points;
    int size;
    float gamma;
    float alpha;
    float beta;
    bool end;
    CImgList<float> extGrad;
    CImg<unsigned char> image;
    void preProcess();
public:
    snake()
    {
        size = 0;
        alpha =0; gamma=0; beta=0;
        end = false;
    }
    void assign(vector<vector<int> > ps, CImg<unsigned char> &img);
    void setParams(float _gamma, float _alpha, float _beta);
    void update();
    void display(CImg<unsigned char> &img, CImgDisplay &disp, const unsigned char *const colorPts, const unsigned char *const colorLines);
    bool stop();
};

struct PointMeasurement //position(2D) of point in the {view} image
{
	Eigen::Vector2f pos;
      int           view, id;

      PointMeasurement()
         : id(-1)
      { }

      PointMeasurement(Eigen::Vector2f const& pos_, int view_, int id_ = -1)
         : pos(pos_), view(view_), id(id_)
      { }

      bool operator==(PointMeasurement const& rhs) const
      {
         // Note: if view number and feature id are the same, we assume that the
         // 2D position is the same.
         return (this->id == rhs.id) && (this->view == rhs.view);
      }


}; // end struct PointMeasurement

struct TriangulatedPoint
{
	Eigen::Vector3d                 pos;
      std::vector<PointMeasurement> measurements;

      TriangulatedPoint() : measurements()
      {
          pos=Eigen::Vector3d::Zero();
      }

      TriangulatedPoint(Eigen::Vector3d const& p, std::vector<PointMeasurement> const& ms)
         : pos(p), measurements(ms)
      { }

}; // end struct TriangulatedPoint



template<typename T1> T1 median(vector<T1> &v)
{
    size_t n = v.size() / 2;
    nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}

inline int computeSIFT(string imagename, Eigen::Vector2d pos, Eigen::VectorXd &outSingleFeatureVector)
{
    cv::Mat B;
    float x = pos(0);
    float y = pos(1);

    const cv::Mat input = cv::imread("data/"+imagename, 0); //Load as grayscale

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
    Eigen::VectorXd outDescriptor(descriptor.cols,1);
    for(int i = 0; i<descriptor.cols;i++)
    {
        outDescriptor(i) = descriptor.at<float>(0,i);
    }

    // pass result to specified Eigen::VectorXf
    outSingleFeatureVector = outDescriptor;

    return 0;
}


inline bool compareSiftFronto(Eigen::Vector3d const &referencePoint, Eigen::Vector3d const &pointToTest,
		Eigen::Vector4d plane,
		//TriangulatedPoint Xref,
		Eigen::Matrix3d K, vector<Eigen::Matrix<double,3,4>> camPoses, vector<int> viewIds,  vector<string> imageNames ){

	//TODO: We use fixed image size (1696x1132). Must read img to check real...
	float const w = 1696;
	float const h = 1132;

	CameraMatrix cam;
	cam.setIntrinsic(K);

	double cosangle = 0;
	double tmpcosangle;
	int bestview;
	Vector2d pbest;
	Vector2d p;

	//for (int k=0; i<Xref.measurements.size(); i++){
	for (int i=0; i<camPoses.size(); i++){

	//get view

		//int view = Xref.measurements[k].view;
		int view = viewIds[i];

		// QUICK AND DIRTY SOLUTION FOR LESS IMAGES

		//if ( (imageNames[view].compare("images/P1010480.JPG") !=0) && (imageNames[view].compare("images/P1010481.JPG") != 0)
		//		&& (imageNames[view].compare("images/P1010482.JPG") !=0)
		//		&& (imageNames[view].compare("images/P1010483.JPG") !=0) )
		if (view < 45 && view > 47){
			continue;
		}

	//check angle between camera-point line and plane normal
		Vector3d line = referencePoint - camPoses[view].block<3,1>(0,3);
	//abs because we dont know the plane orientation
		tmpcosangle = abs(line.dot(plane.head(3)))/(line.norm()*plane.head(3).norm());

		cam.setOrientation(camPoses[i]);
		p = cam.projectPoint(pointToTest);


	//TODO: Handle the case where camera is NOT facing the point.

	//project point into image
		if ((tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
			cosangle = tmpcosangle;
			bestview = i;
			pbest = p;
		}
	}

	Eigen::VectorXd s1;
	computeSIFT(imageNames[bestview],pbest,s1);

	cosangle = 0;
	bestview;
	pbest;
	for (int i=0; i<camPoses.size(); i++){

		int view = viewIds[i];

		// QUICK AND DIRTY SOLUTION FOR LESS IMAGES
	//	if ( (imageNames[view].compare("images/P1010480.JPG") !=0) && (imageNames[view].compare("images/P1010481.JPG") != 0)
	//					&& (imageNames[view].compare("images/P1010482.JPG") !=0)
	//					&& (imageNames[view].compare("images/P1010483.JPG") !=0) )
			if (view < 45 && view > 47){
			continue;
		}

	//check angle between camera-point line and plane normal
		Vector3d line = pointToTest - camPoses[view].block<3,1>(0,3);
	//abs because we dont know the plane orientation
		tmpcosangle = abs(line.dot(plane.head(3)))/(line.norm()*plane.head(3).norm());

	//project point into image
		cam.setOrientation(camPoses[view]);
		p = cam.projectPoint(pointToTest);

		if ((tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
			cosangle = tmpcosangle;
			bestview = view;
			pbest = p;
		}
	}

	Eigen::VectorXd s2;
	computeSIFT(imageNames[bestview],pbest,s2);

	double dotProduct = 0.0;
	for (int i = 0; i<s1.rows();i++)
	{
		dotProduct += s1(i,0)*s2(i,0);
	}

	 //angle (in radians)
	 double theta = acos(dotProduct/(s1.norm()*s2.norm()));


	 if (theta < 2*0.25) // 2*tol_angle from detectRepPoints
	 	 return true;
	 else
		 return false;

}
#endif // TOOLs_H
