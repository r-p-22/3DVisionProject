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
#include "inputManager.h"
//#include "latticeStruct.h"

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

class inputManager;
inline void projectLattice(inputManager inpM, LatticeStructure latt){

//	Vector3d LL = latt.boundary[0];
//	Vector3d TR = latt.boundary[1];
	Vector3d basis1 = latt.basisVectors[0];
	Vector3d basis2 = latt.basisVectors[1];

	/*Matrix<double,3,3> A;
		Vector3d solution;
		A.block<3,1>(0,2) = -(TR-LL);
		A.block<3,1>(0,0) = basis1;
		A.block<3,1>(0,1) = basis2;
		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeFullV);
		solution = svd.matrixV().block<3,1>(0,2);//<sizeRows,sizeCols>(beginRow,beginCol)
		solution = solution/solution[2];
*/
		int k1 = latt.width;
		int k2 = latt.height;

		Vector3d B1 = latt.lowerLeftCorner + k1*basis1;
		Vector3d B2 = latt.lowerLeftCorner + k2*basis2;

		string img = inpM.getImgNames()[45];
		int i=0;
		for (i=0;i<inpM.getCamPoses().size();i++){
			if (inpM.getViewIds()[i] == 45)
				break;
		}

		Eigen::Matrix<double,3,4> P = inpM.getCamPoses()[i];
		float const w = 1696;
		float const h = 1132;
		CImg<unsigned char> image(("data/"+img).c_str());
		const unsigned char color[] = { 0,0,255 };

		CameraMatrix cam;
		cam.setIntrinsic(inpM.getK());

		cam.setOrientation(P);
		Vector3d pa; pa = latt.lowerLeftCorner;
		Vector3d pb; pb = B1;
		Vector2d pa2d, pb2d;
		for (int k=0; k<=k2; k++){
			//project pa, pb into image
       		pa2d = cam.projectPoint(pa);
       		pb2d = cam.projectPoint(pb);
			//plot 2D line into image
	        image.draw_line(pa2d[0],pa2d[1],pb2d[0],pb2d[1],color);

			pa += basis2;
			pb += basis2;

		}
		pa = latt.lowerLeftCorner;
	    pb = B2;
		for (int k=0; k<=k1; k++){
			//project pa, pb into image
			pa2d = cam.projectPoint(pa);
			pb2d = cam.projectPoint(pb);
			//plot 2D line into image
	        image.draw_line(pa2d[0],pa2d[1],pb2d[0],pb2d[1],color);

			pa += basis1;
			pb += basis1;
		}

		CImgDisplay main_disp(image,"");

		image.save("latt_view45.png");
		while (!main_disp.is_closed()){
		    main_disp.wait();
		}
}

inline void projectGroup(inputManager inpM,vector<Vector3d> group){


	string img = inpM.getImgNames()[45];
	int i=0;
	for (i=0;i<inpM.getCamPoses().size();i++){
		if (inpM.getViewIds()[i] == 45)
			break;
	}

	Eigen::Matrix<double,3,4> P = inpM.getCamPoses()[i];
	float const w = 1696;
	float const h = 1132;
	CImg<unsigned char> image(("data/"+img).c_str());
	const unsigned char color[] = { 0,0,255 };

	CameraMatrix cam;
	cam.setIntrinsic(inpM.getK());

	cam.setOrientation(P);

	Vector2d pa2d;
	for (int k=0; k < group.size(); k++){
		pa2d = cam.projectPoint(group[k]);
		image.draw_circle(pa2d[0],pa2d[1],5,color,1);
	}

	CImgDisplay main_disp(image,"");

	image.save("latt_view45.png");
	while (!main_disp.is_closed()){
		main_disp.wait();
	}


}




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

   /* cv::Mat roi(mask, cv::Rect(rectTLx,rectTLy,rectBRx-rectTLx,rectBRy-rectTLy));
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
    }*/

    int KPsize = 9;

    // compute descriptor of desired keypoint location using calculated size
    cv::SIFT siftDetector;
    cv::KeyPoint myKeypoint(x, y, KPsize);
    vector<cv::KeyPoint> myKeyPoints;
    myKeyPoints.push_back(myKeypoint);
    cv::Mat descriptor;

    siftDetector(input,mask,myKeyPoints,descriptor, true);

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
		Eigen::Matrix3d K, vector<Eigen::Matrix<double,3,4>> camPoses, vector<int> viewIds,  vector<string> imageNames ){

	//TODO: We use fixed image size (1696x1132). Must read img to check real...
	float const w = 1696;
	float const h = 1132;

	CameraMatrix cam;
	cam.setIntrinsic(K);

	double cosangle = 0;
	double tmpcosangle=0;
	int bestview=-1;
	Vector2d pbest;
	Vector2d p;

	for (int i=0; i<camPoses.size(); i++){

	//get view
		int view = viewIds[i];

		if ((view < 45) || (view > 47))
		{
			continue;
		}

	//check angle between camera-point line and plane normal
		Vector3d line = referencePoint - camPoses[i].block<3,1>(0,3);
	//abs because we dont know the plane orientation
		tmpcosangle = abs(line.dot(plane.head(3)))/sqrt(line.squaredNorm() * plane.head(3).squaredNorm());

	//project point into image
		cam.setOrientation(camPoses[i]);
		p = cam.projectPoint(referencePoint);

	//TODO: Handle the case where camera is NOT facing the point.
		double d = cam.transformPointIntoCameraSpace(referencePoint)[2];

		if ((d > 0) && (tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
			cosangle = tmpcosangle;
			bestview = view;
			pbest = p;
		}
	}


	//cout << pbest << endl;
	Eigen::VectorXd s1;
	computeSIFT(imageNames[bestview],pbest,s1);

	//==========================

	cosangle = 0;
	pbest[0]=-1;
	pbest[1]=-1;
	bestview = -1;

	for (int i=0; i<camPoses.size(); i++){

		int view = viewIds[i];

		if ((view < 45) || (view > 47))
			{
			continue;
		}

	//check angle between camera-point line and plane normal
		Vector3d line = pointToTest - camPoses[i].block<3,1>(0,3);
	//abs because we dont know the plane orientation
		tmpcosangle = abs(line.dot(plane.head(3)))/(line.norm()*plane.head(3).norm());

	//project point into image
		cam.setOrientation(camPoses[i]);
		p = cam.projectPoint(pointToTest);
		double d = cam.transformPointIntoCameraSpace(pointToTest)[2];
		//cout << p << endl;
		if ((d > 0)&&(tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
			cosangle = tmpcosangle;
			bestview = view;
			pbest = p;

		}
	}

	if (bestview == -1){
		return false;
	}
	//cout << "test bestview: " << bestview << endl;
	//cout << pbest << endl;
	//cout << "==" << endl;

	Eigen::VectorXd s2;
	computeSIFT(imageNames[bestview],pbest,s2);

	double dotProduct = s1.dot(s2);

	/*
	   double dotProduct = 0.0;
	   for (int i = 0; i<s1.rows();i++)
	{
		dotProduct += s1(i,0)*s2(i,0);
	}
	*/

	 //angle (in radians)
	 double theta = acos(dotProduct/(s1.norm()*s2.norm()));


	 if ((theta) < 2*0.25) // 2*tol_angle from detectRepPoints
	 	 return true;
	 else
		 return false;

}
#endif // TOOLs_H
