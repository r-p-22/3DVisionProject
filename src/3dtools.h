#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>

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


inline void get3Dpoints(char* file, vector<Eigen::Vector3d> &points3d,
		vector<TriangulatedPoint> &pointModel){

	std::ifstream instream(file);


	if (!instream){
		cout << "file not opened" << endl;
	}
	int nPoints;
	instream >> nPoints;

	//TODO: CAREFUL: Change nPoints
	nPoints = 1000;

	for (int j = 0; j < nPoints; ++j) {
		TriangulatedPoint X;

		instream >> X.pos[0] >> X.pos[1] >> X.pos[2];
		int nMeasurements = 0;
		instream >> nMeasurements;
		for (int k = 0; k < nMeasurements; ++k) {
			PointMeasurement m;
			instream >> m.view >> m.id >> m.pos[0] >> m.pos[1];
			X.measurements.push_back(m);
		}
		pointModel.push_back(X);
		points3d.push_back(X.pos);
	}

	return ;
}

inline void getCameras(char* file, char* fileK, vector<Eigen::Matrix<double,3,4> > &cameraPoses,
		Eigen::Matrix<double,3,3> &K){

	std::ifstream is(file);
	vector<TriangulatedPoint> pointModel;
	vector<Eigen::Vector3d> points3d;

	if (!is){
		cout << "file not opened" << endl;
	}
	vector<int> viewIds;
	;

	int nViews;
	is >> nViews;
	cout << "Going to read " << nViews << " poses." << endl;

	Eigen::Matrix<double,3,4> P;

	for (int i = 0; i < nViews; ++i)
	{
		int viewId;
		is >> viewId;
		viewIds.push_back(viewId);

		is >> P(0,0) >> P(0,1) >> P(0,2) >> P(0,3);
		is >> P(1,0) >> P(1,1) >> P(1,2) >> P(1,3);
		is >> P(2,0) >> P(2,1) >> P(2,2) >> P(2,3);
		cameraPoses.push_back(P);
	}

	std::ifstream is2(fileK);
	is2 >> K(0,0) >> K(0,1) >> K(0,2) >> K(1,1) >> K(1,2);

}


inline bool compareSiftFronto(Eigen::Vector3d const &referencePoint, Eigen::Vector3d const &pointToTest,
		Eigen::Vector4d plane,
		//TriangulatedPoint Xref,
		vector<Eigen::Matrix<double,3,4>> camPoses, Eigen::Matrix3d K ){


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
		int view = i;
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

	//s1 = computeSIFT(bestview,pbest);

	cosangle = 0;
	bestview;
	pbest;
	for (int i=0; i<camPoses.size(); i++){

	//check angle between camera-point line and plane normal
		Vector3d line = pointToTest - camPoses[i].block<3,1>(0,3);
	//abs because we dont know the plane orientation
		tmpcosangle = abs(line.dot(plane.head(3)))/(line.norm()*plane.head(3).norm());

	//project point into image
		cam.setOrientation(camPoses[i]);
		p = cam.projectPoint(pointToTest);

		if ((tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
			cosangle = tmpcosangle;
			bestview = i;
			pbest = p;
		}
	}

	//s2 = computeSIFT(bestview,pbest);

	return true;

}
#endif // TOOLs_H
