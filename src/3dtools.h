#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>
#include <string>

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

		if ( (imageNames[view].compare("images/P1010480.JPG") !=0) && (imageNames[view].compare("images/P1010481.JPG") != 0)
				&& (imageNames[view].compare("images/P1010482.JPG") !=0)
				&& (imageNames[view].compare("images/P1010483.JPG") !=0) )
		{
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

	//s1 = computeSIFT(bestview,pbest);

	cosangle = 0;
	bestview;
	pbest;
	for (int i=0; i<camPoses.size(); i++){

		int view = viewIds[i];

		// QUICK AND DIRTY SOLUTION FOR LESS IMAGES
		if ( (imageNames[view].compare("images/P1010480.JPG") !=0) && (imageNames[view].compare("images/P1010481.JPG") != 0)
						&& (imageNames[view].compare("images/P1010482.JPG") !=0)
						&& (imageNames[view].compare("images/P1010483.JPG") !=0) )
		{
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

	//s2 = computeSIFT(bestview,pbest);

	return true;

}
#endif // TOOLs_H
