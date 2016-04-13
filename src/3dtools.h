#ifndef TOOLS_H
#define TOOLS_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include <fstream>

#include "CImg.h"

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


inline vector<Eigen::Vector3d> get3Dpoints(char* file){

	std::ifstream instream(file);
	vector<TriangulatedPoint> pointModel;
	vector<Eigen::Vector3d> points3d;

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

	return points3d;
}



#endif // TOOLs_H
