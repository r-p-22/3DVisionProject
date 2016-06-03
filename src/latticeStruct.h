#ifndef LATTICESTRUCT_H
#define LATTICESTRUCT_H


#include <Eigen/Dense>
#include <vector>
#include <math.h>

using namespace std;

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
	  Eigen::Vector3d               pos;
      std::vector<PointMeasurement> measurements;

      TriangulatedPoint() : measurements()
      {
          pos=Eigen::Vector3d::Zero();
      }

      TriangulatedPoint(Eigen::Vector3d const& p, std::vector<PointMeasurement> const& ms)
         : pos(p), measurements(ms)
      { }

}; // end struct TriangulatedPoint



//class inputManager;


struct LatticeStructure
{
	Eigen::Vector4d plane;
	std::vector<Eigen::Vector3d> basisVectors;
	int width; // in the direction of basisVectors[0]
	int height; // in the direction of basisVectors[1]
	Eigen::Vector3d corner;
};

inline int computeNumberOfCells(LatticeStructure latt){
	int numberOfCells = latt.width * latt.height;
	return numberOfCells;
}



#endif
