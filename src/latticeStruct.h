#ifndef LATTICESTRUCT_H
#define LATTICESTRUCT_H


#include <Eigen/Dense>
#include <vector>
#include <math.h>

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
	Eigen::Vector3d lowerLeftCorner;
};

inline int computeNumberOfCells(LatticeStructure latt){
	int numberOfCells = latt.width * latt.height;
	return numberOfCells;
}


//Consolidate/merge the initial lattices, to produce a final lattice list
inline std::vector<LatticeStructure> consolidateLattices(std::vector<LatticeStructure> inputLattices){

	std::vector<LatticeStructure> finalLattices;

	finalLattices.push_back(inputLattices[0]);

	bool matched;

	std::vector<LatticeStructure>::iterator initialLatticeIterator;
	for (initialLatticeIterator = inputLattices.begin()+1; initialLatticeIterator != inputLattices.end(); initialLatticeIterator++){

		LatticeStructure latt = *initialLatticeIterator;
		matched = false;

		std::vector<LatticeStructure>::iterator finalLattIterator;
		for (size_t i = 0; i < finalLattices.size(); i++){

			LatticeStructure lattF = finalLattices[i];

			//if translation vector is less that a threshold, then merge
			double basisVecThresh = sqrt((lattF.basisVectors[0] - lattF.basisVectors[1]).squaredNorm());
			double costheta = latt.plane.dot(lattF.plane)/(latt.plane.norm()*lattF.plane.norm());
			if ( (abs(latt.plane[3] - lattF.plane[3]) < basisVecThresh*0.1 ) && (acos(costheta) <= 0.0349 ) )  {
				matched = true;

				int ncells = computeNumberOfCells(latt);
				int ncellsF = computeNumberOfCells(lattF);

				//the final (merged) lattice will be the one with most cells
				if (ncells > ncellsF){
					finalLattices[i] = latt;
				}

				break;
			}

		}
		if (!matched){
			finalLattices.push_back(latt);
		}

	}
	return finalLattices;
}



#endif
