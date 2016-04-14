/*
 * latticeDetector.h
 *
 *  Created on: Apr 4, 2016
 *      Author: valaina
 */

#ifndef SRC_LATTICEDETECTOR_H_
#define SRC_LATTICEDETECTOR_H_

//#include <3dtools.h>
#include <list>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

class LatticeDetector{

public:
	static constexpr double VECTOR_DISTANCE = 0;

	static constexpr double TRESHOLD1 = 0.1;

	static constexpr double TRESHOLD2 = 0.5;

	vector<Vector3d> points;

	vector<Vector3d> calculateCandidateVectors(bool naive);

	void clusterCandidates(vector<Vector3d> const &candidates, list<list<Vector3d> > &clusteredCandidates);

	void combineCandidates(list<list<Vector3d> > const &clusteredCandidates, vector<Vector3d> &finalCandidates, vector<int> &scores);

	Vector3d translationVector(Vector3d const &point1, Vector3d const &point2);

	bool vectorsAreSimilar(Vector3d const &vector1, Vector3d const &vector2);

	vector<double> validateCandidateVectors(vector<Vector3d> const &candidateVectors);

	double validateVector(Vector3d const &candidateVector);

	double validInvalidRatio(Vector3d const &referencePoint, Vector3d const &candidateVector);

	bool isPointValid(Vector3d const &referencePoint, Vector3d const &pointToTest);

	vector<Vector3d> projectPointsOnLine(Vector3d const &referencePoint, Vector3d const &candidateVector);

	vector<int> getOutermostOnGridPointIndices(vector<Vector3d> const &projectedPoints, Vector3d const &referencePoint, Vector3d const &candidateVector);

	vector<Vector3d> calculateLatticeBoundary(Vector3d const &referencePoint, Vector3d const &latticeVector1, Vector3d const &latticeVector2);

	bool validLine(Vector3d const &referencePoint, Vector3d const &anchorPoint, Vector3d const &directionVector, int length);

//=================================================
//Nektarios's

	Eigen::Vector4d plane;

	Eigen::Matrix3d K;
	vector< Eigen::Matrix<double,3,4> > camPoses;
	vector<Vector3d> getFinalBasisVectors(vector<Vector3d> candidateVectors);
	
	bool isIntegerCombination(int i,vector<Vector3d> candidatesInOrder,vector<bool> valid);

	void setPlane(Eigen::Vector4d plane){
		this->plane = plane;
	}

};



#endif /* SRC_LATTICEDETECTOR_H_ */
