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

#include "inputManager.h"

using namespace Eigen;
using namespace std;

class LatticeDetector{

public:
	static constexpr double VECTOR_DISTANCE = 0.06;

	static constexpr double TRESHOLD1 = 0.1;

	static constexpr double TRESHOLD2 = 0.5;

	vector<Vector3d> reconstructedPoints;

	vector<Vector3d> calculateCandidateVectors(bool naive);

	list<list<Vector3d> > clusterCandidates(vector<Vector3d> const &candidates);

	void combineCandidates(list<list<Vector3d> > const &clusteredCandidates, vector<Vector3d> &finalCandidates, vector<int> &scores);

	Vector3d translationVector(Vector3d const &point1, Vector3d const &point2);

	bool vectorsAreSimilar(Vector3d const &vector1, Vector3d const &vector2);

	vector<double> validateCandidateVectors(vector<Vector3d> const &candidateVectors);

	double validateVector(Vector3d const &candidateVector);

	double validInvalidRatio(Vector3d const &referencePoint, Vector3d const &candidateVector);

	bool isPointValid(Vector3d const &referencePoint, Vector3d const &pointToTest, Vector3d const &vector);

	vector<Vector3d> projectPointsOnLine(Vector3d const &referencePoint, Vector3d const &candidateVector);

	vector<int> getOutermostOnGridPointIndices(vector<Vector3d> const &projectedPoints, Vector3d const &referencePoint, Vector3d const &candidateVector);

	void calculateLatticeBoundary(Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &lowerLeftCornerOut, int &widthOut, int &heightOut);

	void latticeBoundaryForReferencePoint(Vector3d const &referencePoint, Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &lowerLeftCornerOut, int &withOut, int &heightOut);

	bool validLine(Vector3d const &referencePoint, Vector3d const &anchorPoint, Vector3d const &directionVector, int length);

	bool pointEqualsGridPoint(Vector3d point, Vector3d gridPoint, Vector3d vector);

	vector<pair<int,vector<int> > > getOnGridIndices(vector<int> inputIndices, LatticeStructure lattice);

	vector<Vector3d> changeToLatticeBasis(vector<Vector3d> const &points, Vector3d const &latticeVector1, Vector3d const &latticeVector2);



//=================================================
//Nektarios's

private:


	bool isIntegerCombination(int i,vector<Vector3d>& candidatesInOrder,vector<bool>& valid);

public:

	inputManager* inpManager;

	vector<Vector3d> getFinalBasisVectors(vector<Vector3d>& candidateVectors);

	Eigen::Vector4d plane;

};



#endif /* SRC_LATTICEDETECTOR_H_ */
