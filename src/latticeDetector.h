/*
 * latticeDetector.h
 *
 *  Created on: Apr 4, 2016
 *      Author: valaina
 */

#ifndef SRC_LATTICEDETECTOR_H_
#define SRC_LATTICEDETECTOR_H_

#include "3dtools.h"
#include <list>

using namespace Eigen;
using namespace std;

class LatticeDetector{

public:
	static const double VECTOR_DISTANCE=0;

	vector<Vector3d> points;

	vector<Vector3d> calculateCandidateVectors(bool naive);

	void clusterCandidates(vector<Vector3d> &candidates, list<list<Vector3d> > &clusteredCandidates);

	void combineCandidates(list<list<Vector3d> > &clusteredCandidates, vector<Vector3d> &finalCandidates, vector<int> &scores);

	Vector3d translationVector(Vector3d &point1, Vector3d &point2);

	bool vectorsAreSimilar(Vector3d vector1, Vector3d vector2);

};



#endif /* SRC_LATTICEDETECTOR_H_ */
