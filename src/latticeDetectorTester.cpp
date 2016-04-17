#include "latticeDetectorTester.h"
#include <iostream>



LatticeDetectorTester::LatticeDetectorTester(){
	latticeDetector = LatticeDetector();

	vector<Vector3d> validGridPoints = vector<Vector3d>();
	validGridPoints.push_back(Vector3d(0,0,0));
	validGridPoints.push_back(Vector3d(3,0,0));
	validGridPoints.push_back(Vector3d(8,0,0));
	validGridPoints.push_back(Vector3d(9,0,0));

	validGridPoints.push_back(Vector3d(0,1,0));
	validGridPoints.push_back(Vector3d(5,1,0));
	validGridPoints.push_back(Vector3d(6,1,0));
	validGridPoints.push_back(Vector3d(7,1,0));
	validGridPoints.push_back(Vector3d(8,1,0));
	validGridPoints.push_back(Vector3d(9,1,0));

	validGridPoints.push_back(Vector3d(0,2,0));
	validGridPoints.push_back(Vector3d(1,2,0));
	validGridPoints.push_back(Vector3d(3,2,0));
	validGridPoints.push_back(Vector3d(4,2,0));
	validGridPoints.push_back(Vector3d(6,2,0));
	validGridPoints.push_back(Vector3d(7,2,0));
	validGridPoints.push_back(Vector3d(8,2,0));

	validGridPoints.push_back(Vector3d(2,3,0));
	validGridPoints.push_back(Vector3d(5,3,0));
	validGridPoints.push_back(Vector3d(7,3,0));
	validGridPoints.push_back(Vector3d(9,3,0));

	validGridPoints.push_back(Vector3d(0,4,0));
	validGridPoints.push_back(Vector3d(2,4,0));
	validGridPoints.push_back(Vector3d(5,4,0));
	validGridPoints.push_back(Vector3d(6,4,0));
	validGridPoints.push_back(Vector3d(9,4,0));

	validGridPoints.push_back(Vector3d(3,5,0));
	validGridPoints.push_back(Vector3d(4,5,0));
	validGridPoints.push_back(Vector3d(5,5,0));
	validGridPoints.push_back(Vector3d(6,5,0));
	validGridPoints.push_back(Vector3d(7,5,0));

	validGridPoints.push_back(Vector3d(0,6,0));
	validGridPoints.push_back(Vector3d(1,6,0));
	validGridPoints.push_back(Vector3d(2,6,0));
	validGridPoints.push_back(Vector3d(8,6,0));
	validGridPoints.push_back(Vector3d(9,6,0));

	validGridPoints.push_back(Vector3d(0,7,0));
	validGridPoints.push_back(Vector3d(2,7,0));
	validGridPoints.push_back(Vector3d(3,7,0));
	validGridPoints.push_back(Vector3d(5,7,0));
	validGridPoints.push_back(Vector3d(6,7,0));
	validGridPoints.push_back(Vector3d(7,7,0));

	validGridPoints.push_back(Vector3d(3,8,0));
	validGridPoints.push_back(Vector3d(5,8,0));
	validGridPoints.push_back(Vector3d(8,8,0));
	validGridPoints.push_back(Vector3d(9,8,0));

	validGridPoints.push_back(Vector3d(1,9,0));
	validGridPoints.push_back(Vector3d(2,9,0));
	validGridPoints.push_back(Vector3d(6,9,0));
	validGridPoints.push_back(Vector3d(7,9,0));
	validGridPoints.push_back(Vector3d(8,9,0));

	latticeDetector.validGridPoints = validGridPoints;

	vector<Vector3d> reconstructedPoints = vector<Vector3d>();

	reconstructedPoints.push_back(Vector3d(3,0,0));

	reconstructedPoints.push_back(Vector3d(7,1,0));

	reconstructedPoints.push_back(Vector3d(1,2,0));
	reconstructedPoints.push_back(Vector3d(7,2,0));
	reconstructedPoints.push_back(Vector3d(8,2,0));

	reconstructedPoints.push_back(Vector3d(7,3,0));

	reconstructedPoints.push_back(Vector3d(2,4,0));
	reconstructedPoints.push_back(Vector3d(5,4,0));

	reconstructedPoints.push_back(Vector3d(6,5,0));
	reconstructedPoints.push_back(Vector3d(7,5,0));

	reconstructedPoints.push_back(Vector3d(8,6,0));

	reconstructedPoints.push_back(Vector3d(2,7,0));
	reconstructedPoints.push_back(Vector3d(3,7,0));
	reconstructedPoints.push_back(Vector3d(5,7,0));

	reconstructedPoints.push_back(Vector3d(5,8,0));
	reconstructedPoints.push_back(Vector3d(8,8,0));

	reconstructedPoints.push_back(Vector3d(7,9,0));
	reconstructedPoints.push_back(Vector3d(8,9,0));

	latticeDetector.reconstructedPoints = reconstructedPoints;

}

LatticeDetectorTester::~LatticeDetectorTester(){

}

void LatticeDetectorTester::testCalculateCandidateVectors(){

	//vector<Vector3d> candidates = latticeDetector.calculateCandidateVectors(true); //works
	//vector<Vector3d> candidates = latticeDetector.calculateCandidateVectors(false); //works
	// Just add another line to have big enough scope to debug properly
	bool a = true;

}

void LatticeDetectorTester::testClusterCandidates(){

	vector<Vector3d> candidates = vector<Vector3d>();
	candidates.push_back(Vector3d(0,0,1));
	candidates.push_back(Vector3d(0,0,-1));

	list<list<Vector3d> > clusteredCandidates = latticeDetector.clusterCandidates(candidates);

	// Just add another line to have big enough scope to debug properly
		bool a = true;
}

void LatticeDetectorTester::testCombineCandidates(){

	list<list<Vector3d>> clusteredCandidates = list<list<Vector3d> >(0);

	list<Vector3d> candidate1 = list<Vector3d>(0);
	candidate1.push_back(Vector3d(0,0,1));
	candidate1.push_back(Vector3d(0,0,-1));

	clusteredCandidates.push_back(candidate1);

	vector<Vector3d> finalCandidates = vector<Vector3d>();
	vector<int> scores = vector<int>();

	latticeDetector.combineCandidates(clusteredCandidates,finalCandidates,scores);

	// Just add another line to have big enough scope to debug properly
			bool a = true;

}

void LatticeDetectorTester::testTranslationVector(){

}

void LatticeDetectorTester::testVectorsAreSimilar(){

}

void LatticeDetectorTester::testValidateCandidateVectors(){

}

void LatticeDetectorTester::testValidateVector(){

}

void LatticeDetectorTester::testValidInvalidRatio(){

}

void LatticeDetectorTester::testIsPointValid(){

}

void LatticeDetectorTester::testProjectPointsOnLine(){

}

void LatticeDetectorTester::testGetOutermostOnGridPointIndices(){

}

void LatticeDetectorTester::testCalculateLatticeBoundary(){

}

void LatticeDetectorTester::testValidLine(){

}

void LatticeDetectorTester::test(){

	vector<Vector3d> candidates = latticeDetector.calculateCandidateVectors(false);
	//vector<double> scores = latticeDetector.validateCandidateVectors(candidates);
	vector<Vector3d> finalBasisVectors = latticeDetector.getFinalBasisVectors(candidates);

	//vector<Vector3d> latticeBoundaries = latticeDetector.calculateLatticeBoundary(latticeDetector.reconstructedPoints[8],finalBasisVectors[0], finalBasisVectors[1]);

	testCalculateCandidateVectors();

	testClusterCandidates();

	testCombineCandidates();

	testTranslationVector();

	testVectorsAreSimilar();

	testValidateCandidateVectors();

	testValidateVector();

	testValidInvalidRatio();

	testIsPointValid();

	testProjectPointsOnLine();

	testGetOutermostOnGridPointIndices();

	testCalculateLatticeBoundary();

	testValidLine();

}
