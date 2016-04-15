#include "latticeDetectorTester.h"



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

	vector<Vector3d> candidates = latticeDetector.calculateCandidateVectors(true);

}

void LatticeDetectorTester::testClusterCandidates(){

}

void LatticeDetectorTester::testCombineCandidates(){

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
