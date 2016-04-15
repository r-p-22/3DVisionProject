#ifndef SRC_LATTICEDETECTORTESTER_H_
#define SRC_LATTICEDETECTORTESTER_H_

#include "latticeDetector.h"

class LatticeDetectorTester{

public:

	LatticeDetector latticeDetector;

	LatticeDetectorTester();
	~LatticeDetectorTester();

	void testCalculateCandidateVectors();

	void testClusterCandidates();

	void testCombineCandidates();

	void testTranslationVector();

	void testVectorsAreSimilar();

	void testValidateCandidateVectors();

	void testValidateVector();

	void testValidInvalidRatio();

	void testIsPointValid();

	void testProjectPointsOnLine();

	void testGetOutermostOnGridPointIndices();

	void testCalculateLatticeBoundary();

	void testValidLine();

	void test();

};

#endif /* SRC_LATTICEDETECTORTESTER_H_ */
