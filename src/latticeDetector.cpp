#include "latticeDetector.h"

//libs for Nektarios's part
#include <algorithm>
#include <numeric>
#include <math.h>
#include "3dtools.h"

/*
 *
 *
 * # x Find translation between all point pairs (translationVector)
 *
 * # x "Cluster" translation vectors: treshold in difference, weighted average?
 *
 * # x take 2 highest peaks naively
 *
 * # x or take all peaks as candidates
 *
 * # verify basis vector
 *
 *  	x -select a 3D point as reference
 *  	x - predict grid points on line through reference point along vector
 *  	- compare SIFT descriptors of grid points to ref. point:
 *  		if alpha > 2phi -> valid
 *  			# SIFT descriptor is computed in most frontoparallel image
 *  				line connecting camera center and point is closest to plane normal
 *  	x - check every reference point like that
 *  	x - image validation score: SUM of ratios over all reference points
 *
 *  	x PER REFERENCE POINT: number of valid grid points / all grid points
 *
 *  	x - search for two farthest reconstructed on-grid points on both sides: on-grid if
 *  		dist to a certain grid point is smaller than T1=10% of basis vector length
 *  	x - move them away until T2=50% of grid points within are invalid
 *  	x - trim invalid grid points on both ends
 *
 *  # choose two basis vectors:
 *
 *  	- sort in descending order of lengths
 *  	- discard vector if integer combination of rest of the queue
 *  	- vectors along same direction: keep the one with higher score
 *  	- take two highest-score vectors of remaining ones
 */



Vector3d LatticeDetector::translationVector(Vector3d const &point1, Vector3d const &point2){

	Vector3d translation = point1-point2;
	return translation;

}

vector<Vector3d> LatticeDetector::calculateCandidateVectors(bool naive){

	std::vector<Vector3d> rawCandidates = vector<Vector3d>();

	int length = reconstructedPoints.size();

	for (int i=0; i<length; i++){
		for (int j=i+1; j<length; j++){
			Vector3d translation = translationVector(reconstructedPoints[i], reconstructedPoints[j]);
			rawCandidates.push_back(translation);
		}
	}

	// Init data structures
	vector<Vector3d> finalCandidates = vector<Vector3d>();
	vector<int> scores = vector<int>();

	// Cluster candidates that are similar
	list<list<Vector3d> > clusteredCandidates = clusterCandidates(rawCandidates);

	// Average candidates in each cluster to a final candidate, score the clusters based on the number of members
	combineCandidates(clusteredCandidates, finalCandidates, scores);

	// return the two best candidates for naive solution
	if(naive && finalCandidates.size()>2){
		vector<Vector3d> finalNaiveCandidates = vector<Vector3d>();

		std::vector<int>::iterator scoresIt;

		int bestScore = 0;
		int secondBestScore = 0;
		int bestIndex = 0;
		int secondBestIndex = 0;

		for(scoresIt = scores.begin(); scoresIt != scores.end(); ++scoresIt){
			int score = (*scoresIt);

			if(score>=bestScore){

				secondBestScore = bestScore;
				secondBestIndex = bestIndex;

				bestScore = score;
				bestIndex = std::distance(scores.begin(), scoresIt);
			}
			else if(score>=secondBestScore){
				secondBestScore = score;
				secondBestIndex = std::distance(scores.begin(), scoresIt);
			}
		}

		finalNaiveCandidates.push_back(finalCandidates[bestIndex]);
		finalNaiveCandidates.push_back(finalCandidates[secondBestIndex]);

		return finalNaiveCandidates;
	}

	else{
		return finalCandidates;
	}
}

void LatticeDetector::combineCandidates(list<list<Vector3d> > const &clusteredCandidates, vector<Vector3d> &finalCandidates, vector<int> &scores){

	std::list<list<Vector3d> >::const_iterator clusterItOuter;
	std::list<Vector3d>::const_iterator clusterItInner;

	for (clusterItOuter = clusteredCandidates.begin(); clusterItOuter != clusteredCandidates.end(); ++clusterItOuter){

		int score = 0;
		Vector3d finalCandidateSum = Vector3d(0,0,0);


		for(clusterItInner = (*clusterItOuter).begin(); clusterItInner != (*clusterItOuter).end(); ++clusterItInner){
			Vector3d nextCandidate = *clusterItInner;

			// If not the first candidate, align the candidate with the orientation of the candidates processed until now
			if (clusterItInner!=(*clusterItOuter).begin()){
				Vector3d finalCandidatePreview = finalCandidateSum/score;

				//if pointing into opposite direction, revert the candidate
				if((nextCandidate+finalCandidatePreview).norm() < (nextCandidate-finalCandidatePreview).norm()){
					nextCandidate = nextCandidate*(-1);
				}
			}

			// Sum the candidates in a cluster
			finalCandidateSum = finalCandidateSum + nextCandidate;
			// Count the candidates of the cluster
			score++;
		}

		// Average the candidates in a cluster
		Vector3d finalCandidate = finalCandidateSum/score;

		scores.push_back(score);
		finalCandidates.push_back(finalCandidate);
	}
}

bool LatticeDetector::vectorsAreSimilar(Vector3d const &vector1, Vector3d const &vector2){

	// Note: Difference of two vectors is the intuitive vector to check. However, if the vectors
	// are similar but point in opposite directions, the difference will be huge only for the orientation
	// reasons. So check as well the sum, which is the difference of the two vectors with one being turned 180°
	// by mirroring it in the origin (vector1-(vector2*(-1))

	if(((vector1-vector2).norm()<=VECTOR_DISTANCE)){
		return true;
	}
	else if(((vector1+vector2).norm()<=VECTOR_DISTANCE)){
		return true;
	}
	else{
		return false;
	}
}

// Clusters candidates together if they are similar. Clusters are lists of candidate vectors, themselves stored in a list
list<list<Vector3d> > LatticeDetector::clusterCandidates(vector<Vector3d> const &candidates){

	std::vector<Vector3d>::const_iterator candidateIt;
	std::list<list<Vector3d> >::iterator clusterItOuter;
	std::list<Vector3d>::iterator clusterItInner;

	list<list<Vector3d> > clusteredCandidates = list<list<Vector3d> >(0);

	for(candidateIt = candidates.begin(); candidateIt!=candidates.end(); ++candidateIt){

		// Make a new cluster for the candidate
		list<Vector3d> cluster = list<Vector3d>();
		cluster.push_back(*candidateIt);

		for (clusterItOuter = clusteredCandidates.begin(); clusterItOuter!=clusteredCandidates.end();){ // iterator is increased manually!

			bool inCluster = false;

			// Search all the vectors in the current cluster. If one of them is similar to the candidate vector,
			// the current cluster shall be merged into the candidate's cluster.
			// Note: More than one existing cluster can be merged into the newly formed cluster.
			for(clusterItInner = (*clusterItOuter).begin(); clusterItInner != (*clusterItOuter).end(); ++clusterItInner){
				if (vectorsAreSimilar(*candidateIt, *clusterItInner)){
					inCluster = true;
					break;
				}
			}
			if (inCluster){
				// Merge the old cluster into the new cluster
				cluster.splice(cluster.end(), *clusterItOuter);
				// Remove the old cluster from the list. Advances iterator automatically.
				clusterItOuter=clusteredCandidates.erase(clusterItOuter);
			}
			else{
				// increase iterator
				++clusterItOuter;
			}
		}

		// append the new cluster to the list
		clusteredCandidates.push_back(cluster);
	}

	return clusteredCandidates;

}

vector<double> LatticeDetector::validateCandidateVectors(vector<Vector3d> const &candidateVectors){

	vector<double> scores = vector<double>();

	std::vector<Vector3d>::const_iterator candidateIt;

	// store the score of every candidate vector
	for(candidateIt = candidateVectors.begin(); candidateIt != candidateVectors.end(); ++candidateIt){
		double score = validateVector(*candidateIt);
		scores.push_back(score);
	}

	return scores;

}

double LatticeDetector::validateVector(Vector3d const &candidateVector){

	std::vector<Vector3d>::iterator pointsIt;

	double imageValidationScore = 0;

	// sum up the score (ration between valid and invalid grid points) of every reference point
	for(pointsIt = reconstructedPoints.begin(); pointsIt != reconstructedPoints.end(); ++pointsIt){
		double pointScore = validInvalidRatio((*pointsIt), candidateVector);
		imageValidationScore = imageValidationScore + pointScore;
	}

	return imageValidationScore;
}

double LatticeDetector::validInvalidRatio(Vector3d const &referencePoint, Vector3d const &candidateVector){

	vector<Vector3d> projectedPoints = projectPointsOnLine(referencePoint, candidateVector);

	vector<int> outermostOnGridPointIndices = getOutermostOnGridPointIndices(projectedPoints, referencePoint, candidateVector);

	int minIndex = outermostOnGridPointIndices[0];
	int maxIndex = outermostOnGridPointIndices[1];

	int smallestValidIndex = minIndex;
	int highestValidIndex = maxIndex;

	int validCount = 0;
	int totalCount = 0;

	// check whether points between the outermost on grid points are valid
	for (int index = minIndex; index < maxIndex + 1; index++){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
		}
		totalCount++;
	}

	//TODO: Shall the on-grid-points be re-checked? We probably should leave them out and take them as valid

	// expand into negative direction if treshold wasn't met yet
	int index = minIndex - 1;
	while((((double)validCount) / ((double)totalCount)) >= TRESHOLD2){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
			smallestValidIndex = index;
		}
		totalCount++;
		index--;
	}

	// remove outermost invalid points
	totalCount = highestValidIndex - smallestValidIndex + 1;

	// expand into positive direction if treshold wasn't met yet
	index = maxIndex + 1;
	while((((double)validCount) / ((double)totalCount)) >= TRESHOLD2){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
			highestValidIndex = index;
		}
		totalCount++;
		index++;
	}

	// remove outermost invalid points
	totalCount = highestValidIndex - smallestValidIndex + 1;

	double ratio = ((double)validCount) / ((double)totalCount);

	return ratio;
}

// project points on the line through referencePoint in the direction of candidateVector
vector<Vector3d> LatticeDetector::projectPointsOnLine(Vector3d const &referencePoint, Vector3d const &candidateVector){

	vector<Vector3d> projectedPoints = vector<Vector3d>();
	projectedPoints.reserve(reconstructedPoints.size());

	std::vector<Vector3d>::iterator pointsIt;

	for(pointsIt = reconstructedPoints.begin(); pointsIt != reconstructedPoints.end(); ++pointsIt){

		// point P, reference point R
		Vector3d point = (*pointsIt);

		// calculate RP
		Vector3d vectorRefToPoint = point-referencePoint;

		// project RP onto candidateVector C: RP' = (RP*C/(C*C))*C
		Vector3d projectionOnCandidateVector = (vectorRefToPoint.dot(candidateVector) / candidateVector.dot(candidateVector) ) * candidateVector;

		// find P' by adding R and RP'
		Vector3d projectedPoint = referencePoint+projectionOnCandidateVector;

		projectedPoints.push_back(projectedPoint);
	}

	return projectedPoints;
}

// returns the indices of the two outermost grid points that have reconstructed points on them
vector<int> LatticeDetector::getOutermostOnGridPointIndices(vector<Vector3d> const &projectedPoints, Vector3d const &referencePoint, Vector3d const &candidateVector){

	int minIndex = 0;
	int maxIndex = 0;

	double treshold = candidateVector.norm() * TRESHOLD1;

	std::vector<Vector3d>::const_iterator projectedPointsIt;

	// This iterator is advanced in the same way as projectedPointsIt, but explicitly
	std::vector<Vector3d>::iterator originalPointsIt;

	originalPointsIt = reconstructedPoints.begin();

	// Iterate through the (projected) points
	for (projectedPointsIt = projectedPoints.begin(); projectedPointsIt != projectedPoints.end(); ++projectedPointsIt){

		Vector3d projectedPoint = *projectedPointsIt;
		Vector3d originalPoint = *originalPointsIt;

		// projectedPoint = referencePoint + candidateVector*factor
		// factor = |(projectedPoint - referencePoint)| / |candidateVector|
		double factor = (projectedPoint - referencePoint).norm() / candidateVector.norm();

		// the sign of the factor is missing, so test for it separately
		double diffWithPlus = ((referencePoint + candidateVector*factor) - projectedPoint).norm();
		double diffWithMinus = ((referencePoint - candidateVector*factor) - projectedPoint).norm();
		if (diffWithMinus < diffWithPlus){
			factor = -factor;
		}

		// Find the most adjacent two grid points
		int previousGridPointIndex = floor(factor);
		int nextGridPointIndex = ceil(factor);

		Vector3d previousGridPoint = referencePoint + candidateVector*previousGridPointIndex;
		Vector3d nextGridPoint = referencePoint + candidateVector*nextGridPointIndex;

		int finalGridPointIndex;

		// Check whether the original point is on-grid
		if((originalPoint - previousGridPoint).norm() < treshold){ // Point is on-grid on the previous grid point
			finalGridPointIndex = previousGridPointIndex;
		}
		else if((originalPoint - nextGridPoint).norm() < treshold){ // Point is on-grid on the next grid point
			finalGridPointIndex = nextGridPointIndex;
		}
		else{
			++originalPointsIt;
			continue;
		}

		// Check whether it is the outermost found point in one direction
		if(finalGridPointIndex < minIndex){
			minIndex = finalGridPointIndex;
		}
		else if(finalGridPointIndex > maxIndex){
			maxIndex = finalGridPointIndex;
		}

		++originalPointsIt;
	}

	vector<int> outermostOnGridPointIndices = vector<int>();

	outermostOnGridPointIndices.push_back(minIndex);
	outermostOnGridPointIndices.push_back(maxIndex);

	return outermostOnGridPointIndices;
}

vector<Vector3d> LatticeDetector::calculateLatticeBoundary(Vector3d const &latticeVector1, Vector3d const &latticeVector2){

	std::vector<Vector3d>::iterator referencePointsIt;

	int finalArea = -1;
	vector<Vector3d> finalLatticeBoundary;

	for(referencePointsIt = reconstructedPoints.begin(); referencePointsIt != reconstructedPoints.end(); ++referencePointsIt){
		Vector3d referencePoint = (*referencePointsIt);
		vector<Vector3d> latticeBoundary = vector<Vector3d>();
		int area = 0;
		latticeBoundaryForReferencePoint(referencePoint, latticeVector1, latticeVector2, latticeBoundary, area);

		if (area > finalArea){
			finalLatticeBoundary = latticeBoundary;
			finalArea = area;
		}
	}

	return finalLatticeBoundary;
}

void LatticeDetector::latticeBoundaryForReferencePoint(Vector3d const &referencePoint, Vector3d const &latticeVector1, Vector3d const &latticeVector2, vector<Vector3d> &latticeBoundaryOut, int &areaOut){

	// assume latticeVector1 to point towards right, latticeVector2 to point towards up

	Vector3d lowerLeft = referencePoint;
	Vector3d upperRight = referencePoint;

	int width = 0;
	int height = 0;

	int unexpandableCount = 0;

	while(unexpandableCount<4){
		// expand right

		if(validLine(referencePoint, upperRight + latticeVector1, -latticeVector2, height)){
			width++;
			upperRight = upperRight + latticeVector1;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}

		// expand top
		if(validLine(referencePoint, upperRight + latticeVector2, -latticeVector1, width)){
			height++;
			upperRight = upperRight + latticeVector2;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}

		// expand left
		if(validLine(referencePoint, lowerLeft - latticeVector1, latticeVector2, height)){
			width++;
			lowerLeft = lowerLeft - latticeVector1;
			unexpandableCount = 0;
		}
			else{
			unexpandableCount++;
		}

		// expand bottom
		if(validLine(referencePoint, lowerLeft - latticeVector2, latticeVector1, width)){
			height++;
			lowerLeft = lowerLeft - latticeVector2;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}
	}

	latticeBoundaryOut.push_back(lowerLeft);
	latticeBoundaryOut.push_back(upperRight);

	areaOut = (width + 1) * (height + 1);
}

bool LatticeDetector::validLine(Vector3d const &referencePoint, Vector3d const &anchorPoint, Vector3d const &directionVector, int length){

	int validCount = 0;
	int totalCount = length+1;

	for (int i=0; i<=length; i++){
		Vector3d pointToTest = anchorPoint + directionVector*i;
		if(isPointValid(referencePoint, pointToTest)){
			validCount++;
		}
	}

	bool valid = ((double)validCount)/totalCount >= TRESHOLD2;

	return valid;
}

//=================================================================================
// Nektarios's

bool LatticeDetector::isPointValid(Vector3d const &referencePoint, Vector3d const &pointToTest){
	bool valid = compareSiftFronto(referencePoint, pointToTest, this->plane,
			//TriangulatedPoint Xref,
			this->inpManager->getK(), this->inpManager->getCamPoses(), this->inpManager->getViewIds(),
			this->inpManager->getImgNames());

	//TODO Change, only for testing purposes
	//bool valid = (std::find(validGridPoints.begin(), validGridPoints.end(), pointToTest) != validGridPoints.end());

	return valid;
}

bool LatticeDetector::isIntegerCombination(int i,vector<Vector3d> candidatesInOrder,vector<bool> valid){


	//int validElems = std::accumulate(valid.begin()+i+1, valid.end(), 0);

	/*
	Matrix<double,3,Eigen::Dynamic> A(3,validElems+1);
	int v = 0;
	for (int j=i+1; j<candidatesInOrder.size(); j++){
		if (valid.at(j)){
			A.block<3,1>(0,v) = candidatesInOrder.at(j);
			v++;
		}
	}
	A.block<3,1>(0,v) = -candidatesInOrder.at(i);

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeFullV);
	VectorXd combination = svd.matrixV().block<4,1>(0,3);//<sizeRows,sizeCols>(beginRow,beginCol)

	combination = combination/combination[validElems];


	//only works for >=3 elements
	//VectorXd combination = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(candidatesInOrder[i]);

	bool isIntComb = true;
	for (int i=0; i<candidatesInOrder.size()-i-1; i++){
		if ((fmod(combination[i], 1) > 0.00001) &&  (fmod(combination[i],1) < 0.99999)){
			isIntComb = false;
		}
	}

	std::cout << "is integer comb" << std::endl;

	return isIntComb;
*/

	Matrix<double,3,Eigen::Dynamic> A(3,3);
	VectorXd solution;
	A.block<3,1>(0,2) = -candidatesInOrder.at(i);

	for (int j=i+1;j<candidatesInOrder.size();j++){
		if (!valid.at(j))
			continue;
//		Vector3d n =  candidatesInOrder[i].cross(candidatesInOrder[j]);
		A.block<3,1>(0,0) = candidatesInOrder.at(j);

		for (int k=j+1;j<candidatesInOrder.size();k++){
			if (!valid.at(k))
				continue;

			A.block<3,1>(0,1) = candidatesInOrder.at(k);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeFullV);

			solution = svd.matrixV().block<3,1>(0,2);//<sizeRows,sizeCols>(beginRow,beginCol)
			solution = solution/solution[2];
			if ((fmod(solution[0], 1) > 0.001) ||  (fmod(solution[0],1) < 0.999) ||
					(fmod(solution[1], 1) > 0.001) ||  (fmod(solution[1],1) < 0.999))
				return true;

		}
	}

	return false;
}


//given the set of candidate vectors, find the best two basis vectors
vector<Vector3d> LatticeDetector::getFinalBasisVectors(vector<Vector3d> candidateVectors){

	// Get the scores
	vector<double> scores = this->validateCandidateVectors(candidateVectors);

	std::vector<int> indices(candidateVectors.size());
	std::iota(indices.begin(), indices.end(), 0); //0 is the starting number.


	//Sort the indices according to vector's length
	std::sort(indices.begin(), indices.end(),
			[&](const int& a, const int& b) {
	return (candidateVectors.at(a).squaredNorm() > candidateVectors.at(b).squaredNorm());
	}
	       );

	//reorder
	vector<double> scoresInOrder(scores.size());
	vector<Vector3d> candidatesInOrder(scores.size());
	for (int i = 0; i < scores.size(); i++){
		scoresInOrder[i] = scores[indices[i]];
		candidatesInOrder[i] = candidateVectors[indices[i]];
	}

	vector<bool> valid(scores.size());
	std::fill(valid.begin(),valid.end(),true);
	for (int i = 0; i < scores.size(); i++){

		//1st check: if integer combination
		if (scores.size() - (i+1) >= 2){
			if (isIntegerCombination(i,candidatesInOrder,valid)){
				valid.at(i) = false;
				continue;
			}
		}
		//2nd check: if inline with other and lower score
		for (int j=i+1;j<candidatesInOrder.size();j++){
			if (!valid[j])
				continue;
			double cosangle = (candidatesInOrder[i].dot(candidatesInOrder[j]))/
					( sqrt(candidatesInOrder[i].squaredNorm())*sqrt(candidatesInOrder[j].squaredNorm()) ) ;
			if (abs(cosangle - 1) < 0.01){
				if (scoresInOrder[i] < scoresInOrder[j])
					valid[i] = false;
				else
					valid[j] = false;
			}
		}
	}

	//get 2 remaining valid-with max score
	vector<Vector3d> basis;
	for (int i=scoresInOrder.size()-1;i>=0;i--){
		if (valid[i]){
		    basis.push_back(candidatesInOrder[i]);
		}
		if (basis.size()==2)
			break;
	}

	// Hopefully 2 basis vectors are returned

	return basis;
}

