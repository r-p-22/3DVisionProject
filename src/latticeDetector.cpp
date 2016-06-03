#include "latticeDetector.h"

#include <algorithm>
#include <numeric>
#include <math.h>
#include "3dtools.h"


LatticeDetector::LatticeDetector(vector<Vector3d> const &aPoints, Vector4d const &aPlane, inputManager* aInputManager){
	points = aPoints;
	plane = aPlane;
	inpManager = aInputManager;
}

LatticeDetector::~LatticeDetector(){

}

Vector3d LatticeDetector::translationVector(Vector3d const &point1, Vector3d const &point2){

	Vector3d translation = point1-point2;
	return translation;

}

vector<Vector3d> LatticeDetector::calculateCandidateVectors(bool naive){

	vector<Vector3d> rawCandidates = vector<Vector3d>();

	int length = points.size();

	for (int i=0; i<length; i++){
		for (int j=i+1; j<length; j++){
			Vector3d translation = translationVector(points[i], points[j]);
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

	//remove very small finalcandidates
	vector<Vector3d>::iterator i = finalCandidates.begin();
	while (i != finalCandidates.end())
	{
	    bool isActive = (*i).norm() <= this->VECTOR_DISTANCE;
	    if (isActive)
	    {
	    	//finalCandidates.erase(i++);  // alternatively,
	    	i = finalCandidates.erase(i);
	    }
	    else
	    {
	        ++i;
	    }
	}

	// return the two best candidates for naive solution
	if(naive && finalCandidates.size() > 2){
		vector<Vector3d> finalNaiveCandidates = vector<Vector3d>();

		vector<int>::iterator scoresIt;

		int bestScore = 0;
		int secondBestScore = 0;
		int bestIndex = 0;
		int secondBestIndex = 0;

		for(scoresIt = scores.begin(); scoresIt != scores.end(); ++scoresIt){
			int score = (*scoresIt);

			if(score >= bestScore){

				secondBestScore = bestScore;
				secondBestIndex = bestIndex;

				bestScore = score;
				bestIndex = std::distance(scores.begin(), scoresIt);
			}
			else if(score >= secondBestScore){
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

void LatticeDetector::combineCandidates(list<list<Vector3d> > const &clusteredCandidates, vector<Vector3d> &combinedCandidates, vector<int> &scores){

	list<list<Vector3d> >::const_iterator clusterItOuter;
	list<Vector3d>::const_iterator clusterItInner;

	for (clusterItOuter = clusteredCandidates.begin(); clusterItOuter != clusteredCandidates.end(); ++clusterItOuter){

		int score = 0;
		Vector3d combinedCandidateSum = Vector3d(0,0,0);


		for(clusterItInner = (*clusterItOuter).begin(); clusterItInner != (*clusterItOuter).end(); ++clusterItInner){
			Vector3d nextCandidate = *clusterItInner;

			// If not the first candidate, align the candidate with the orientation of the candidates processed until now
			if (clusterItInner!=(*clusterItOuter).begin()){
				Vector3d combinedCandidatePreview = combinedCandidateSum/score;

				//if pointing into opposite direction, revert the candidate
				if((nextCandidate+combinedCandidatePreview).norm() < (nextCandidate-combinedCandidatePreview).norm()){
					nextCandidate = nextCandidate*(-1);
				}
			}

			// Sum the candidates in a cluster
			combinedCandidateSum = combinedCandidateSum + nextCandidate;
			// Count the candidates of the cluster
			score++;
		}

		// Average the candidates in a cluster
		Vector3d combinedCandidate = combinedCandidateSum/score;

		scores.push_back(score);
		combinedCandidates.push_back(combinedCandidate);
	}
}

// If similar and in line: 1
// If similar and opposite orientation: 2
// If not similar: -1
int LatticeDetector::vectorsAreSimilar(Vector3d const &vector1, Vector3d const &vector2, double treshold){

	// Note: Difference of two vectors is the intuitive vector to check. However, if the vectors
	// are similar but point in opposite directions, the difference will be huge only for the orientation
	// reasons. So check as well the sum, which is the difference of the two vectors with one being turned 180°
	// by mirroring it in the origin (vector1-(vector2*(-1))

	if(((vector1 - vector2).norm() <= treshold)){
		return 1;
	}
	else if(((vector1 + vector2).norm() <= treshold)){
		return 2;
	}
	else{
		return -1;
	}
}

// Clusters candidates together if they are similar. Clusters are lists of candidate vectors, themselves stored in a list
list<list<Vector3d> > LatticeDetector::clusterCandidates(vector<Vector3d> const &candidates){

	vector<Vector3d>::const_iterator candidateIt;
	list<list<Vector3d> >::iterator clusterItOuter;
	list<Vector3d>::iterator clusterItInner;

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
				if (vectorsAreSimilar(*candidateIt, *clusterItInner, VECTOR_DISTANCE) > 0){
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

	vector<Vector3d>::const_iterator candidateIt;

	// store the score of every candidate vector
	for(candidateIt = candidateVectors.begin(); candidateIt != candidateVectors.end(); ++candidateIt){
		double score = validateVector(*candidateIt);
		scores.push_back(score);
	}

	return scores;

}

double LatticeDetector::validateVector(Vector3d const &candidateVector){

	vector<Vector3d>::iterator pointsIt;

	double imageValidationScore = 0;

	// sum up the score (ration between valid and invalid grid points) of every reference point
	for(pointsIt = points.begin(); pointsIt != points.end(); ++pointsIt){

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

	int smallestValidIndex = 0;
	int highestValidIndex = 0;

	// initialize with -1 as the reference point will raise it to 0
	int validCount = -1;

	// check whether points between the outermost on grid points are valid
	for (int index = minIndex; index <= maxIndex; index++){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
			if(index < smallestValidIndex){
				smallestValidIndex = index;
			}
			if(index > highestValidIndex){
				highestValidIndex = index;
			}
		}
	}

	// remove outermost invalid points
	int totalCount = highestValidIndex - smallestValidIndex;

	// expand into negative direction if treshold wasn't met yet
	int index = smallestValidIndex - 1;
	while((totalCount == 0) || ((((double)validCount) / ((double)totalCount)) >= TRESHOLD2)){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
			smallestValidIndex = index;
		}
		totalCount++;
		index--;
	}

	// remove outermost invalid points
	totalCount = highestValidIndex - smallestValidIndex;

	// expand into positive direction if treshold wasn't met yet
	index = maxIndex + 1;
	while((totalCount == 0) || ((((double)validCount) / ((double)totalCount)) >= TRESHOLD2)){
		Vector3d pointToTest = referencePoint + candidateVector*index;
		if (isPointValid(referencePoint, pointToTest)){
			validCount++;
			highestValidIndex = index;
		}
		totalCount++;
		index++;
	}

	// remove outermost invalid points
	totalCount = highestValidIndex - smallestValidIndex;

	double ratio = 0;

	if (totalCount > 0){
		ratio = ((double)validCount) / ((double)totalCount);
	}

	return ratio;
}

// project points on the line through referencePoint in the direction of candidateVector
vector<Vector3d> LatticeDetector::projectPointsOnLine(Vector3d const &referencePoint, Vector3d const &candidateVector){

	vector<Vector3d> projectedPoints = vector<Vector3d>();
	projectedPoints.reserve(points.size());

	vector<Vector3d>::iterator pointsIt;

	for(pointsIt = points.begin(); pointsIt != points.end(); ++pointsIt){

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

	vector<Vector3d>::const_iterator projectedPointsIt;

	// This iterator is advanced in the same way as projectedPointsIt, but explicitly
	vector<Vector3d>::iterator originalPointsIt;

	originalPointsIt = points.begin();

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
		if(pointEqualsGridPoint(originalPoint, previousGridPoint, candidateVector)){ // Point is on-grid on the previous grid point
			finalGridPointIndex = previousGridPointIndex;
		}
		else if(pointEqualsGridPoint(originalPoint, nextGridPoint, candidateVector)){ // Point is on-grid on the next grid point
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

bool LatticeDetector::pointEqualsGridPoint(Vector3d const &point, Vector3d const &gridPoint, Vector3d const &vector){

	double treshold = vector.norm() * TRESHOLD1;

	bool equals = (point - gridPoint).norm() < treshold;

	return equals;
}

void LatticeDetector::calculateLatticeBoundary(Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &cornerOut, int &widthOut, int &heightOut){

	vector<Vector3d>::iterator referencePointsIt;

	int finalArea = -1;
	int finalWidth = 0;
	int finalHeight = 0;
	Vector3d finalCorner = Vector3d(0,0,0);

	for(referencePointsIt = points.begin(); referencePointsIt != points.end(); ++referencePointsIt){
		Vector3d referencePoint = (*referencePointsIt);

		int width;
		int height;
		Vector3d corner;
		latticeBoundaryForReferencePoint(referencePoint, latticeVector1, latticeVector2, corner, width, height);

		int area = (width + 1) * (height + 1);

		if (area > finalArea){
			finalArea = area;
			finalWidth = width;
			finalHeight = height;
			finalCorner = corner;
		}
	}

	widthOut = finalWidth;
	heightOut = finalHeight;
	cornerOut = finalCorner;

}

void LatticeDetector::latticeBoundaryForReferencePoint(Vector3d const &referencePoint, Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &cornerOut, int &widthOut, int &heightOut){

	// assume latticeVector1 to point towards "right", latticeVector2 to point towards "up"

	Vector3d lowerLeft = referencePoint;
	Vector3d upperRight = referencePoint;

	widthOut = 0;
	heightOut = 0;

	int unexpandableCount = 0;

	while(unexpandableCount<4){

		// expand right
		if(validLine(referencePoint, upperRight + latticeVector1, -latticeVector2, heightOut)){
			widthOut++;
			upperRight = upperRight + latticeVector1;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}

		// expand top
		if(validLine(referencePoint, upperRight + latticeVector2, -latticeVector1, widthOut)){
			heightOut++;
			upperRight = upperRight + latticeVector2;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}

		// expand left
		if(validLine(referencePoint, lowerLeft - latticeVector1, latticeVector2, heightOut)){
			widthOut++;
			lowerLeft = lowerLeft - latticeVector1;
			unexpandableCount = 0;
		}
			else{
			unexpandableCount++;
		}

		// expand bottom
		if(validLine(referencePoint, lowerLeft - latticeVector2, latticeVector1, widthOut)){
			heightOut++;
			lowerLeft = lowerLeft - latticeVector2;
			unexpandableCount = 0;
		}
		else{
			unexpandableCount++;
		}
	}

	cornerOut = lowerLeft;

}

bool LatticeDetector::validLine(Vector3d const &referencePoint, Vector3d const &anchorPoint, Vector3d const &directionVector, int length){

	int validCount = 0;
	int totalCount = length+1;

	for (int i = 0; i <= length; i++){
		Vector3d pointToTest = anchorPoint + directionVector*i;
		if(isPointValid(referencePoint, pointToTest)){
			validCount++;
		}
	}

	bool valid = ((double)validCount)/(double(totalCount)) >= TRESHOLD2;

	return valid;
}


vector<pair<int, vector<int> > > LatticeDetector:: getOnGridIndices(vector<int> const &inputIndices, LatticeStructure const &lattice){

	// assume latticeVector1 to point in positive "width" direction, latticeVector2 to point in positive "height" direction

	vector<pair<int, vector<int> > > pointsToIndices = vector<pair<int, vector<int> > >();

	vector<Vector3d> gridPoints = vector<Vector3d>();
	vector<vector<int> > gridPointIndices = vector<vector<int> >();

	// get lattice parameters
	int width = lattice.width;
	int height = lattice.height;
	Vector3d latticeVector1 = lattice.basisVectors[0];
	Vector3d latticeVector2 = lattice.basisVectors[1];
	Vector3d corner = lattice.corner;

	// get all grid points together with their indices
	for (int i = 0; i <= width; i++){
		for (int j = 0; j <= height; j++){
			Vector3d gridPoint = corner + latticeVector1*i + latticeVector2*j;

			vector<int> indices = vector<int>();
			indices.push_back(i);
			indices.push_back(j);

			gridPoints.push_back(gridPoint);
			gridPointIndices.push_back(indices);
		}
	}

	// get all reconstructed points in basis of lattice
	vector<Vector3d> reconstructedPointsInLatticeBasis = changeToLatticeBasis(points, latticeVector1, latticeVector2);
	vector<Vector3d> gridPointsInLatticeBasis = changeToLatticeBasis(gridPoints, latticeVector1, latticeVector2);

	int numberOfReconstructedPoints = points.size();
	int numberOfGridPoints = gridPoints.size();

	// For every grid point, find the reconstructed points lying on it and select the closest one
	for (int i = 0; i < numberOfGridPoints; i++){

		Vector3d gridPoint = gridPoints[i];
		Vector3d gridPointInLatticeBasis = gridPointsInLatticeBasis[i];

		double minDistance ;
		vector<int> minIndices;

		bool pointFound = false;

		int jmin = -1;
		for(int j = 0; j < numberOfReconstructedPoints; j++){

			Vector3d reconstructedPointInLatticeBasis = reconstructedPointsInLatticeBasis[j];
			Vector3d reconstructedPoint = points[j];

			// Determine distance to the grid point
			double distance = (gridPoint-reconstructedPoint).norm();

			// Determine whether the point is close to a grid point
			// NOTE: "Close" has to be determined here wrt two different basis vectors, so the trivial check is not applicable.
			// To combine both tresholds, we consider points as "close" if they lie in a parallelogram around the grid point,
			// in a manner similar as depicted (on the left, depiction of the two basis vectors, parallelogram has the same shape,
			// but only expands to TRESHOLD1*correspondingLatticeVector into each direction around the grid point, so the "valid"
			// parallelogram area is bound by two vectors that correspond to TRESHOLD1*2*latticeVector1 and TRESHOLD1*2*latticeVector2
			//
			//         /     ______
			//        /     /     /
			//       /     /  .  /
			//      /     /_____/
			//	   /
			//    /
			//   /___________
			//
			// Thanks to a change to lattice basis coordinate system, the check whether a point is in that area can be done easily by
			// checking each of the first two coordinates of the distance vector separately.

			Vector3d distanceVector = reconstructedPointInLatticeBasis-gridPointInLatticeBasis;
			bool onGridInWidth = abs(distanceVector[0]) < TRESHOLD1;
			bool onGridInHeight = abs(distanceVector[1]) < TRESHOLD1;
			bool onGrid = onGridInWidth && onGridInHeight;

			// if the point is on the grid and it is either the first found point or it is closer than all other found points, update
			if(onGrid && (!pointFound || (distance < minDistance))){
			//	if((distance < minDistance)){
				pointFound = true;
				minDistance = distance;
				minIndices = gridPointIndices[i];
				jmin = j;
			}
		}

		// If there was a reconstructed point found for that particular grid point, save it together with the grid indices
		if(pointFound){
			pair<int,vector<int> > indexPair = pair<int,vector<int> >();
			indexPair.first = inputIndices[jmin];
			indexPair.second = minIndices;
			pointsToIndices.push_back(indexPair);
		}
	}

	return pointsToIndices;
}


vector<Vector3d> LatticeDetector::changeToLatticeBasis(vector<Vector3d> const &inputPoints, Vector3d const &latticeVector1, Vector3d const &latticeVector2){

	vector<Vector3d> coordinatesInLatticeBasis = vector<Vector3d>();

	// Compose the transformationMatrix from lattice basis to canonical basis
	Matrix3d transformationToCanonical = Matrix3d();
	transformationToCanonical.col(0) = latticeVector1;
	transformationToCanonical.col(1) = latticeVector2;
	transformationToCanonical.col(2) = latticeVector1.cross(latticeVector2);

	// Inverse is the transformationMatrix from canonical basis to lattice basis
	Matrix3d transformationToLatticeBasis = transformationToCanonical.inverse();

	vector<Vector3d>::const_iterator pointsIt;

	for(pointsIt = inputPoints.begin(); pointsIt < inputPoints.end(); ++pointsIt){
		Vector3d transformedCoordinates = transformationToLatticeBasis * (*pointsIt);
		coordinatesInLatticeBasis.push_back(transformedCoordinates);
	}

	return coordinatesInLatticeBasis;
}


//=================================================================================
// Nektarios's

bool LatticeDetector::isPointValid(Vector3d const &referencePoint, Vector3d const &pointToTest){

	bool valid = compareSiftFronto(referencePoint, pointToTest, this->plane,
			this->inpManager->getK(), this->inpManager->getCamPoses(), this->inpManager->getViewIds(),
			this->inpManager->getImgNames());

	return valid;
}

bool LatticeDetector::isIntegerCombination(int i,vector<Vector3d>& candidatesInOrder,vector<bool>& valid){

	Matrix<double,3,3> A;
	Vector3d solution;
	A.block<3,1>(0,2) = -candidatesInOrder.at(i);

	for (size_t j=i+1; j<candidatesInOrder.size(); j++){
		if (!valid.at(j))
			continue;

//		Vector3d n =  candidatesInOrder[i].cross(candidatesInOrder[j]);
		A.block<3,1>(0,0) = candidatesInOrder.at(j);

		for (size_t k=j+1; k<candidatesInOrder.size(); k++){
			if (!valid.at(k))
				continue;

			A.block<3,1>(0,1) = candidatesInOrder.at(k);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeFullV);

			solution = svd.matrixV().block<3,1>(0,2);//<sizeRows,sizeCols>(beginRow,beginCol)
			solution = solution/solution[2];

			solution[1] = abs(solution[1]);
			solution[0] = abs(solution[0]);
			//check if integer combination of self, i.e. all elems close to zero
			if ( ((solution[0]) < 0.03) && ((solution[1]) < 0.03) )
				continue;
			if ( ( (fmod(solution[0], 1) < 0.012) ||  (fmod(solution[0],1) > 0.988) ) &&
					( (fmod(solution[1], 1) < 0.012) ||  (fmod(solution[1],1) > 0.988) ) )
				return true;
		}
	}

	return false;
}


//given the set of candidate vectors, find the best two basis vectors
vector<Vector3d> LatticeDetector::getFinalBasisVectors(vector<Vector3d> const &candidateVectors){


	int N = candidateVectors.size();

	cout << "initial candvecs: "<< N << endl;
	vector<int> indices(N);
	std::iota(indices.begin(), indices.end(), 0); //0 is the starting number.


	//Sort the indices according to vector's length
	std::sort(indices.begin(), indices.end(),
			[&](const int& a, const int& b) {
	return (candidateVectors.at(a).squaredNorm() > candidateVectors.at(b).squaredNorm());
	}
		   );

	//reorder
	vector<Vector3d> candidatesInOrder(N);
	for (int i = 0; i < N; i++){
		candidatesInOrder[i] = candidateVectors[indices[i]];
	}
	vector<bool> valid(N);
	std::fill(valid.begin(),valid.end(),true);

	for (int i=0; i<N;i++){
		cout << candidatesInOrder[i] << endl;
		cout <<"--"<<endl;
	}

	for (int i = 0; i < N; i++){

		//1st check: if integer combination
		if (N - (i+1) >= 2){
			bool a;
			a = isIntegerCombination(i,candidatesInOrder,valid);
			if (a){
				valid.at(i) = false;
				continue;
			}
		}
	}


	//remove invalid candidatesInOrder
	vector<Vector3d>::iterator i = candidatesInOrder.begin();
	vector<bool>::iterator v = valid.begin();
	while (v != valid.end())
	{
		bool isActive = !(*v);
		if (isActive)
		{
			//finalCandidates.erase(i++);  // alternatively,
			i = candidatesInOrder.erase(i);
			v = valid.erase(v);
		}
		else
		{
			++i;
			++v;
		}
	}

	N = candidatesInOrder.size();

	cout << "remaining candvecs after int.comb: "<< N << endl;


	// Get the scores
	vector<double> scoresInOrder = this->validateCandidateVectors(candidatesInOrder);

	for (int i=0; i<N;i++){
			cout << candidatesInOrder[i] << endl;
			cout << scoresInOrder[i] << endl;
			cout << valid[i] << endl;
		}

	//cout << "score calculated" << endl;

	for (int i = 0; i < N; i++){

		//2nd check: if inline with other and lower score
		for (int j=i+1;j<N;j++){
			if (!valid[j])
				continue;
			double cosangle = (candidatesInOrder[i].dot(candidatesInOrder[j]))/
					( sqrt(candidatesInOrder[i].squaredNorm())*sqrt(candidatesInOrder[j].squaredNorm()) ) ;
			if ( (1-abs(cosangle)) < 0.04){
				if (scoresInOrder[i] < scoresInOrder[j]) {
					valid[i] = false;
					break;
				} else {
					valid[j] = false;
				}

			}
		}
	}

	vector<Vector3d> basis;

	int num = 0;
	int a;
	while (num < 2){
		//get maximum element
		a = std::distance(scoresInOrder.begin(), std::max_element(scoresInOrder.begin(), scoresInOrder.end()));
		if (scoresInOrder[a] <= 0){
			cout << "==================================" << endl;
			cout << "===ERROR: NO TWO BASIS VECS FOUND===" << endl;
			cout << "==================================" << endl;
			break;
		}
		if (valid[a]){
			basis.push_back(candidatesInOrder[a]);
			num++;
		}
		scoresInOrder[a] = -1;

	}

	// Hopefully 2 basis vectors are returned

	return basis;
}

