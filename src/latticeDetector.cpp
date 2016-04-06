#include "latticeDetector.h"

/*
 *
 *
 * # Find translation between all point pairs (translationVector)
 *
 * # "Cluster" translation vectors: treshold in difference, weighted average?
 *
 * # take 2 highest peaks naively
 *
 * # or take all peaks as candidates
 *
 * # verify basis vector
 *
 *  	-select a 3D point as reference
 *  	- predict grid points on line through reference point along vector
 *  	- compare SIFT descriptors of grid points to ref. point:
 *  		if alpha > 2phi -> valid
 *  			# SIFT descriptor is computed in most frontoparallel image
 *  				line connecting camera center and point is closest to plane normal
 *  	- check every reference point like that
 *  	- image validation score: SUM of ratios over all reference points
 *
 *  	PER REFERENCE POINT: number of valid grid points / all grid points
 *
 *  	- search for two farthest reconstructed on-grid points on both sides: on-grid if
 *  		dist to a certain grid point is smaller than T1=10% of basis vector length
 *  	- move them away until T2=50% of grid points within are invalid
 *  	- trim invalid grid points on both ends
 *
 *  # choose two basis vectors:
 *
 *  	- sort in descending order of lengths
 *  	- discard vector if integer combination of rest of the queue
 *  	- vectors along same direction: keep the one with higher score
 *  	- take two highest-score vectors of remaining ones
 */



Vector3d LatticeDetector::translationVector(Vector3d &point1, Vector3d &point2){

	Vector3d translation = point1-point2;
	return translation;

}

vector<Vector3d> LatticeDetector::calculateCandidateVectors(bool naive){

	std::vector<Vector3d> rawCandidates = vector<Vector3d>();

	int length = points.size();

	for (int i=0; i<length; i++){
		for (int j=i+1; j<length; j++){
			Vector3d translation = translationVector(rawCandidates[i], rawCandidates[j]);
			rawCandidates.push_back(translation);
		}
	}

	// Init data structures
	list<list<Vector3d> > clusteredCandidates = list<list<Vector3d> >();
	vector<Vector3d> finalCandidates = vector<Vector3d>();
	vector<int> scores = vector<int>();

	// Cluster candidates that are similar
	clusterCandidates(rawCandidates, clusteredCandidates);

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

void LatticeDetector::combineCandidates(list<list<Vector3d> > &clusteredCandidates, vector<Vector3d> &finalCandidates, vector<int> &scores){

	std::list<list<Vector3d> >::iterator clusterItOuter;
	std::list<Vector3d>::iterator clusterItInner;

	for (clusterItOuter = clusteredCandidates.begin(); clusterItOuter != clusteredCandidates.end(); ++clusterItOuter){

		int score = 0;
		Vector3d finalCandidate = Vector3d(0,0,0);

		for(clusterItInner = (*clusterItOuter).begin(); clusterItInner != (*clusterItOuter).end(); ++clusterItInner){
			// Average the candidates in a cluster
			finalCandidate = finalCandidate + (*clusterItInner);
			// Count the candidates of the cluster
			score++;
		}

		// Average the candidates in a cluster
		finalCandidate = finalCandidate/score;

		scores.push_back(score);
		finalCandidates.push_back(finalCandidate);
	}
}

bool LatticeDetector::vectorsAreSimilar(Vector3d vector1, Vector3d vector2){

	// Note: Difference of two vectors is the intuitive vector to check. However, if the vectors
	// are similar but point in opposite directions, the difference will be huge only for the orientation
	// reasons. So check as well the sum, which is the difference of the two vectors with one being turned 180°
	// by mirroring it in the origin (vector1-(vector2*(-1))

	if(((vector1-vector2).norm()<=VECTOR_DISTANCE) || ((vector1+vector2).norm()<=VECTOR_DISTANCE)){
		return true;
	}
	else{
		return false;
	}
}

// Clusters candidates together if they are similar. Clusters are lists of candidate vectors, themselves stored in a list
void LatticeDetector::clusterCandidates(vector<Vector3d> &candidates, list<list<Vector3d> > &clusteredCandidates){

	std::vector<Vector3d>::iterator candidateIt;
	std::list<list<Vector3d> >::iterator clusterItOuter;
	std::list<Vector3d>::iterator clusterItInner;

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

}
