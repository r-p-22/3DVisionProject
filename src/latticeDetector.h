/*
 * latticeDetector.h
 *
 *  Created on: Apr 4, 2016
 *      Author: valaina
 */

#ifndef SRC_LATTICEDETECTOR_H_
#define SRC_LATTICEDETECTOR_H_

#include <list>
#include <Eigen/Dense>

#include "inputManager.h"

using namespace Eigen;
using namespace std;

class LatticeDetector{

public:

	/*!
	 * The constructor for the lattice detector
	 *
	 * @param[in] aPoints		The points to fit a lattice.
	 * @param[in] aPlane		The plane on which the points lie.
	 * @param[in] aInputManager	The input manager that gives access to all given initial data
	 */
	LatticeDetector(vector<Vector3d> const &aPoints, Vector4d const &aPlane, inputManager* aInputManager);

	/*!
	 * The destructor for the lattice detector
	 */
	~LatticeDetector();

	/*!
	 * Calculates the candidate basis vectors for the given points. If naive is true, only the two candidate vectors that occurred most often are returned.
	 *
	 * @param[in] naive		If false, all candidate vectors are returned. If true, only the two candidate vectors that occurred most often are returned.
	 * @return	The candidate basis vectors for the given points, either all or only the two that occurred most often, depending on "naive".
	 */
	vector<Vector3d> calculateCandidateVectors(bool naive);

	/*!
	 * Selects the best two basis vectors from a vector of candidates. A message is shown if not two vectors were found. Then the resulting vector will have < 2 entries, so an extra check outside the function must be made.
	 *
	 * @param[in] candidateVectors	The candidate vectors from which the best two basis vectors should be chosen
	 * @return	The best two basis vectors, selected from the available candidates.
	 */
	vector<Vector3d> getFinalBasisVectors(vector<Vector3d> const &candidateVectors);

	/*!
	 * Calculates the lattice boundary in terms of corner, width and height.
	 * The lattice is thought to span a parallelogram, originating at corner, with sides latticeVector1*width and latticeVector2*height
	 *
	 * @param[in] latticeVector1	The first basis vector of the lattice.
	 * @param[in] latticeVector2	The second basis vector of the lattice.
	 * @param[out] cornerOut		The corner of the lattice.
	 * @param[out] widthOut			The width of the lattice (in the direction of latticeVector1).
	 * @param[out] heightOut		The height of the lattice (in the direction of latticeVector2).
	 */
	void calculateLatticeBoundary(Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &cornerOut, int &widthOut, int &heightOut);

	/*!
	 * Determines which of the points that the lattice was fit into are actually lying on the given lattice.
	 * For every such on-grid point, its global index (in terms of all model points) together with its grid position in the lattice
	 * (in terms of integer-valued (i,j)-coordinates, where position = latticeCorner + i*latticeVector1 + j*latticeVector2) is returned.
	 *
	 * @param[in] inputIndices	The global indices (in terms of all model points) of the points the lattice was fit into.
	 * @param[in] lattice		The lattice for wich the on grid points shall be determined. Note that the lattice needs to
	 * 							have been generated with this LatticeDetector instance to achieve the desired behavior.
	 * @return		A vector containing a pair of (global index, (i,j)-coordinates) for each on-grid point of the lattice.
	 * 				The global index of such a point is its index in terms of all model points. The (i,j)-coordinates are integer-valued,
	 * 				and determine the position of the corresponding lattice grid point, such that position = latticeCorner + i*latticeVector1 + j*latticeVector2.
	 */
	vector<pair<int,vector<int> > > getOnGridIndices(vector<int> const &inputIndices, LatticeStructure const &lattice);

	// GENERAL HELPER FUNCTION

	/*!
	 * Determines the similarity of two vectors.
	 * Case 1: 	The two vectors have similar length and direction.
	 * Case 2: 	The two vectors have similar length and similar but opposite direction (= similar orientation in space but with different sign).
	 * 			If this is the case, vector1 and (-vector2) are similar in the sense of Case 1 and the method will check for that.
	 * Case -1: The two vectors are not similar in terms of the former two cases.
	 * Similarity in the sense of Case 1 is defined in terms of a treshold, where two vectors are said to be similar if the norm of their
	 * difference vector is smaller than the treshold.
	 *
	 * @param[in] vector1	The first vector to compare.
	 * @param[in] vector2	The second vector to compare.
	 * @param[in] treshold	The treshold used to determine similarity.
	 * @return	The number of case that applies (1, 2 or -1). For details, see the method description.
	 */
	static int vectorsAreSimilar(Vector3d const &vector1, Vector3d const &vector2, double treshold);

	// CONSTANTS

	static constexpr double VECTOR_DISTANCE = 0.06;	/*!< Constant that is used
															a) as treshold for vectorsAreSimilar in this class, and
															b) as a treshold for lattice basis vector lengths, such that shorter basis vector candidates
																are discarded. 	*/

	static constexpr double TRESHOLD1 = 0.1; 		/*!< Constant that is used to determine whether a point is said to lie on a lattice grid point.
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	When determining whether a point is said to lie on a lattice grid point, the basis vector on which
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	the lattice grid point is approached is taken into consideration. So the on-grid check is
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	norm(point - lattice grid point) < TRESHOLD1 * norm(basisvector)  */

	static constexpr double TRESHOLD2 = 0.5; 		/*!< Constant that is used
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	a) to guide the lattice boundary expansion, by determining the ratio of valid points that needs
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 to be met in order to expand the lattice to another column / row
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	b) in the calculation of the valid-invalid ratio of points on a line, to determine the length of the
	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 	 line to be searched. */

	static constexpr double ANGLETRESHOLD = 0.0349;

private:

	vector<Vector3d> points; 	/*!< The points to fit a lattice. All points lie on a plane. */

	inputManager* inpManager; 	/*!< The input manager that gives access to all given initial data. */

	Eigen::Vector4d plane; 		/*!< The plane all the points lie on. */

	// *** GENERAL HELPER FUNCTIONS

	/*!
	 * Calculates the translation vector between two points
	 *
	 * @param[in] point1 The first point
	 * @param[in] point2 The second point
	 * @return	The translation vector from point2 to point1 (Care about that behavior!)
	 */
	Vector3d translationVector(Vector3d const &point1, Vector3d const &point2);


	// *** HELPER FUNCTIONS TO GET CANDIDATE VECTORS

	/*!
	 * Clusters candidate vectors together that are similar. For a vector to be added to a cluster, it is enough if the vector is similar to one of the vectors in the cluster.
	 *
	 * @param[in] candidates	Candidates to be clustered.
	 * @return	A list of clusters, each cluster itself represented by a list of candidates.
	 */
	list<list<Vector3d> > clusterCandidates(vector<Vector3d> const &candidates);

	/*!
	 * Combines clustered candidates by averaging them to a single candidate vector and saves the original number of vectors in the cluster as a score for the
	 * averaged candidate vector.  NOTE: This score is only important for the naive calculation of candidate vectors, where only the two most probable candidates
	 * are returned. The score mentioned in the paper, regarding the validation of candidates, is a different score that is calculated at a later stage, when
	 * searching for the final basis vectors.
	 *
	 * @param[in] 	clusteredCandidates	The clustered candidate vectors to be combined.
	 * @param[out] 	combinedCandidates	The combined candidate vectors. Each combined candidate vector is the average of one cluster of original candidate vectors.
	 * @param[out]	scores				The "scores" of the combined candidate vectors. The score of a combined candidate vector is the number of vectors in the
	 * 										cluster that gave rise to the combined candidate vector.
	 */
	void combineCandidates(list<list<Vector3d> > const &clusteredCandidates, vector<Vector3d> &combinedCandidates, vector<int> &scores);


	// *** HELPER FUNCTIONS TO VALIDATE CANDIDATE VECTORS

	/*!
	 * Validates the candidate basis vectors. The higher the score, the more probable it is that the vector is true basis vector.
	 * @param[in] candidateVectors	The candidate vectors.
	 * @return	The score for every candidate vector.
	 */
	vector<double> validateCandidateVectors(vector<Vector3d> const &candidateVectors);

	/*!
	 * Validates a candidate basis vector. The higher the score, the more probable it is that the vector is a true basis vector.
	 *
	 * @param[in]	candidateVector	The candidate basis vector.
	 * @return	The score for the candidate basis vector.
	 */
	double validateVector(Vector3d const &candidateVector);

	/*!
	 * Calculates the ratio of valid points to invalid points on the line through a reference point, spanned by a candidate basis vector.
	 *
	 * @param[in] referencePoint	The reference point.
	 * @param[in] candidateVector	The candidate basis vector.
	 * @return	The ratio of valid points to invalid points on the line.
	 */
	double validInvalidRatio(Vector3d const &referencePoint, Vector3d const &candidateVector);

	/*!
	 * Projects the points that the lattice was fit into onto the line through a reference point, in the direction of a candidate vector.
	 *
	 * @param[in] referencePoint	The reference point.
	 * @param[in] candidateVector	The candidate vector.
	 * @return	The projected points, lying on the line through the reference point in the direction of the candidate vector.
	 */
	vector<Vector3d> projectPointsOnLine(Vector3d const &referencePoint, Vector3d const &candidateVector);

	/*!
	 * Determines the outermost on-grid point indices (coordinates) for a given grid line; A grid line goes through a reference point and
	 * expands in space by adding integer multiples of the candidate vector to the reference point. Note that these integers can also be negative.
	 * Each point reachable in this way is called a "grid point", with its index/coordinate the integer that the candidate vector was multiplied
	 * with to reach it. One of the points that the lattice shall be fit into is called "on-grid" if it is close to a grid point. The method returns
	 * the outermost (minimum and maximum) indices of all grid points of the given grid line that have a close "on-grid" point.
	 *
	 * @param[in] projectedPoints	The projections of the points that the lattice shall be fit into onto the line through the reference point,
	 * 								in the direction of the candidate vector.
	 * @param[in] referencePoint	The reference point.
	 * @param[in] candidateVector	The candidate vector.
	 * @return	A vector containing the minimum and maximum on-grid point indices: [minIndex, maxIndex].
	 */
	vector<int> getOutermostOnGridPointIndices(vector<Vector3d> const &projectedPoints, Vector3d const &referencePoint, Vector3d const &candidateVector);


	/*!
	* Checks whether the pointToTest is a valid lattice point, i.e.
	* it's SIFT descriptor in the most frontoparallel view is 
	* similar to the referencePoint's SIFT.
	* calls the function computeSiftFronto, located in 3dtools.h.
	*
	* @param[in] referencePoint the 3D point from which the lattice expands
	* @param[in] pointToTest    the 3D point to check if it belongs to the lattice
	* return boolean variable according to whether the point is valid lattice point	
	*/
	bool isPointValid(Vector3d const &referencePoint, Vector3d const &pointToTest);

/*!
	* Checks whether the vector with index i is an integer combination of a subset of the vectors with index i+1:end.
	* @param[in] i index of the vector to test if it is an integer combination of the rest
	* @param[in] candidatesInOrder    the set of vectors to check the co-linearity	
	* @param[in] valid a validation index of each vector. Invalid the test vector will not be checked against invalid vectors
	* return boolean if the vector is an int.comb.
	*/
	bool isIntegerCombination(int i,vector<Vector3d>& candidatesInOrder,vector<bool>& valid);


	// *** HELPER FUNCTIONS TO GET LATTICE BOUNDARY

	/*!
	 * Calculates the lattice boundary in terms of corner, width and height, expanded around a certain reference point.
	 * The lattice is thought to span a parallelogram, originating at corner, with sides latticeVector1*width and latticeVector2*height
	 *
	 * @param[in] referencePoint	The reference point around which the lattice boundary is expanded.
	 * @param[in] latticeVector1	The first basis vector of the lattice.
	 * @param[in] latticeVector2	The second basis vector of the lattice.
	 * @param[out] cornerOut		The corner of the lattice.
	 * @param[out] widthOut			The width of the lattice (in the direction of latticeVector1).
	 * @param[out] heightOut		The height of the lattice (in the direction of latticeVector2).
	 */
	void latticeBoundaryForReferencePoint(Vector3d const &referencePoint, Vector3d const &latticeVector1, Vector3d const &latticeVector2, Vector3d &cornerOut, int &widthOut, int &heightOut);

	/*!
	 * Determines whether a certain grid line (originating at anchorPoint, spanning in the direction of directionVector, for [length] steps)
	 * is "valid", meaning that it contains a certain ratio of valid grid points. A grid point's validity is determined by comparing its appearance
	 * to a reference point.
	 *
	 * @param[in] referencePoint	The reference point. Grid points are compared with that reference point to determine whether they are valid or not.
	 * @param[in] anchorPoint		The anchor point, which is the origin / start point of the grid line in examination.
	 * @param[in] directionVector	The direction vector to build up the grid line originating at the anchor point. Grid points lie at
	 * 									integer multiples of the direction vector.
	 * @param[in] length			The length of the grid line in terms of direction vector multiples - the end of the grid is at anchorPoint+length*directionVector.
	 * @return	If true, the line is considered valid.
	 */
	bool validLine(Vector3d const &referencePoint, Vector3d const &anchorPoint, Vector3d const &directionVector, int length);

	/*!
	 * Determines whether a point is close enough to a certain grid point to be called "on-grid". The additional vector is needed to
	 * establish a proper treshold for the distance.
	 *
	 * @param[in] point		The point for which to decide whether it is "on-grid"
	 * @param[in] gridPoint	The grid point that the point should be checked against. If they are close enough, the point in question is called "on-grid".
	 * @param[in] vector	A vector to establish a proper distance treshold based on which the "on-grid" decision is made. See implementation for details.
	 * @return	If true, the point is called "on-grid". If false, the point is not considered to be "on-grid".
	 */
	bool pointEqualsGridPoint(Vector3d const &point, Vector3d const &gridPoint, Vector3d const &vector);

	// *** HELPER FUNCTIONS TO GET ON GRID POINTS

	/*!
	 * Changes the coordinate system for a set of points from canonical base to a base spanned by latticeVector1, latticeVector2 and latticeVector1xlatticeVector2
	 * @param[in] points			Points for which the coordinate system should be changed.
	 * @param[in] latticeVector1	Lattice vector 1.
	 * @param[in] latticeVector2	Lattice vector 2.
	 * @return	The points in the lattice basis coordinate system.
	 */
	vector<Vector3d> changeToLatticeBasis(vector<Vector3d> const &points, Vector3d const &latticeVector1, Vector3d const &latticeVector2);

};



#endif /* SRC_LATTICEDETECTOR_H_ */
