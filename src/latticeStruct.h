
#include <Eigen/Dense>
#include <vector>
#include <math.h>
struct LatticeStructure
{
	Eigen::Vector4d plane;
	std::vector<Eigen::Vector3d> basisVectors;
	std::vector<Eigen::Vector3d> boundary;

};

int computeNumberOfCells(LatticeStructure latt){
	Vector3d LL = latt.boundary[0];
	Vector3d TR = latt.boundary[1];
	Vector3d basis1 = latt.basisVectors[0];
	Vector3d basis2 = latt.basisVectors[1];

	Matrix<double,3,3> A;
	Vector3d solution;
	A.block<3,1>(0,2) = -(TR-LL);
	A.block<3,1>(0,0) = basis1;
	A.block<3,1>(0,1) = basis2;
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeFullV);
	solution = svd.matrixV().block<3,1>(0,2);//<sizeRows,sizeCols>(beginRow,beginCol)
	solution = solution/solution[2];

	int k1 = round(solution[0]);
	int k2 = round(solution[1]);

	return k1*k2;
}


//Consolidate/merge the initial lattices, to produce a final lattice list
vector<LatticeStructure> consolidateLattices(vector<LatticeStructure> inputLattices){

	vector<LatticeStructure> finalLattices;

	finalLattices.push_back(inputLattices[0]);

	bool matched;

	vector<LatticeStructure>::iterator initialLatticeIterator;
	for (initialLatticeIterator = inputLattices.begin()+1; initialLatticeIterator != inputLattices.end(); initialLatticeIterator++){

		LatticeStructure latt = *initialLatticeIterator;
		matched = false;

		vector<LatticeStructure>::iterator finalLattIterator;
		for (int i = 0; i < finalLattices.size(); i++){

			LatticeStructure lattF = finalLattices[i];

			//if translation vector is less that a threshold, then merge
			double basisVecThresh = sqrt((lattF.basisVectors[0] - lattF.basisVectors[1]).squaredNorm());
			double costheta = latt.plane.dot(lattF.plane)/(latt.plane.norm()*lattF.plane.norm());
			if ( (abs(latt.plane[3] - lattF.plane[3]) < basisVecThresh*0.1 ) && (acos(costheta) <= 0.0349 ) )  {
				matched = true;

				int ncells = computeNumberOfCells(latt);
				int ncellsF = computeNumberOfCells(lattF);

				//the final (merged) lattice will be the one with most cells
				if (ncells > ncellsF){
					finalLattices[i] = latt;
				}

				break;
			}

		}
		if (!matched){
			finalLattices.push_back(latt);
		}


	}

	return finalLattices;
}

