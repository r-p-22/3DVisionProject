
#include <Eigen/Dense>
#include <vector>

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

	double cos_phi1 = (TR-LL).dot(basis1)/( sqrt((TR-LL).squaredNorm())*sqrt(basis1.squaredNorm()));
	Vector3d B1; B1 = sqrt((TR-LL).squaredNorm())*cos_phi1*basis1/sqrt(basis1.squaredNorm()) + LL;

	double cos_phi2 = (TR-LL).dot(basis2)/( sqrt((TR-LL).squaredNorm())*sqrt(basis2.squaredNorm()));
			Vector3d B2; B2 = sqrt((TR-LL).squaredNorm())*cos_phi2*basis2/sqrt(basis2.squaredNorm()) + LL;

	//number of lattices in the 2 axes
	int k1 = round( sqrt((B1-LL).squaredNorm())/sqrt(basis1.squaredNorm()) );
	int k2 = round( sqrt((B2-LL).squaredNorm())/sqrt(basis2.squaredNorm()) );

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
			if ( abs(latt.plane[3] - lattF.plane[3]) < basisVecThresh*0.1 ) {
				matched = true;

				int ncells = computeNumberOfCells(latt);
				int ncellsF = computeNumberOfCells(lattF);

				//the final (merged) lattice will be the one with more cells
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

