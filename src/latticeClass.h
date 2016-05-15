#ifndef LATTICECLASS_H
#define LATTICECLASS_H


#include "CImg.h"
#include "camera.h"

#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"
#include "my_v3d_vrmlio.h" // already imported in main_test2



using namespace std;

class LatticeClass {

	LatticeDetector LattDetector;
	PlaneFitter pf;

public:

	// the 3d points in the group and their indices
	vector<Vector3d> pointsInGroup;
	vector<int> groupPointsIdx;

	// the inlier 3d points in the group
	vector<Vector3d> planeInliersProjected;
	vector<int> planeInlierIdx;

	LatticeStructure LattStructure;

	vector<pair<int, vector<int> > > latticeGridIndices;

	inputManager* inpM;

	/*constructor: input:
		the inputManager (i.e. image names, cameras, points, etc.)
		the group points (3D points)
		the respected point indices
	*/
	LatticeClass(inputManager& inpm, vector<Vector3d>& _groupPoints, vector<int>& _groupPointsIndices){
		this->inpM = &inpm;
		this->pointsInGroup = _groupPoints;
		this->groupPointsIdx = _groupPointsIndices;

	}

	/* 2nd constructor: load the LatticeStructure and latticeGridIndices from file
	 *  */
	LatticeClass(inputManager& inpm, vector<Vector3d> _groupPoints, vector<int> _groupPointsIndices,
			const char* file){
		this->inpM = &inpm;
		this->pointsInGroup  = _groupPoints;
		this->groupPointsIdx = _groupPointsIndices;
		loadFromFile(file);
	}
	~LatticeClass(){
		inpM=NULL;

	}

	// copy constructor
	LatticeClass(const LatticeClass& cSource) {
		inpM = cSource.inpM;

		groupPointsIdx = vector<int>(cSource.groupPointsIdx);
		pointsInGroup = vector<Vector3d>(cSource.pointsInGroup);

		LattStructure = cSource.LattStructure;

		planeInliersProjected = vector<Vector3d>(cSource.planeInliersProjected);
		planeInlierIdx = vector<int>(cSource.planeInlierIdx);

		latticeGridIndices = vector<pair<int, vector<int> > >(cSource.latticeGridIndices);

	}

	//Method to compute the lattice end-to-end,
	//It will populate the fields of LatticeStructure.
	void fitLattice(){

		//-----Fit the plane--------------
		this->planeInlierIdx = pf.ransacFit(pointsInGroup,groupPointsIdx);
    	planeInliersProjected = pf.getProjectedInliers();

    	LattStructure.plane = pf.getFittedPlane();
//    	cout << "plane:" << endl;
//    	cout << LattStructure.plane << endl;

    	//-----Fit lattice----------------

    		//.0 initialize the class
    	LattDetector.reconstructedPoints = planeInliersProjected;
    	LattDetector.plane = LattStructure.plane;
    	LattDetector.inpManager = inpM;

    		//.1 calculate candidate basis vectors
    	vector<Eigen::Vector3d> candidateBasisVecs = LattDetector.calculateCandidateVectors(0);

    		//.2 calculate final basis vectors
    	LattStructure.basisVectors = LattDetector.getFinalBasisVectors(candidateBasisVecs);

//    		cout << "calculated final bvecs " << endl;
//	    	cout << candidateBasisVecs[0] << endl;
//	    	cout << candidateBasisVecs[1] << endl;

    		//.3 calculate boundaries
		LattDetector.calculateLatticeBoundary(LattStructure.basisVectors[0], LattStructure.basisVectors[1], LattStructure.lowerLeftCorner, LattStructure.width, LattStructure.height);

			cout << "plane: " << endl;
			cout << LattStructure.plane << endl;
			cout << "calculated final bvecs: " << endl;
			cout << LattStructure.basisVectors[0] << endl;
			cout << "---- " << endl;
			cout << LattStructure.basisVectors[1] << endl;
			cout << "---- " << endl;
			cout << "boundary computed. " << endl;
			cout << LattStructure.lowerLeftCorner << endl;
			cout << "width: " << LattStructure.width << ". height: " << LattStructure.height << "." << endl;
			cout << "---- " << endl;

		//----get the indices of the on-grid points
        this->latticeGridIndices = LattDetector.getOnGridIndices(planeInlierIdx, LattStructure);


	}


	void saveLatticeToFile(const char* file){

			ofstream os;

			os.open(file,ios::out);

			os << LattStructure.basisVectors[0].x() << endl;
			os << LattStructure.basisVectors[0].y() << endl;
			os << LattStructure.basisVectors[0].z() << endl;

			os << LattStructure.basisVectors[1].x() << endl;
			os << LattStructure.basisVectors[1].y() << endl;
			os << LattStructure.basisVectors[1].z() << endl;

			os << LattStructure.width << endl;
			os << LattStructure.height << endl;

			os << LattStructure.lowerLeftCorner.x() << endl;
			os << LattStructure.lowerLeftCorner.y() << endl;
			os << LattStructure.lowerLeftCorner.z() << endl;

			os << LattStructure.plane[0] << endl;
			os << LattStructure.plane[1] << endl;
			os << LattStructure.plane[2] << endl;
			os << LattStructure.plane[3] << endl;

			int indicesCount = latticeGridIndices.size();

			os << indicesCount << endl;

			for (int i = 0; i < indicesCount; i++){
				pair<int,vector<int> > gridIndex = latticeGridIndices[i];
				os << gridIndex.first << endl;
				os << gridIndex.second[0] << endl;
				os << gridIndex.second[1] << endl;
			}

			os.close();
		}

		void loadFromFile(const char* file){

			ifstream is;

			is.open(file);

			this->LattStructure = LatticeStructure();
			this->LattStructure.basisVectors = vector<Vector3d>();

			double x,y,z,w;

			is >> x;
			is >> y;
			is >> z;
			this->LattStructure.basisVectors.push_back(Vector3d(x,y,z));

			is >> x;
			is >> y;
			is >> z;
			this->LattStructure.basisVectors.push_back(Vector3d(x,y,z));

			is >> LattStructure.width;
			is >> LattStructure.height;

			is >> x;
			is >> y;
			is >> z;
			this->LattStructure.lowerLeftCorner = Vector3d(x,y,z);

			is >> x;
			is >> y;
			is >> z;
			is >> w;
			this->LattStructure.plane = Vector4d(x,y,z,w);

			this->latticeGridIndices = vector<pair<int,vector<int> > >();

			int indicesCount;

			is >> indicesCount;

			for (int i = 0; i < indicesCount; i++){
				pair<int,vector<int> > gridIndex = pair<int, vector<int> >();
				is >> gridIndex.first;
				gridIndex.second = vector<int>();
				int width, height;
				is >> width;
				is >> height;
				gridIndex.second.push_back(width);
				gridIndex.second.push_back(height);

				this->latticeGridIndices.push_back(gridIndex);
			}

			is.close();
		}

	void projectLatticeToImage(){

		LatticeStructure latt = this->LattStructure;
		Vector3d basis1 = latt.basisVectors[0];
		Vector3d basis2 = latt.basisVectors[1];

		int k1 = latt.width;
		int k2 = latt.height;

		Vector3d B1 = latt.lowerLeftCorner + k1*basis1;
		Vector3d B2 = latt.lowerLeftCorner + k2*basis2;

		CameraMatrix cam;
		cam.setIntrinsic(inpM->getK());


		//selects the 1st image that the 1st point is visible
		int pointidx  = latticeGridIndices[0].first;
		int imgview = inpM->getPointModel()[pointidx].measurements[0].view;
		string img = inpM->getImgNames()[imgview];
		int i=0;
		for (i=0;i<inpM->getCamPoses().size();i++){
			if (inpM->getViewIds()[i] == imgview)
				break;
		}

		Eigen::Matrix<double,3,4> P = inpM->getCamPoses()[i];
		//float const w = 1696;
		//float const h = 1132;
		cimg_library::CImg<unsigned char> image(("data/"+img).c_str());
		const unsigned char color[] = { 0,0,255 };

		cam.setOrientation(P);
		Vector3d pa; pa = latt.lowerLeftCorner;
		Vector3d pb; pb = B1;
		Vector2d pa2d, pb2d;
		for (int k=0; k<=k2; k++){
			//project pa, pb into image
			pa2d = cam.projectPoint(pa);
			pb2d = cam.projectPoint(pb);
			//plot 2D line into image
			image.draw_line(pa2d[0],pa2d[1],pb2d[0],pb2d[1],color);

			pa += basis2;
			pb += basis2;

		}
		pa = latt.lowerLeftCorner;
		pb = B2;
		for (int k=0; k<=k1; k++){
			//project pa, pb into image
			pa2d = cam.projectPoint(pa);
			pb2d = cam.projectPoint(pb);
			//plot 2D line into image
			image.draw_line(pa2d[0],pa2d[1],pb2d[0],pb2d[1],color);

			pa += basis1;
			pb += basis1;
		}

		cimg_library::CImgDisplay main_disp(image,"");

		//image.save("latt_view45.png");
		while (!main_disp.is_closed()){
			main_disp.wait();
		}
	}

	void projectGroupToImage(){

		vector<Vector3d> group = this->pointsInGroup;

		CameraMatrix cam;
		cam.setIntrinsic(inpM->getK());
		const unsigned char color[] = { 0,0,255 };
		Eigen::Matrix<double,3,4> P;

		int pointidx  = groupPointsIdx[0];
		cout << pointidx << endl;
		for (size_t kk = 0; kk <inpM->getPointModel()[pointidx].measurements.size(); kk+=3 ){

			//get the view id and the respected image name for this point
			//The m.view refer to the image index
			// Pose index points to the image index

			int imgview = inpM->getPointModel()[pointidx].measurements[kk].view;
			cout << imgview << endl;

			string img = inpM->getImgNames()[imgview];

			//get the camera pose for this viewid
			size_t i = 0;
			for (i = 0; i < inpM->getCamPoses().size(); i++){
				if (inpM->getViewIds()[i] == imgview)
					break;
			}
			cout << i << endl;

			P = inpM->getCamPoses()[i];

			cimg_library::CImg<unsigned char> image(("data/"+img).c_str());

			cam.setOrientation(P);
			cout << "gsize: "<<group.size() << endl;

			Vector2f pa2d;
			for (size_t k=0; k < group.size(); k++){
				pa2d = cam.projectPoint(inpM->getPointModel()[groupPointsIdx[k]].pos).cast<float>();
				//pa2d = cam.projectPoint(group[k]).cast<float>();
				image.draw_circle(pa2d[0],pa2d[1],5,color,1);

				if (k == 0){
					cout << " projection pos diff:" << endl;
					cout << (pa2d - inpM->getPointModel()[pointidx].measurements[kk].pos) << endl;
				}

			}

			cimg_library::CImgDisplay main_disp(image,"");

			//image.save("latt_view45.png");
			while (!main_disp.is_closed()){
				main_disp.wait(500);
			}

		}
	}


	void writeToVRML(const char* filename, const bool append = true){
		writeLatticeToVRML(this->LattStructure.plane,this->LattStructure.basisVectors,
				this->LattStructure.lowerLeftCorner,this->LattStructure.width,this->LattStructure.height,
				filename, append);
	}

};



#endif
