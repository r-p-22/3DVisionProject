#ifndef LATTICECLASS_H
#define LATTICECLASS_H


#include "CImg.h"
#include "camera.h"

#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"

//#include "my_v3d_vrmlio.h" // already imported in main_test2



using namespace std;

class LatticeClass {

	LatticeDetector LattDetector;
	PlaneFitter pf;

	LatticeStructure LattStructure;


	inputManager* inpM;

	vector<Vector3d> pointsInGroup;
	vector<int> groupPointsIdx;

	vector<Vector3d> planeInliersProjected;
	vector<int> planeInlierIdx;

	vector<pair<int, vector<int> > > latticeGridIndices;

public:


	/*constructor: input:
		the inputManager (i.e. image names, cameras, points, etc.)
		the group points (3D points)
		the respected point indices */
	LatticeClass(inputManager& inpm, const vector<Vector3d>& _groupPoints, const vector<int>& _groupPointsIndices){
		this->inpM = &inpm;
		this->pointsInGroup = _groupPoints;
		this->groupPointsIdx = _groupPointsIndices;

	}

	/* 2nd constructor: load the LatticeStructure and latticeGridIndices from file */
	LatticeClass(inputManager& inpm, char* file){
			this->inpM = &inpm;
			loadFromFile(file);
	}

	~LatticeClass(){}

	//Method to compute the lattice end-to-end,
	//It will populate the fields of LatticeStructure.
	void fitLattice(){

		//-----Fit the plane--------------
		this->planeInlierIdx = pf.ransacFit(pointsInGroup,groupPointsIdx);
    	planeInliersProjected = pf.getProjectedInliers();

    	LattStructure.plane = pf.getFittedPlane();

    	//-----Fit lattice----------------

    		//.0 initialize the class
    	LattDetector.reconstructedPoints = planeInliersProjected;
    	LattDetector.plane = LattStructure.plane;
    	LattDetector.inpManager = inpM;

    		//.1 calculate candidate basis vectors
    	vector<Eigen::Vector3d> candidateBasisVecs = LattDetector.calculateCandidateVectors(0);

    		//.2 calculate final basis vectors
    	LattStructure.basisVectors = LattDetector.getFinalBasisVectors(candidateBasisVecs);

    		cout << "calculated final bvecs " << endl;

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


	void saveLatticeToFile(char* file){

	}

	void loadFromFile(char* file){

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


		//selects the 5th image that the 1st point is visible
		int pointidx  = latticeGridIndices[0].first;
		int view = inpM->getPointModel()[pointidx].measurements[5].view;
		string img = inpM->getImgNames()[view];
		int i=0;
		for (i=0;i<inpM->getCamPoses().size();i++){
			if (inpM->getViewIds()[i] == view)
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

		//select the 1st point of the group is visible
		int pointidx  = groupPointsIdx[0];

		for (size_t kk = 0; kk <inpM->getPointModel()[pointidx].measurements.size(); kk+=1 ){

			//get the view id and the respected image name for this point
			//TODO: Does the m.view refer to the camPose index?
			int poseview = inpM->getPointModel()[pointidx].measurements[kk].view;

			string img = inpM->getImgNames()[inpM->getViewIds()[poseview]];

			//get the camera pose for this viewid
			int i = poseview;
			/*for (i = 0; i < inpM->getCamPoses().size(); i++){
				if (inpM->getViewIds()[i] == poseview)
					break;
			}*/

			Eigen::Matrix<double,3,4> P = inpM->getCamPoses()[i];

			cimg_library::CImg<unsigned char> image(("data/"+img).c_str());


			cam.setOrientation(P);

			Vector2d pa2d;
			for (size_t k=0; k < group.size(); k++){
				pa2d = cam.projectPoint(group[k]);
				image.draw_circle(pa2d[0],pa2d[1],5,color,1);

				if (k == 0){
					cout << " projected point pos:" << endl;
					cout << pa2d << endl;
				}

			}

			cimg_library::CImgDisplay main_disp(image,"");

			cout << "view: ";
			cout << poseview << endl;
			cout << "pose: ";
			cout << i << endl;
			cout << "registered point pos:" << endl;
			cout << inpM->getPointModel()[pointidx].measurements[kk].pos << endl;

			//image.save("latt_view45.png");
			while (!main_disp.is_closed()){
				main_disp.wait(500);
			}

		}
	}


	void writeToVRML(const char* filename, const bool append = true){
		writeLatticeToVRML(this->LattStructure.plane,this->LattStructure.basisVectors,
				this->LattStructure.lowerLeftCorner,this->LattStructure.width,this->LattStructure.height,
				filename,append);
	}

};



#endif