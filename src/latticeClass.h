#ifndef LATTICECLASS_H
#define LATTICECLASS_H


#include "CImg.h"
#include "camera.h"

#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"

#include <iomanip>


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

    bool computeDensifyingPoints = true;
    vector<TriangulatedPoint> densifyingPoints;
    string densifyingPointsFile = "data/densifyingPoints.txt";


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

	void loadFromFile(char* file){

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


    int densifyStructure(){

        if(computeDensifyingPoints)
        {
            cout << "Adding Lattice grid points as new structure 3d points." << endl;
            // loop through lattice structure and look for missing points
            for(size_t i = 0; i<LattStructure.width; i+=1)               // gridpoint horizontal coordinate
            {
                for(size_t j = 0; j<LattStructure.height; j+=1)          // gridpoint vertical coordinates
                {
                    bool point_exists = false;
                    for(size_t k = 0; k<latticeGridIndices.size(); k+=1)
                    {
                        // check if a point already exists for the given grid point
                        if(latticeGridIndices[k].second[0] == i && latticeGridIndices[k].second[1] == j)
                            point_exists = true;
                    }
                    if(!point_exists)
                    {
                        // add new 3d point
                        Eigen::Vector3d pos;
                        pos = LattStructure.lowerLeftCorner
                                +i*LattStructure.basisVectors[0]
                                +j*LattStructure.basisVectors[1];

                        // add views of new 3d point
                        CameraMatrix cam;
                        cam.setIntrinsic(inpM->getK());

                        vector<PointMeasurement> ms;
                        for(size_t view = 0; view < inpM->getViewIds().size(); view+=1)
                        {
                            string img = inpM->getImgNames()[inpM->getViewIds()[view]];
                            cimg_library::CImg<unsigned char> image(("data/"+img).c_str());

                            Eigen::Matrix<double,3,4> P = inpM->getCamPoses()[view];

                            cam.setOrientation(P);
                            Vector2d pa2d = cam.projectPoint(pos);

                            // check if 3d point projects onto image -> add new measurement
                            if(pa2d[0] < image.width() && pa2d[1] < image.height())
                            {
                                Eigen::Vector2f pos2d;
                                pos2d << pa2d[0], pa2d[1];

                                PointMeasurement newMeasurement(pos2d,view);
                                ms.push_back(newMeasurement);
                            }

                        } // end checking views to add

                        // add new point to densifying points container
                        TriangulatedPoint newPoint(pos,ms);
                        densifyingPoints.push_back(newPoint);
                        // TODO: add to inpM->pointModel as well?
                        // inpM->getPointModel().push_back(newPoint);
                    }
                } // end vertical search of grid points
            } // end horizontal search of grid points

            // save densifying points to separate file
            ofstream os(densifyingPointsFile);
            if(!os.good())
            {
                cout << "Problem opening filestream for saving densifying points" << endl;
                return -1;
            }

            // save total number of densifying points
            os << densifyingPoints.size() << endl;

            for(size_t i = 0; i<densifyingPoints.size(); i+=1)
            {
                // output 3d coordinates
                os << setprecision(10)
                   << densifyingPoints[i].pos[0] << " " << densifyingPoints[i].pos[1] << " " << densifyingPoints[i].pos[2] << " ";

                // output measurements
                os << setprecision(10)
                   << densifyingPoints[i].measurements.size() << " ";

                for(size_t m = 0; m<densifyingPoints[i].measurements.size(); m+=1)
                {
                    os << setprecision(10)
                       << densifyingPoints[i].measurements[m].view << " "
                       << densifyingPoints[i].measurements[m].id << " "
                       << densifyingPoints[i].measurements[m].pos[0] << " "
                       << densifyingPoints[i].measurements[m].pos[1];
                }
                os << endl;
            }
            os.close();

        } // end recompute densifying points

        else
        {
            // read densifying points from file
            ifstream is(densifyingPointsFile);
            if(!is.good())
            {
                cout << "Problem opening file to read densifying points from" << endl;
                return -1;
            }

            int n_points;
            is >> n_points;

            for(size_t i = 0; i<n_points; i+=1)
            {
                //read 3d point coordinates
                double pos3d_x, pos3d_y, pos3d_z;
                is >> pos3d_x >> pos3d_y >> pos3d_z;
                Eigen::Vector3d pos;
                pos << pos3d_x, pos3d_y, pos3d_z;

                // read measurements
                int n_measurements;
                is >> n_measurements;
                vector<PointMeasurement> ms;

                for(size_t m = 0; m<n_measurements; i+=1)
                {
                    double view, id, pos_x, pos_y;
                    is >> view >> id >> pos_x >> pos_y;
                    Eigen::Vector2f pos2d;
                    pos2d << pos_x,pos_y;
                    PointMeasurement newMeasurement(pos2d,view);
                    ms.push_back(newMeasurement);
                }

                // add new point to densifying points vector
                TriangulatedPoint newPoint(pos,ms);
                densifyingPoints.push_back(newPoint);
                // TODO: add to inpM->pointModel as well?
                // inpM->getPointModel().push_back(newPoint);

            } // end reading points
        } // end case read from file

        return 0;
    } // end densifying points method
};



#endif
