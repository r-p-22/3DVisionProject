#ifndef LATTICECLASS_H
#define LATTICECLASS_H


#include "CImg.h"
#include "camera.h"

#include "planeFitter.h"
#include "latticeDetector.h"
#include "latticeStruct.h"
#include "my_v3d_vrmlio.h" // already imported in main_test2

#include <iomanip>

using namespace std;

class LatticeClass {

	LatticeDetector* LattDetector;
	PlaneFitter pf;

	  bool computeDensifyingPoints = false;
	        vector<TriangulatedPoint> densifyingPoints;
	        vector<pair<int, vector<int> > > densifyingLatticeGridIndices;
	        string densifyingPointsFile = "data/densifyingPoints.txt";
	        string densifyingPointsIndicesFile = "data/densifyingPointsIndices.txt";


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

    	// TESTED - OK: planeInlierIdx aligned with planeInliersProjected


    	LattStructure.plane = pf.getFittedPlane();
//    	cout << "plane:" << endl;
//    	cout << LattStructure.plane << endl;

    	//-----Fit lattice----------------

    		//.0 initialize the class
    	LattDetector = new LatticeDetector(planeInliersProjected,LattStructure.plane,inpM);
      		//.1 calculate candidate basis vectors
    	vector<Eigen::Vector3d> candidateBasisVecs = LattDetector->calculateCandidateVectors(0);

    		//.2 calculate final basis vectors
    	LattStructure.basisVectors = LattDetector->getFinalBasisVectors(candidateBasisVecs);

//    		cout << "calculated final bvecs " << endl;
//	    	cout << candidateBasisVecs[0] << endl;
//	    	cout << candidateBasisVecs[1] << endl;

    		//.3 calculate boundaries
		LattDetector->calculateLatticeBoundary(LattStructure.basisVectors[0], LattStructure.basisVectors[1], LattStructure.lowerLeftCorner, LattStructure.width, LattStructure.height);

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
        this->latticeGridIndices = LattDetector->getOnGridIndices(planeInlierIdx, LattStructure);


        //LattStructure = this->inpM->getPointModel()

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

	float calculateReprojectionError(){
		Vector2f pa;
		CameraMatrix cam;
		cam.setIntrinsic(inpM->getK());

		float reprojectionError = 0;
		for (int p=0; p < this->groupPointsIdx.size(); p++){
			int pointidx  = groupPointsIdx[p];

			for (int j=0; j<inpM->getPointModel()[pointidx].measurements.size();j++){
				int imgview = inpM->getPointModel()[pointidx].measurements[j].view;

				int i=0;
				for (i=0;i<inpM->getCamPoses().size();i++){
					if (inpM->getViewIds()[i] == imgview)
						break;
				}
				Eigen::Matrix<double,3,4> P = inpM->getCamPoses()[i];
				pa = cam.projectPoint(inpM->getPointModel()[pointidx].pos).cast<float>();

				reprojectionError += (pa-inpM->getPointModel()[pointidx].measurements[j].pos).squaredNorm();
			}
		}

		return reprojectionError;
	}

	void projectLatticeToImage(bool debug = false){

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
		Vector2d pa2d, pb2d;

		for (size_t k=0; k < this->pointsInGroup.size(); k++){
			pa2d = cam.projectPoint(inpM->getPointModel()[groupPointsIdx[k]].pos);//.cast<float>();
			//pa2d = cam.projectPoint(group[k]).cast<float>();
			image.draw_circle(float(pa2d[0]),float(pa2d[1]),4,color,1);

		}

		Vector3d pa; pa = latt.lowerLeftCorner;
		Vector3d pb; pb = B1;
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

			if (debug){
				if (k==0){
					int c[] = {255,0,0};
					image.draw_circle(pa2d[0],pa2d[1],6,c,1);
				}
			}
			pa += basis1;
			pb += basis1;
		}

		if (debug){
			for (int j=0; j< this->latticeGridIndices.size();j++){
				pa2d = cam.projectPoint(inpM->getPointModel()[latticeGridIndices[j].first].pos);
				image.draw_circle(float(pa2d[0]),float(pa2d[1]),2+1*j,color,1);
			}
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
		for (size_t kk = 0; kk <inpM->getPointModel()[pointidx].measurements.size(); kk+=3 ){

			//get the view id and the respected image name for this point
			//The m.view refer to the image index
			// Pose index points to the image index

			int imgview = inpM->getPointModel()[pointidx].measurements[kk].view;

			string img = inpM->getImgNames()[imgview];

			//get the camera pose for this viewid
			size_t i = 0;
			for (i = 0; i < inpM->getCamPoses().size(); i++){
				if (inpM->getViewIds()[i] == imgview)
					break;
			}

			P = inpM->getCamPoses()[i];

			cimg_library::CImg<unsigned char> image(("data/"+img).c_str());

			cam.setOrientation(P);

			Vector2f pa2d;
			for (size_t k=0; k < group.size(); k++){
				pa2d = cam.projectPoint(inpM->getPointModel()[groupPointsIdx[k]].pos).cast<float>();
				//pa2d = cam.projectPoint(group[k]).cast<float>();
				image.draw_circle(pa2d[0],pa2d[1],5,color,1);
				if (k==0)
					cout << pa2d.cast<float>() - inpM->getPointModel()[pointidx].measurements[kk].pos << endl;
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


	int densifyStructure(int start_index){

		int new_point_count = 0;            // counter variable to handle indexing

		cout << "Adding Lattice grid points as new structure 3d points:" << endl;

		// toggle console output off
		streambuf *old = cout.rdbuf(0);

		if(computeDensifyingPoints)
		{

			// loop through lattice structure and look for missing points
			for(size_t i = 0; i<=LattStructure.width; i+=1)               // gridpoint horizontal coordinate
			{
				for(size_t j = 0; j<=LattStructure.height; j+=1)          // gridpoint vertical coordinates
				{
					bool point_exists = false;

					cout << "Checking lattice grid point (" << i << "," << j << ")" << " : " << endl;
					for(size_t k = 0; k<latticeGridIndices.size(); k+=1)
					{
						// check if a point already exists for the given grid point
						if(latticeGridIndices[k].second[0] == i && latticeGridIndices[k].second[1] == j)
						{
							cout << "3d point already exists." << endl;
							point_exists = true;
						}
					}
					if(!point_exists)
					{
						cout << "adding new 3d point with views "<< endl;

						Eigen::Vector3d pos;
						pos = LattStructure.lowerLeftCorner
								+i*LattStructure.basisVectors[0]
								+j*LattStructure.basisVectors[1];

						// add most frontoparallel view as only measurement
						CameraMatrix cam;
						cam.setIntrinsic(inpM->getK());

						vector<PointMeasurement> ms;
						double cosangle = 0;
						double tmpcosangle=0;
						int bestview=-1;
						Vector4d plane = LattStructure.plane;
						Vector2d pbest(0,0);
						Vector2d p;

						for(size_t i=0; i<inpM->getCamPoses().size(); i++)
						{
							//get view
							int view = inpM->getViewIds()[i];

							if ((view < 45) || (view > 47))
							{
								continue;
							}

							string img = inpM->getImgNames()[view];
							cimg_library::CImg<unsigned char> image(("data/"+img).c_str());
							float const w = image.width();
							float const h = image.height();

							//check angle between camera-point line and plane normal
							Vector3d line = pos - inpM->getCamPoses()[i].block<3,1>(0,3);
							//abs because we dont know the plane orientation
							tmpcosangle = abs(line.dot(plane.head(3)))/sqrt(line.squaredNorm() * plane.head(3).squaredNorm());

							//project point into image
							cam.setOrientation(inpM->getCamPoses()[i]);
							p = cam.projectPoint(pos);

							double d = cam.transformPointIntoCameraSpace(pos)[2];

							if ((d > 0) && (tmpcosangle > cosangle) && (p[0]>=0)&&(p[1]>=0)&& (p[0]<w)&&(p[1]<h)){
								cosangle = tmpcosangle;
								bestview = view;
								pbest = p;
								//if ( abs(cosangle) > 0.9)
								//	break;
							   }
						}

						// bestview ist the one we add as measurement
						string img = inpM->getImgNames()[bestview];
						cimg_library::CImg<unsigned char> image(("data/"+img).c_str());

						/*
						int bestViewIndex;
						for(int l = 0; l<inpM->getViewIds().size(); l++)
						{
							if(inpM->getViewIds()[l] == bestview)
								bestViewIndex = l;
						}

						cam.setOrientation(inpM->getCamPoses()[bestViewIndex]);
						p = cam.projectPoint(pos);


						const unsigned char color[] = { 0,0,255 };
						image.draw_circle(p[0],p[1],5,color,1);

						cimg_library::CImgDisplay main_disp(image,"");

						while (!main_disp.is_closed()){
							main_disp.wait(500);
						}
						*/

						Eigen::Vector2f p2f;
						p2f << p(0), p(1);
						PointMeasurement newMeasurement(p2f,bestview);
						ms.push_back(newMeasurement);

						// add new point to densifying points container
						TriangulatedPoint newPoint(pos,ms);
						densifyingPoints.push_back(newPoint);

						// update lattice grid indices
						vector<int> gridPosition;
						pair<int,vector<int> > newPointPair;

						gridPosition.push_back(i);
						gridPosition.push_back(j);
						newPointPair.first = start_index+new_point_count;
						newPointPair.second = gridPosition;

						densifyingLatticeGridIndices.push_back(newPointPair);
						latticeGridIndices.push_back(newPointPair);

						new_point_count++;

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

			// save updated lattGridPoint indices
			ofstream os2(densifyingPointsIndicesFile);
			if(!os2.good())
			{
				cout << "Problem opening filestream for saving densifying points indices" << endl;
				return -1;
			}

			//output total number of added points

			// save total number of grid points
			os2 << densifyingLatticeGridIndices.size() << endl;

			for (size_t i = 0; i<densifyingLatticeGridIndices.size(); i++)
			{
				// output index of 3d point corresponding to grid point
				os2 << densifyingLatticeGridIndices[i].first << " ";

				// output grid coordinates
				os2 << densifyingLatticeGridIndices[i].second[0] << " ";
				os2 << densifyingLatticeGridIndices[i].second[1] << endl;
			}
			os2.close();

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

			int n_points, n_measurements;
			double pos3d_x, pos3d_y, pos3d_z, view, id, pos_x, pos_y, index;
			is >> n_points;

			for(size_t i = 0; i<n_points; i+=1)
			{
				//read 3d point coordinates
				is >> pos3d_x >> pos3d_y >> pos3d_z;
				Eigen::Vector3d pos;
				pos << pos3d_x, pos3d_y, pos3d_z;

				// read measurements
				is >> n_measurements;
				vector<PointMeasurement> ms;

				for(size_t m = 0; m<n_measurements; m+=1)
				{
					is >> view >> id >> pos_x >> pos_y;
					Eigen::Vector2f pos2d;
					pos2d << pos_x,pos_y;

					PointMeasurement newMeasurement(pos2d,view);
					ms.push_back(newMeasurement);
				}

				// add new point to densifying points vector
				TriangulatedPoint newPoint(pos,ms);
				densifyingPoints.push_back(newPoint);
				cout << "Read densifying point" << pos.transpose() << endl;

			} // end reading points

			// read extra lattGridIndices
			ifstream is2(densifyingPointsIndicesFile);
			if(!is2.good())
			{
				cout << "Problem opening file to read densifying points indices from file" << endl;
				return -1;
			}

			is2 >> n_points;

			for(size_t i = 0; i<n_points; i+=1)
			{
				//read information to fill lattGridIndices
				is2 >> index >> pos_x >> pos_y;

				vector<int> gridPosition;
				pair<int,vector<int> > newPointPair;

				gridPosition.push_back(pos_x);
				gridPosition.push_back(pos_y);
				newPointPair.first = start_index+new_point_count;
				newPointPair.second = gridPosition;

				new_point_count++;

				densifyingLatticeGridIndices.push_back(newPointPair);
				latticeGridIndices.push_back(newPointPair);

				cout << "Read latticeGridIndices ( " << pos_x <<  "," << pos_y << ") with 3dpoint index" << newPointPair.first << endl;
			}
			is2.close();

		} // end case read from file

		// toggle cout on
		cout.rdbuf(old);

		// print information to consoles
		if(new_point_count == densifyingPoints.size())
			cout << "Added " << new_point_count << " new points." << endl;
		else
			cerr << "Error counting new points from densification process." << endl;

		// update end index: end_index = last index used+1.
		int end_index = start_index+new_point_count;

		return end_index;       // returns end index of added points in pointModel container.
	}



};



#endif
