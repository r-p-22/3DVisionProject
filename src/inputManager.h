#ifndef INPUTMANAGER_H
#define INPUTMANAGER_H

#include <vector>
#include <Eigen/Dense>
#include <fstream>
#include <iostream>


#include "CImg.h"
#include "camera.h"
#include "latticeStruct.h" 

using namespace std;

/**
 * \class inputManager
 *
 *
 * This class provides an input interface and a container for the data available in a model. 
 * The data are loaded in the container, given the name of the files.
 * Contains getter and setter methods to provide an interface to the model data.
 * It should be created in the beginning, and imported to LatticeClass object.
 * 
 */

class inputManager{

	Eigen::Matrix<double,3,3> camK;
	vector<string> imageNames;
	vector<int> viewIds;

	vector<Eigen::Vector3d> allPoints;

	vector<Eigen::Matrix<double,3,4>> camPoses;


public:

	vector<TriangulatedPoint> pointModel; /*!< The struct that contains all the information related to a 3d point: position, cameras observed, 2d positions in these cameras. */



	/*!
	 * Method to update the allPoints field with the changed pointModel.
	 * The pointModel is provided to other methods, by reference, so it is possible to have its 
	 * values changed (e.g. in Bundle adjustment). This method should be called after BA.
	 */
	void updatePointsWithModel(){
		for (int i = 0; i < pointModel.size(); i++){
			allPoints[i] = pointModel[i].pos;
			//cout << allPoints[i] << endl;
			//cout << pointModel[i].pos << endl;
			//cout << "----"<<endl;
		}

	}

	vector<string> getImgNames(){
		return this->imageNames;
	}

	Eigen::Matrix<double,3,3> getK(){
			return this->camK;
	}

	vector<Eigen::Matrix<double,3,4>> getCamPoses(){
			return this->camPoses;
	}

	void setCamPoses(vector<Eigen::Matrix<double,3,4>> newcams ){
			this->camPoses = vector<Eigen::Matrix<double,3,4>>(newcams);
	}

	vector<int> getViewIds(){
		return this->viewIds;
	}
	vector<Eigen::Vector3d> getPoints(){
		return this->allPoints;
	}

	vector<TriangulatedPoint> getPointModel(){
		return this->pointModel;
	}

	inputManager(char** argv){
		read3Dpoints(argv[2],allPoints,pointModel);
		readCameras(argv[3],argv[4],camPoses,camK, viewIds);
		readImgNames(argv[1],imageNames);
	}

private:

	void  readImgNames(char* file, vector<string>& imageNames){

		ifstream is(file);
		if (!is){
					std::cout << "Images file not opened" << endl;
				}
		string name;
		while (is >> name)
		{
			imageNames.push_back(name);
		}
		is.close();

	}

	void read3Dpoints(char* file, vector<Eigen::Vector3d> &points3d,
			vector<TriangulatedPoint> &pointModel){

		std::ifstream instream(file);


		if (!instream){
			std::cout << "file not opened" << endl;
		}
		int nPoints;
		instream >> nPoints;

		for (int j = 0; j < nPoints; ++j) {
			TriangulatedPoint X;

			instream >> X.pos[0] >> X.pos[1] >> X.pos[2];
			int nMeasurements = 0;
			instream >> nMeasurements;
			for (int k = 0; k < nMeasurements; ++k) {
				PointMeasurement m;
				instream >> m.view >> m.id >> m.pos[0] >> m.pos[1];
				X.measurements.push_back(m);
			}
			pointModel.push_back(X);
			points3d.push_back(X.pos);
		}

		instream.close();

		return;
	}

	void readCameras(char* file, char* fileK, vector<Eigen::Matrix<double,3,4> >& cameraPoses,
			Eigen::Matrix<double,3,3> &K, vector<int>& viewIds){

		std::ifstream is(file);
		vector<TriangulatedPoint> pointModel;
		vector<Eigen::Vector3d> points3d;

		if (!is){
			std::cout << "file not opened" << endl;
		}


		int nViews;
		is >> nViews;
		std::cout << "Going to read " << nViews << " poses." << endl;

		Eigen::Matrix<double,3,4> P;

		for (int i = 0; i < nViews; ++i)
		{
			int viewId;
			is >> viewId;
			viewIds.push_back(viewId);

			is >> P(0,0) >> P(0,1) >> P(0,2) >> P(0,3);
			is >> P(1,0) >> P(1,1) >> P(1,2) >> P(1,3);
			is >> P(2,0) >> P(2,1) >> P(2,2) >> P(2,3);
			cameraPoses.push_back(P);
		}

		K = Eigen::Matrix3d::Identity(3,3);
		std::ifstream is2(fileK);
		is2 >> K(0,0) >> K(0,1) >> K(0,2) >> K(1,1) >> K(1,2);

		is2.close();
	}

};


#endif

