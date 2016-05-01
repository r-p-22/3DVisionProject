/*
 * PlaneFitter.cpp
 *
 *  Created on: Mar 10, 2016
 *      Author: lianos91
 */

#include "planeFitter.h"

#include <iostream>

using namespace std;

PlaneFitter::PlaneFitter() {
	this->fittedplane = new Eigen::Vector4d();
	this->Dinliers = new Eigen::Matrix<double,4,Eigen::Dynamic>(4,1); //dummy initialization
}

PlaneFitter::~PlaneFitter() {
	delete this->fittedplane;
}

Eigen::Vector4d PlaneFitter::getFittedPlane(){
	return *this->fittedplane;
}

Eigen::Matrix<double,4,Eigen::Dynamic> PlaneFitter::getInlierPoints(){
	return *this->Dinliers;
}

std::vector<int> PlaneFitter::ransacFit(std::vector<Eigen::Vector3d> points3d, std::vector<int> inputIndices){

	Eigen::Matrix<double,4,Eigen::Dynamic> D(4,points3d.size());
	this->vecToEigenMat(points3d,D);

	int N = D.cols();
	cout << N << endl;

	// RANSAC variables
	float p = 0.99;
	float testRansiter = 0;
	int bestInlNum = 0;
	int best_it=-1;
	int inl_num;
	srand (time(NULL));
	Eigen::Vector4d *testplane = new Eigen::Vector4d();


	//keep the indices of the final inliers
	std::vector<int> bestInlierIds;

	// RANSAC loop
	for (int iter=0;iter<1000;iter++){

		//sample 3 different random points
		std::vector<int> samplevec;
		Eigen::Matrix<double,3,4> D_ransac_sample;
		int x;
		for (int s=0; s<3; s++){
			do{
				x = rand() % N;
			}while (std::find(samplevec.begin(), samplevec.end(), x) != samplevec.end());
			samplevec.push_back(x);
			D_ransac_sample.block<1,4>(s,0) = D.block<4,1>(0,x).transpose();//<sizeRows,sizeCols>(beginRow,beginCol)
		}

		//fit model
		this->fitPlane(D_ransac_sample,testplane);

		std::vector<int> inl_id;
		//inlier test
		inl_num = 0;
		for (int i=0;i<N;i++){
			if ( abs(D.block<4,1>(0,i).dot(*testplane)/testplane->head(3).norm()) < this->RANSAC_THRESH ){
				inl_num++;
				inl_id.push_back(i);
			}
		}

		//cout << "iter,inliers: " << iter <<","<<inl_num<< endl;

		//test (adaptive ransac)
		testRansiter = log(1-p)/log(1-pow((float)inl_num/N,3));
		if (inl_num > bestInlNum){
			bestInlNum = inl_num;
			best_it = iter;
			bestInlierIds.swap(inl_id);
		}

	}

	//bestplane contains the fitted plane
	cout << "best iter,inliers: " << best_it <<","<<bestInlNum<< endl;


	//Refine best plane
	delete this->Dinliers;
	this->Dinliers = new Eigen::Matrix<double,4,Eigen::Dynamic>(4,bestInlNum);

	this->vecToEigenMat(points3d,*this->Dinliers,bestInlierIds);

	//get the indices of the inlier points for the output
	vector<int> outputIndices;
	for (size_t i=0; i< bestInlierIds.size(); i++){
		outputIndices.push_back(inputIndices[bestInlierIds[i]]);
	}

	this->fitPlane(Dinliers->transpose(),testplane);

	delete this->fittedplane;
	this->fittedplane = new Eigen::Vector4d();
	this->fittedplane = testplane;
	return outputIndices;
}


void PlaneFitter::fitPlane(Eigen::Matrix<double,Eigen::Dynamic,4> Points,Eigen::Vector4d *plane){

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(Points, Eigen::ComputeThinU | Eigen::ComputeFullV);
	*plane = svd.matrixV().block<4,1>(0,3);//<sizeRows,sizeCols>(beginRow,beginCol)

	return;
}

std::vector<Eigen::Vector3d> PlaneFitter::getProjectedInliers(){
		vector<Eigen::Vector3d> latticePoints;
		Eigen::Matrix<double,4,Eigen::Dynamic> Dbest = *this->Dinliers;
		Eigen::Vector4d bestplane = *this->fittedplane;
		for (int i=0; i<Dbest.cols(); i++){
			Eigen::Vector3d p = Dbest.block<3,1>(0,i);
			float dist = p.dot(bestplane.head<3>()) + bestplane(3);
			//cout << dist << endl;

			p = p - bestplane.head<3>()*(dist)/sqrt(bestplane.head<3>().squaredNorm());
			latticePoints.push_back(p);
		}
		return latticePoints;

	}

