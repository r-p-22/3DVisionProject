/*
 * PlaneFitter.h
 *
 *  Created on: Mar 10, 2016
 *      Author: lianos91
 */

#ifndef PLANEFITTER_H_
#define PLANEFITTER_H_

#include <Eigen/Dense>
#include <vector>

class PlaneFitter {

public:

	PlaneFitter();

	virtual ~PlaneFitter();

	float RANSAC_THRESH = 0.01;

	Eigen::Vector4d* fittedplane;

	Eigen::Matrix<double,4,Eigen::Dynamic>* Dinliers;

	Eigen::Vector4d getFittedPlane();
	Eigen::Matrix<double,4,Eigen::Dynamic> getInlierPoints();

	void fitPlane(Eigen::Matrix<double,Eigen::Dynamic,4> Points, Eigen::Vector4d *plane);

	void ransacFit(std::vector<Eigen::Vector3d>);

	std::vector<Eigen::Vector3d> getProjectedInliers();

	static void vecToEigenMat(std::vector<Eigen::Vector3d> in,
			Eigen::Matrix<double,4,Eigen::Dynamic> &out, std::vector<int> idvec = std::vector<int>()){

		if (idvec.size()==0){
			for (size_t i = 0;i < in.size();i++){
				out.block<3,1>(0,i) = in.at(i);
				out(3,i) = 1.;
			}
		}else{
			for (size_t i = 0;i < idvec.size();i++){
				out.block<3,1>(0,i) = in.at(idvec[i]);
				out(3,i) = 1.;
			}
		}
		return;
	}


};

#endif /* PLANEFITTER_H_ */
