/*
 * SphereProvider.hpp
 *
 *  Created on: Nov 10, 2012
 *      Author: petr
 */
#pragma once
#ifndef SPHEREPROVIDER_HPP_
#define SPHEREPROVIDER_HPP_

#include <vector>


class SphereProvider {

	// [ [x,y,z], [x,y,z], [x,y,z], ... ]
	std::vector<std::vector<double> > coords;

	std::vector<double> center;

	double radius;

	int dim;

	int numberOfPoints;

	std::vector<double> minX;

	std::vector<double> maxX;

public:

	SphereProvider() {};
	SphereProvider(std::vector<double> center, double radius, int nPoints);

	void initializeGeometry();
	void initializeBoundingBox();


	double computeDistance(const std::vector<double>& point);

	const double getMinX(int dim) const;
	const double getMaxX(int dim) const;

	const std::vector<std::vector<double> > & getCoord() const;



};


///===================================
///          Implementation
///===================================


SphereProvider::SphereProvider(std::vector<double> center, double radius, int nPoints) {
	this->dim            = center.size();
	this->center         = center;
	this->radius         = radius;
	this->numberOfPoints = nPoints;

	initializeGeometry();
	initializeBoundingBox();

}


void SphereProvider::initializeGeometry() {
	// creates the sphere points in Nd

	if (dim == 1) {

		// for 1D case only two points
		std::vector<double> coordRow;
		double x = center[0] - radius;
		coordRow.push_back(x);
		coords.push_back(coordRow);

		x = center[0] + radius;
		coordRow.push_back(x);
		coords.push_back(coordRow);

	} else if (dim == 2) {

		// x(t) = center(1) - radius*cos(angle_step * t);
		// y(t) = center(2) - radius*sin(angle_step * t);

		double angle_step = (2*3.1415926535897) / numberOfPoints;

		for (int i = 0; i < numberOfPoints; ++i) {

			std::vector<double> coordRow;
			double x = center[0] - radius * cos(angle_step * i);
			double y = center[1] - radius * sin(angle_step * i);

			coordRow.push_back(x);
			coordRow.push_back(y);


			coords.push_back(coordRow);
		}

	} else if (dim == 3) {

//		x(ind) = center(1) - radius * sin(angleStep1(i))*cos(angleStep2(j));
//		y(ind) = center(2) - radius * sin(angleStep1(i))*sin(angleStep2(j));
//		z(ind) = center(3) - radius * cos(angleStep1(i));

		double angle_step_1 = (2*3.1415926535897) /  numberOfPoints;
		double angle_step_2 = (3.1415926535897)   / (numberOfPoints / 2);


		for (int i = 0; i <= numberOfPoints; ++i) {

			for (int j = 0; j <= numberOfPoints / 2; ++j) {
				std::vector<double> coordRow;

				coordRow.push_back( center[0] - radius*sin(angle_step_1 * i)*cos(angle_step_2*j) );
				coordRow.push_back( center[1] - radius*sin(angle_step_1 * i)*sin(angle_step_2*j) );
				coordRow.push_back( center[2] - radius*cos(angle_step_1 * i) );

				coords.push_back(coordRow);
			}

		}

	}
}

void SphereProvider::initializeBoundingBox() {

	if (dim >= 1) {
		minX.push_back( center[0] - radius );
		maxX.push_back( center[0] + radius );
	}

	if (dim >= 2) {
		minX.push_back( center[1] - radius );
		maxX.push_back( center[1] + radius );
	}

	if (dim == 3) {
		minX.push_back( center[2] - radius );
		maxX.push_back( center[2] + radius );
	}

}

double SphereProvider::computeDistance(const std::vector<double>& point) {
	// this function really speaks for itself. It is just sphere equation

	double dist = 0;

	switch (dim) {
		case 1:
			dist = sqrt( (point[0] - center[0])*(point[0] - center[0]) ) - radius;
			break;
		case 2:
			dist = sqrt( (point[0] - center[0])*(point[0] - center[0]) +
					     (point[1] - center[1])*(point[1] - center[1])) - radius;
			break;
		case 3:
			dist = sqrt( (point[0] - center[0])*(point[0] - center[0]) +
				         (point[1] - center[1])*(point[1] - center[1]) +
				         (point[2] - center[2])*(point[2] - center[2])) - radius;
			break;
		default:
			break;
	}

	return dist;

}

inline const double SphereProvider::getMinX(int dim) const{
	return minX[dim];
}

inline const double SphereProvider::getMaxX(int dim) const{
	return maxX[dim];
}

inline const std::vector<std::vector<double> > & SphereProvider::getCoord() const {
	return coords;
}



#endif /* SPHEREPROVIDER_HPP_ */
