/*
 * spatialTest.cpp
 *
 *  Created on: Feb 5, 2014
 *      Author: petr
 */

#include <iostream>
#include <Eigen/Dense>

#include "Geometry.hpp"

int main (int argc, char* argv[]) {

	Geometry geom = Geometry::fromFile("/home/petr/Dropbox/testGeometries/dragon.stl", 1.0);
	geom.initSearchAccelerator();


	std::cout << Eigen::Vector3d({ 0, 0, 0 }).transpose() << " = " << geom.computeDistance(Eigen::Vector3d({ 0, 0, 0 }), true) << std::endl;
	std::cout << Eigen::Vector3d({ 0, 0, 0 }).transpose() << " = " << geom.computeDistance(Eigen::Vector3d({ 0, 0, 0 }), false) << std::endl;

	return 0;
}





