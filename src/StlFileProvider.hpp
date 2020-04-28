/*
 * StlFileProvider.hpp
 *
 *  Created on: Nov 15, 2012
 *      Author: petr
 */
#pragma once
#ifndef STLFILEPROVIDER_HPP_
#define STLFILEPROVIDER_HPP_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <limits>
#include <Eigen/Dense>
// #include <boost/iostreams/device/back_inserter.hpp>
// #include <boost/iostreams/stream.hpp>
// #include <boost/archive/binary_oarchive.hpp>
// #include <boost/archive/binary_iarchive.hpp>
// #include <boost/serialization/vector.hpp>
#include "TriangleElement.hpp"
#include "BoundingBox.hpp"
#include "FileProvider.hpp"


/// Loads the data from the STL file and creates geometry class which stores the elements inside.
template <typename GeometryManager, typename type>
class StlFileProvider : public FileProvider {



public:
	StlFileProvider(){};
	StlFileProvider(const char* file);

	virtual void readFile(const char* iname) {};

	friend std::ostream& operator<<(std::ostream& os, const StlFileProvider& provider);

};


///====================================
///          Implementation
///====================================


StlFileProvider::StlFileProvider(const char* file) {

	int rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	//readInputFile(file);

}


void StlFileProvider::readFile(const char* iname) {

	std::ifstream iFile;
	char          line[120];
	char          name[20];

	iFile.open(iname);

	// read first line and the name of the body from it
	if (!iFile.is_open()) {
		// file not existent or some other stuff
		// throw error? :)
	}

	// read first line and the name of the body from it
	iFile.getline(line, 120);

	sscanf(line, "%*s %s", name);

	long int nElements     = 0;
	long int nNormals      = 0;
	long int nVertices     = 0;
	int      nFacet        = 0;
	int      nOuter        = 0;
	int      vertexCounter = 0;

	/*aabb.maxX[0] = -std::numeric_limits<double>::max();
	aabb.maxX[1] = -std::numeric_limits<double>::max();
	aabb.maxX[2] = -std::numeric_limits<double>::max();

	aabb.minX[0] = std::numeric_limits<double>::max();
	aabb.minX[1] = std::numeric_limits<double>::max();
	aabb.minX[2] = std::numeric_limits<double>::max();
*/

	TriangleElement<double> workingElement;
	elements.reserve(11000000);

	std::vector<Eigen::Vector3d > vertices;


	// go through the rest of the lines and load the data
	while ( iFile.getline(line, 120)) {

		char flag[120];

		sscanf(line, "%s", flag);

		if ( strstr(line, "facet normal") ) {
			// load the normal, we do not use them, but lets read them
			double normal[3];
			sscanf(line, "%*s %*s %lf %lf %lf", &normal[0], &normal[1], &normal[2]);

			workingElement.ID = nElements; // simply element number

			nFacet++;
			nElements++;
			nNormals++;
		}
		if ( strstr(flag, "outer") ) {
			nOuter++;
			vertexCounter = 0;
		}
		if ( strstr(flag, "vertex")  != 0 ) {
			Eigen::Vector3d vertex;

			sscanf(line, "%*s %lf %lf %lf", &vertex[0],
										    &vertex[1],
										    &vertex[2]);

			if (aabb.maxX[0] < vertex[0])
				aabb.maxX[0] = vertex[0];
			if (aabb.maxX[1] < vertex[1])
				aabb.maxX[1] = vertex[1];
			if (aabb.maxX[2] < vertex[2])
				aabb.maxX[2] = vertex[2];

			if (aabb.minX[0] > vertex[0])
				aabb.minX[0] = vertex[0];
			if (aabb.minX[1] > vertex[1])
				aabb.minX[1] = vertex[1];
			if (aabb.minX[2] > vertex[2])
				aabb.minX[2] = vertex[2];

			vertices.push_back(vertex);

			nVertices++;
			vertexCounter++;
		}
		if ( strstr(flag, "endloop")  != 0 ) {
			nOuter--;
		}
		if ( strstr(line, "endfacet")) {

			// compute dx, dy, dz
			double dx;

			TriangleElement<double> tr(vertices[0], vertices[1], vertices[2], nElements);

			elements.push_back( tr );

//			std::cout << tr << std::endl;

			nFacet--;
			vertices.clear();
		}

		// read next line
		line[0] = 'X';

	}

	if ( (nFacet != 0) && (nOuter != 0) && ((nVertices/3) != nElements) ) {
		// file was loaded incorectly, some shit happened
		std::cout << "Error" << std::endl;
	}

	numberOfElements = nElements;
	dim              = 3;

//	std::cout << aabb << std::endl;

	iFile.close();

	// done loading elemenents

}



// ostream operators


std::ostream& operator<<(std::ostream& os, const StlFileProvider& provider) {

	os << "Dimension : " << provider.dim << std::endl;
	os << "Number of elements: " << provider.numberOfElements << std::endl;
	os << "Min DX : " << provider.mindx[0] << ", " << provider.mindx[1] << ", " << provider.mindx[2] << std::endl;
	os << "AABB : " << provider.aabb << std::endl;

    return os;
}

#endif /* STLFILEPROVIDER_HPP_ */
