#pragma once
#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include <iostream>
#include <fstream>
#include <vector>
#include <list>
#include <map>
#include <limits>
#include <spatial/point_multiset.hpp>
#include <spatial/metric.hpp>
#include <spatial/bits/spatial_region.hpp>
#include <spatial/bits/spatial_neighbor.hpp>
#include <spatial/bits/spatial_overlap_region.hpp>
#include "mpi.h"
#include "utility.h"
#include "TriangleElement.hpp"
#include "BoundingBox.hpp"
#include "SearchGrid.hpp"



struct cDistance{
	double signedDistance;
	double pDistance;

	cDistance() : signedDistance(std::numeric_limits<double>::max()), pDistance(std::numeric_limits<double>::min()) {};
	cDistance(double sd, double pd) : signedDistance(sd), pDistance(pd) {};

};

class SphereGeometry {
private:

	typedef BoundingBox<double, 3, TriangleElement<double> > BB_type;
	typedef spatial::point_multiset<3, TriangleElement<double>, TriangleCompare> triagTree_type;
//	typedef spatial::box_multiset<6, BB_type, BoundingBoxCompare> boxTree_type;
	typedef std::vector< TriangleElement<double> > elementsList_type;

	/// The dimension of the mesh 
	int dim;

	/// The number of elements
	long int numberOfElements;

	int _numberOfCellsPerTask;

	elementsList_type elements;

	bool searchInit;
	triagTree_type triagIndex;

	friend std::ostream& operator<<(std::ostream& os, const Geometry& geom);


public:

	Box<double, 3>   aabb; // global aabb
	std::vector<int> localElements;
	std::vector<int> groupElements;

	Geometry() : dim(3),
                 numberOfElements(0),
                 _numberOfCellsPerTask(0),
                 searchInit(false) {}


	template <typename iterator>
	Geometry(iterator it_begin, iterator it_end, const BB_type& span) : dim(span.dims()),
																		aabb(span),
																		_numberOfCellsPerTask(0),
																		searchInit(false) {
		this->elements.assign(it_begin, it_end);
		numberOfElements = elements.size();
	}

	Geometry(const Geometry& geom): dim(geom.dim),
			                        numberOfElements(geom.numberOfElements),
			                        aabb(geom.aabb),
			                        _numberOfCellsPerTask(geom._numberOfCellsPerTask),
			                        searchInit(false) {
		elements = geom.elements;
	}

	~Geometry() {
		if ( searchInit ) {

		}
	}


	void appendElement(const TriangleElement<double> & element);

	void preSortElements(const Box<double, 3>& localRegion, double dx, int numberOfSubGroups, int groupID);

	int getNumberOfCellsPerDomainLocal() const;


	void initSearchAccelerator();


	double getMaxX(int i);
	double getMinX(int i);
	long int getNumberOfElements() const {
		return numberOfElements;
	};

	const TriangleElement<double>& getElement(int ind) const {
		return elements[ind];
	}



	const elementsList_type::const_iterator getElementsIterator() const;
	const elementsList_type::const_iterator getElementsEndIterator() const;

	double computeDistance(const Eigen::Vector3d& point, bool acc);
	void computeDistance(const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points,
							   cDistance distances[]);


	// Factory method to load geometry from STL file
	static Geometry fromFile(const char* iname, double enlargeBoundingBox);

};



inline void Geometry::appendElement(const TriangleElement<double> & element) {

	elements.push_back(element);
	numberOfElements++;

	// expand the bounding box - will be used as root node
	aabb.expand( element.getBoundingBox() );

}

/// This method is designed to provide labeling elements by colour of their group.
/// This means that you can split elements into coloured groups and use them in paralellization.
/// This method also assing indices to local part of geometry for local geometry processing.
/// Local geometry processing will be superseded by kd tree division.
void Geometry::preSortElements(const Box<double, 3>& localRegion, double dx, int numberOfSubGroups, int groupID) {

	typedef spatial::region_iterator<triagTree_type, TrianglePredicate> regionIterator_type;

	Box<double, 3> lb = localRegion;
//	lb.grow(1.2);

	regionIterator_type localIter = spatial::region_begin(triagIndex, TrianglePredicate(lb));
	regionIterator_type endIter = spatial::region_end(triagIndex, TrianglePredicate(lb));

	int elementGroupWidth = ceil( double(elements.size()) / double(numberOfSubGroups) );


	for (; localIter != endIter; ++localIter) {
		localElements.push_back( (*localIter).ID );
	}


//	_numberOfCellsPerTask = 0;
//	// this can be done efficiently using tree search.
//	for (size_t i = 0; i < elements.size(); ++i) {
////		BoundingBox<double, 3, TriangleElement<double> > elementBox = elements[i].getBoundingBox();
//
////		if ( localRegion.intersect( elementBox ) ) {
////			localElements.push_back(i);
////			_numberOfCellsPerTask += ((elementBox.tr() - elementBox.bl()) / dx).prod();
////		}
//
//		if ( (i/elementGroupWidth) == groupID ) {
//			groupElements.push_back(i);
//		}

//	}

	std::cout << "local elements : " << localElements.size() << std::endl;

}

int Geometry::getNumberOfCellsPerDomainLocal() const {
	return _numberOfCellsPerTask;
}


// Initialize hierarchichal spatial search
// int his case it is KD-TREE implemented in FLANN library
void Geometry::initSearchAccelerator() {
	for (auto& element: elements) {
		triagIndex.insert( element );
		std::cout << "inserted" << std::endl;
	}

	searchInit = true;

}


inline double Geometry::getMaxX(int i) {
	return aabb.maxX(i);
}

inline double Geometry::getMinX(int i) {
	return aabb.minX(i);
}

inline const std::vector< TriangleElement<double> >::const_iterator Geometry::getElementsIterator() const {
	return elements.begin();
}

inline const std::vector< TriangleElement<double> >::const_iterator Geometry::getElementsEndIterator() const {
	return elements.end();
}

double Geometry::computeDistance(const Eigen::Vector3d& point, bool acc) {
	double eps           =  my_eps;
	double perpDistance  = -std::numeric_limits<double>::max();

	SignedDistance<double> distance;
	distance.dist = std::numeric_limits<double>::max();

	if ( searchInit && acc ) {

		// std::cout << "accelerated search" << std::endl;

		spatial::neighbor_iterator<triagTree_type, TriangleMetric> knnIt =
						spatial::neighbor_begin(triagIndex, TriangleMetric(), TriangleElement<double>(point));

		long int id = -1;
		for (int k = 0; k < 10; ++k, ++knnIt) {

			// std::cout << (*knnIt).ID << std::endl;
			// std::cout << (*knnIt) << std::endl;

			// for each nearest element compute distance
			SignedDistance<double> currentDistance =
					TriangleElement<double>::computeDistance((*knnIt), point);
			double currentPerpDistance =
					TriangleElement<double>::computePerpendicularDistance((*knnIt), point);

			// std::cout << knnIt.distance() << std::endl;

			if ( std::abs( currentDistance.dist - distance.dist ) < eps  ) {
				if ( perpDistance <= currentPerpDistance ) {
					distance     = currentDistance;
					perpDistance = currentPerpDistance;
					id = (*knnIt).ID;
				}
			} else if ( currentDistance.dist < distance.dist ) {
				distance     = currentDistance;
				perpDistance = currentPerpDistance;
				id = (*knnIt).ID;
			}
		}
		// std::cout << "++++++++++++++++++++" << std::endl;
		// std::cout << "ACC ID = " << id << std::endl;
		return distance.dist * distance.sign;

	} else {
		long int id = -1;
		for (auto element: elements) {
			// for all elements in the geometry compute closest match
			SignedDistance<double> currentDistance = TriangleElement<double>::computeDistance(element, point);
			double             currentPerpDistance = TriangleElement<double>::computePerpendicularDistance(element, point);

			if ( std::abs( currentDistance.dist - distance.dist ) < eps  ) {
				if ( perpDistance <= currentPerpDistance ) {
					distance     = currentDistance;
					perpDistance = currentPerpDistance;
					id = element.ID;
				}
			} else if ( currentDistance.dist < distance.dist ) {
				distance     = currentDistance;
				perpDistance = currentPerpDistance;
				id = element.ID;
			}
		}
		// std::cout << "ID = " << id << std::endl;
		// std::cout << elements[id].getBoundingBox() << std::endl;
		return distance.dist * distance.sign;

	}
}


void Geometry::computeDistance(
		const std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >& points,
		cDistance distances[]) {

	// used to compute minimal distance to surface

	double eps           =  my_eps;
	double perpDistance  = -std::numeric_limits<double>::max();

	SignedDistance<double> distance;
	distance.dist = std::numeric_limits<double>::max();

	std::cout << "compute distance" << std::endl;

	if ( searchInit ) {

		int i = 0;
		std::cout << "accelerated search" << std::endl;
		for (auto point: points) {

			spatial::neighbor_iterator<triagTree_type, TriangleMetric> knnIt =
							spatial::neighbor_begin(triagIndex, TriangleMetric(), TriangleElement<double>(point));

			for (int k = 0; k < 10; ++k, ++knnIt) {

					// for each nearest element compute distance
				SignedDistance<double> currentDistance =
						TriangleElement<double>::computeDistance((*knnIt), points[i]);
				double currentPerpDistance =
						TriangleElement<double>::computePerpendicularDistance((*knnIt), points[i]);

				if ( std::abs( currentDistance.dist - distance.dist ) < eps  ) {
					if ( perpDistance <= currentPerpDistance ) {
						distance     = currentDistance;
						perpDistance = currentPerpDistance;
					}
				} else if ( currentDistance.dist < distance.dist ) {
					distance     = currentDistance;
					perpDistance = currentPerpDistance;
				}
			}

			distances[i].signedDistance = distance.dist * distance.sign;
			distances[i].pDistance = perpDistance;

			i++;

		}

	} else {

		for (size_t i = 0; i < points.size(); ++i) {
			// for each point loop and compute accorging distance

			for (auto element: elements) {
				// for all elements in the geometry compute closest match
				SignedDistance<double> currentDistance = TriangleElement<double>::computeDistance(element, points[i]);
				double             currentPerpDistance = TriangleElement<double>::computePerpendicularDistance(element, points[i]);

				if ( std::abs( currentDistance.dist - distance.dist ) < eps  ) {
					if ( perpDistance <= currentPerpDistance ) {
						distance     = currentDistance;
						perpDistance = currentPerpDistance;
					}
				} else if ( currentDistance.dist < distance.dist ) {
					distance     = currentDistance;
					perpDistance = currentPerpDistance;
				}
			}

			distances[i].signedDistance = distance.dist * distance.sign;
			distances[i].pDistance = perpDistance;

		}
	}

}








////////////////////////////////////////////////////////////////////////////////
//                            FACTORY METHOD
////////////////////////////////////////////////////////////////////////////////






/// factory method to load data from stl file
Geometry Geometry::fromFile(const char* iname, double enlargeBoundingBox) {


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


	std::vector<TriangleElement<double> > loadedElements;
	TriangleElement<double>               workingElement;
	Box<double, 3>                        aabb;

	typedef std::vector<TriangleElement<double> >::iterator elit;

	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > vertices;
	Eigen::Vector3d                                                          vertex, normal;

	Geometry geom;

	// go through the rest of the lines and load the data
	while ( iFile.getline(line, 120) ) {

		char flag[120];

		sscanf(line, "%s", flag);

		if ( strstr(line, "facet normal") ) {
			// load the normal, we do not use them, but lets read them
			sscanf(line, "%*s %*s %lf %lf %lf", &normal[0], &normal[1], &normal[2]);

			nFacet++;
			nElements++;
			nNormals++;
		}
		if ( strstr(flag, "outer") ) {
			nOuter++;
			vertexCounter = 0;
		}
		if ( strstr(flag, "vertex")  != 0 ) {

			sscanf(line, "%*s %lf %lf %lf", &vertex[0],
										    &vertex[1],
										    &vertex[2]);

			vertices.push_back(vertex);

			nVertices++;
			vertexCounter++;
		}
		if ( strstr(flag, "endloop")  != 0 ) {
			nOuter--;
		}
		if ( strstr(line, "endfacet")) {
		
			geom.appendElement( TriangleElement<double>(vertices, (nElements - 1)) );

			nFacet--;
			vertices.clear();
		}

		// read next line
		line[0] = 'X';

	}

	if ( (nFacet != 0) && (nOuter != 0) && ((nVertices/3) != nElements) ) {
		// file was loaded incorectlly, some shit happened
		std::cout << "Error" << std::endl;
	}

	iFile.close();
	// done loading elements

	// this part of code is used to split elements into coloured groups
//	int rank, nprocs;
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
//
//	int elements_per_proc = loadedElements.size() / nprocs;
//
//	elit begin = loadedElements.begin() + ( rank * elements_per_proc);
//	elit end   = (rank < (nprocs-1)) ? loadedElements.begin() + ( (rank+1) * elements_per_proc) : loadedElements.end();
//	Geometry geom(begin, end, aabb);


	geom.aabb.grow(enlargeBoundingBox);

	return geom;

};




//==========================================
// ostream operators used for debug printing
//==========================================


std::ostream& operator<<(std::ostream& os, const Geometry& geom) {
	
	os << "Geometry : " << std::endl;
	os << "Bounding box : " << std::endl;
	os << geom.aabb << std::endl;

	os << "--------------" <<std::endl;

	os << "Number of elements : " << geom.numberOfElements << std::endl;

    return os;
}

#endif


