/*
 * Initializer.hpp
 *
 *  Created on: Nov 8, 2012
 *      Author: petr
 */
#pragma once
#ifndef INITIALIZER_NEW_HPP_
#define INITIALIZER_NEW_HPP_

#include <vector>
#include <limits>
#include <petscdmadda.h>
#include "TriangleElement.hpp"
#include "BoundingBox.hpp"


template <typename DataProvider, int dim>
class Initializer {

	// Actually all data should come from DataProvider class, even the distance function

	DataProvider* dataProvider;


	// utility function, tranform coordinates to indices
	Box<int, 3> locateIndices(const Box<double, 3>& coord, Grid<double, dim>* gr);

	// select which boundaries need to be initialized
	int selectBoundaries(Grid<double, dim>* gr);

	// check if node ows the whole geometry
	bool haveAllGeometry(Grid<double, dim>* gr);

	// check if have some geometry
	bool haveSomeGeometry(Grid<double, dim>* gr);

	// check if coordinate is on this node
	bool haveElement(TriangleElement<double>& el, Grid<double, dim>* gr) const;

	// check if coordinate is on this node
	bool haveCoordinate(const TriangleElement<double>& coord, Grid<double, dim>* gr) const;


public:

	long int cellCounter;

	Initializer();
	Initializer(DataProvider* provider): dataProvider(provider) {};
	~Initializer(){};

	// this method will do the magic and initialize the local part of the data accordingly to the
	// defined input stuff
	void updateArray(PetscReal* data, Grid<double, dim>* gr);

};


///================================================
///               Implementation
///================================================


template <typename DataProvider, int dim> void Initializer<DataProvider, dim>::updateArray(PetscReal* data, Grid<double, dim>* gr) {

	double eps        = 1E-05;
	double narrowBand = ((gr->getDx(0) + gr->getDx(1) + gr->getDx(2)) / 3) * 10;

	cellCounter = 0;

	DMDALocalInfo* info = gr->getLocalInfo();

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// perpendicular distance computed towards the triangle, it is used when we are no longer sure, which normal shall we use
	long int gridCellNum = gr->getNumberOfLocalCells();

	double* perpendicualDistance = new double[gridCellNum];

	for (long int i = 0; i < gridCellNum; ++i) {
		data[i] = std::numeric_limits<double>::max();
		perpendicualDistance[i] = -std::numeric_limits<double>::max();
	}

	// place starting geometry
	if ( haveSomeGeometry(gr) ) {

		long int ackCounter = 0;
		long int rejCounter = 0;

		for (auto it = dataProvider->getElementsIterator(); it != dataProvider->getElementsEndIterator(); ++it) {
			// loop over the elements and locate theirs place in the grid

			// double dereferencing: 1)* to get item stored in the iterator 2)* to get Element stored within this pointer
			// it is derefereni
			if ( haveElement((**it), gr) ) {

				ackCounter++;

				Box<int, 3> ind = gr->locateBoxIndices( (*it)->getBoundingBox() );

				for ( int k = ind.minX(2); k <= ind.maxX(2); ++k ) {
					for ( int j = ind.minX(1); j <= ind.maxX(1); ++j ) {
						for ( int i = ind.minX(0); i <= ind.maxX(0); ++i ) {

							int                    inds[3] = {i,j,k};
							int                    linearI = gr->getLinearIndex(inds);
							std::vector<PetscReal> co      = gr->getLocalCoord(inds);
//							double                 p[3]    = {co[0], co[1], co[2]};
							Eigen::Vector3d        p;
							p << co[0], co[1], co[2];

							SignedDistance<double> v       = TriangleElement<double>::computeDistance((**it), p);
							double                 pd      = TriangleElement<double>::computePerpendicularDistance((**it), p);

							if ( v.dist > narrowBand ) {
								continue;
							}


							if (v.dist < eps) {
								data[linearI]                 = eps;
								perpendicualDistance[linearI] = eps;
								cellCounter++;

							}  else if ( fabs(v.dist - fabs(data[linearI])) < eps ) {
								// equal not funny
								if ( perpendicualDistance[linearI] <= pd ) {
									perpendicualDistance[linearI] = pd;
									data[linearI]                 = v.dist * v.sign;
									cellCounter++;
								}

							} else if ( v.dist < fabs(data[linearI]) ) {
								// less so ok than
								data[linearI]                 = v.dist * v.sign;
								perpendicualDistance[linearI] = pd;
								cellCounter++;
							}


						}
					}
				}

			} else {
				rejCounter++;
			}
		}



	} else {
//		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "RANK : %d, got nothing\n", rank);
//		PetscSynchronizedFlush(PETSC_COMM_WORLD);
	}



	// initialize boundaries, or not
	int selectedBoundaries = selectBoundaries(gr);

	int ghostLenX = info->gxm - info->xm;
	int ghostLenY = info->gym - info->ym;
	int ghostLenZ = info->gzm - info->zm;


	if ( selectedBoundaries & 1 ) {
		// this means lets initialize left side
		for (int i = 0; i < ghostLenX; ++i) {
			for (int j = 0; j < info->gym; ++j) {
				for (int k = 0; k < info->gzm; ++k) {

					PetscInt inds[dim]           = {i, j, k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

	if ( selectedBoundaries & 2 ) {
		// this means lets initialize right side
		for (int i = 0; i < ghostLenX; ++i) {
			for (int j = 0; j < info->gym; ++j) {
				for (int k = 0; k < info->gzm; ++k) {

					PetscInt inds[dim]           = {(info->gxm - 1) - i, j, k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

	if ( selectedBoundaries & 4 ) {
		// this means initialize the bottom side
		for (int j = 0; j < ghostLenY; ++j) {
			for (int i = 0; i < info->gxm; ++i) {
				for (int k = 0; k < info->gzm; ++k) {

					PetscInt inds[dim]           = {i, j, k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

	if ( selectedBoundaries & 8 ) {
		// this means initialize the top side
		for (int j = 0; j < ghostLenY; ++j) {
			for (int i = 0; i < info->gxm; ++i) {
				for (int k = 0; k < info->gzm; ++k) {

					PetscInt inds[dim]           = {i, (info->gym-1) - j, k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

	if ( selectedBoundaries & 16 ) {
		// this means initialize the front side
		for (int k = 0; k < ghostLenZ; ++k) {
			for (int i = 0; i < info->gxm; ++i) {
				for (int j = 0; j < info->gym; ++j) {

					PetscInt inds[dim]           = {i, j, k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

	if ( selectedBoundaries & 32 ) {
		// this means initialize the far side
		for (int k = 0; k < ghostLenZ; ++k) {
			for (int i = 0; i < info->gxm; ++i) {
				for (int j = 0; j < info->gym; ++j) {

					PetscInt inds[dim]           = {i, j, (info->gzm - 1) - k};
					long int linearInd           = gr->getLinearIndex(inds);
					std::vector<PetscReal> coord = gr->getLocalCoord(inds);
					Eigen::Vector3d p            = {coord[0], coord[1], coord[2]};

					if ( data[linearInd] < std::numeric_limits<double>::max() )
						continue;

					double dist                  = dataProvider->computeDistance(p);

					data[linearInd] = dist;
					cellCounter++;

				}
			}
		}
	}

// for debug purpose, set everything to 10, this is for post visualization
//	 for (long int i = 0; i < gr->getNumberOfLocalCells(); ++i)
//	 	if (data[i] == std::numeric_limits<double>::max())
//	 		data[i] = -10;
	
	delete perpendicualDistance;

}


template <typename DataProvider, int dim> int Initializer<DataProvider, dim>::selectBoundaries(Grid<double, dim>* gr) {
	// chooses which boundaries of the domain have to be initialized
	// ordering : [left, right, bottom, top, front, back]

	int selectedBoundaries_int = 0;

	for (int i = 0; i < dim; ++i) {
		if ( gr->getLocalGhostedMax(i) <= dataProvider->getMaxX(i) ) {
			selectedBoundaries_int |= ( 1 << (i*2 + 1) );
		}
		if ( gr->getLocalGhostedMin(i) >= dataProvider->getMinX(i) ) {
			selectedBoundaries_int |= ( 1 << (i*2) );
		}
	}

	return selectedBoundaries_int;

}


//template <typename DataProvider, int dim> std::vector<int> Initializer<DataProvider, dim>::selectBoundaries(Grid<double, dim>* gr) {
//	// chooses which boundaries of the domain have to be initialized
//	// ordering : [left, right, bottom, top, front, back]
//
//	std::vector<int> selectedBoundaries;
//	Box<double, 3> bb = dataProvider->getRegion();
//
//	if (dim >= 1) {
//		if ( bb.minX[0] <= gr->getLocalGhostedMin(0) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//
//		if (  bb.maxX[0] >= gr->getLocalGhostedMax(0) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//	}
//	if (dim >= 2) {
//		if ( bb.minX[1] <= gr->getLocalGhostedMin(1) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//
//		if ( bb.maxX[1] >= gr->getLocalGhostedMax(1) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//	}
//	if (dim == 3) {
//		if ( bb.minX[2] <= gr->getLocalGhostedMin(2) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//
//		if ( bb.maxX[2] >= gr->getLocalGhostedMax(2) ) {
//			selectedBoundaries.push_back(1);
//		} else {
//			selectedBoundaries.push_back(0);
//		}
//	}
//
//	return selectedBoundaries;
//
//}

template <typename DataProvider, int dim> bool Initializer<DataProvider, dim>::haveAllGeometry(Grid<double, dim>* gr) {
	// does the node ows the whole geometry
	// basically one needs to loop thourgh all the data and find bounding box and than check if
	// it fits inside the node.
	// for sphere : simple to do

	bool           owns = false;
	DMDALocalInfo* info = gr->getLocalInfo();

	if (dim == 1) {

	} else if (dim == 2) {

	}else if (dim == 3) {
		if ( !(gr->getLocalMin(0) <= dataProvider->getMinX(0) && dataProvider->getMaxX(0) <=  gr->getLocalMax(0)) ) {
			return false;
		}
		if ( !(gr->getLocalMin(1) <= dataProvider->getMinX(1) && dataProvider->getMaxX(1) <=  gr->getLocalMax(1)) ) {
			return false;
		}
		if ( !(gr->getLocalMin(2) <= dataProvider->getMinX(2) && dataProvider->getMaxX(2) <=  gr->getLocalMax(2)) ) {
			return false;
		}
	}

	return true;

}

template <typename DataProvider, int dim> bool Initializer<DataProvider, dim>::haveSomeGeometry(Grid<double, dim>* gr) {
	// does the node ows the whole geometry
	// basically one needs to loop thourgh all the data and find bounding box and than check if
	// it fits inside the node.
	// for sphere : simple to do

	if (dim == 1) {

	} else if (dim == 2) {

	} else if (dim == 3) {

		if ( !(gr->getLocalGhostedMin(0) <= dataProvider->getMaxX(0) && dataProvider->getMinX(0) <=  gr->getLocalGhostedMax(0)) ) {
			return false;
		}

		if ( !(gr->getLocalGhostedMin(1) <= dataProvider->getMaxX(1) && dataProvider->getMinX(1) <=  gr->getLocalGhostedMax(1)) ) {
			return false;
		}

		if ( !(gr->getLocalGhostedMin(2) <= dataProvider->getMaxX(2) && dataProvider->getMinX(2) <=  gr->getLocalGhostedMax(2)) ) {
			return false;
		}
	}

	return true;

}

template <typename DataProvider, int dim> inline bool
Initializer<DataProvider, dim>::haveElement(TriangleElement<double>& el, Grid<double, dim>* gr) const {

	bool           owns = false;

	if ( Box<double, 3>::intersects( el.getBoundingBox(), gr->getLocalGhostedRegion() ) )
		owns = true;
	else
		owns = false;

	return owns;

}

template <typename DataProvider, int dim> inline bool
 Initializer<DataProvider, dim>::haveCoordinate(const TriangleElement<double>& coord, Grid<double, dim>* gr) const {

	bool           owns = false;
	DMDALocalInfo* info = gr->getLocalInfo();

	if (dim == 1) {

	} else if (dim == 2) {

	}else if (dim == 3) {

		if ( coord.isPartlyInRegion( gr->getLocalGhostedRegion() ) )
			owns = true;
		else
			owns = false;

	}

	return owns;

}


#endif /* INITIALIZER_NEW_HPP_ */
