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


template <typename DataProvider, int dim>
class Initializer {

	// Actually all data should come from DataProvider class, even the distance function

	DataProvider* dataProvider;


public:
	Initializer();
	Initializer(DataProvider* provider): dataProvider(provider) {};
	~Initializer(){};

	// this method will do the magic and initialize the local part of the data accordingly to the
	// defined input stuff
	PetscReal* updateArray(PetscReal* data, Grid<double, dim>* gr);

};


///================================================
///               Implementation
///================================================


template <typename DataProvider, int dim> PetscReal* Initializer<DataProvider, dim>::updateArray(PetscReal* data, Grid<double, dim>* gr) {

	DMDALocalInfo* info = gr->getLocalInfo();

	int rank;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	long int sz      = gr->getNumberOfLocalCells();
	long int counter = 0;

	for (long int i = 0; i < sz; ++i)
		data[i] = std::numeric_limits<double>::max();

	// TODO: update grid so this is no longer needed
	for (int k = 0; k < info->gzm; ++k) {
		for (int j = 0; j < info->gym; ++j) {
			for (int i = 0; i < info->gxm; ++i) {



				int ii[3] = {i, j, k};
				int ind   = gr->getLinearIndex(ii);

				std::vector<PetscReal> co = gr->getLocalCoord(ii);

				Eigen::Vector3d p            = {co[0], co[1], co[2]};


//				double v = rank*10 + 1;

				double v = dataProvider->computeDistance(p, gr->getDx(0));

				std::cout << "Progress: " << double(counter) / double(sz) << std::endl;
				counter++;



				if (fabs(v) < fabs(data[ind]))
					data[ind]  = v;

			}
		}
	}

	return data;

}

#endif /* INITIALIZER_NEW_HPP_ */
