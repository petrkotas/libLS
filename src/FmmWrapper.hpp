/*
 * FmmWrapper.hpp
 *
 *  Created on: Nov 13, 2012
 *      Author: petr
 */
#pragma once
#ifndef FMMWRAPPER_HPP_
#define FMMWRAPPER_HPP_

#include <limits>
#include "LSMLIB_config.h"
#include "lsm_fast_marching_method.h"
#include "Interface.hpp"
#include "Grid.hpp"

class FmmWrapper {
public:

	template <int dim>
	static bool solveEikonalEquation(const Grid<double, dim>& gr, Interface<double, dim>* phi, double maxShift);

};


///===================================
///          Implementation
///===================================

template <int dim> bool FmmWrapper::solveEikonalEquation(const Grid<double, dim>& gr, Interface<double, dim>* phi, double maxShift) {
	// wrappes solving of local part of eikonal equation

	LSMLIB_REAL* arr                      = (LSMLIB_REAL*) phi->getArray();
	LSMLIB_REAL* mask                     = 0;
	int          spatial_derivative_order = 1;
	LSMLIB_REAL* speed                    = new LSMLIB_REAL[gr.getNumberOfLocalCells()];
//	LSMLIB_REAL* dist                     = new LSMLIB_REAL[grid->getNumberOfLocalCells()];
	int          grid_dims[dim];
	double       dx[dim];
	double       mindx                    = std::numeric_limits<double>::max();
	double       shift                    = maxShift;

	// switch the dimension so liblsm computes it corectly
//	dx[0] = grid->getDx(1);
//	dx[1] = grid->getDx(0);
//	dx[2] = grid->getDx(2);
//
//	grid_dims[0] = grid->getLocalM(1);
//	grid_dims[1] = grid->getLocalM(0);
//	grid_dims[2] = grid->getLocalM(2);

	for (int i = 0; i < dim; ++i) {
		dx[i] = gr.getDx(i);
		grid_dims[i] = gr.getLocalM(i);
		if (mindx > dx[i])
			mindx = dx[i];

//		PetscPrintf(PETSC_COMM_WORLD, "dx[%d] = %f\n", i, dx[i]);
	}

//	PetscPrintf(PETSC_COMM_WORLD, "mindx = %f\n", mindx);

	for (long int i = 0; i < gr.getNumberOfLocalCells(); ++i) {
		if (arr[i] >= 1E+15) {
			arr[i] = -1;
		} else if (arr[i] <= 0)
			arr[i] = -arr[i];
		else
			arr[i] += shift;
	}

	for (int k = 0; k < grid_dims[2]; ++k) {
		for (int j = 0; j < grid_dims[1]; ++j) {
			for (int i = 0; i < grid_dims[0]; ++i) {
				speed[i + j*grid_dims[0] + k*grid_dims[0]*grid_dims[1]] = 1;
			}
		}
	}

	if (dim == 2) {
		//2d case
		int ierr = solveEikonalEquation2d(arr, speed, mask, spatial_derivative_order, grid_dims, dx);

	} else if(dim == 3){
		// 3d case
		int ierr = solveEikonalEquation3d(arr, speed, mask, spatial_derivative_order, grid_dims, dx);
		if (ierr != 0)
			std::cout << "error" << std::endl;

	}

	for (long int i = 0; i < gr.getNumberOfLocalCells(); ++i) {
		if (arr[i] < shift)
			arr[i] = -arr[i];
		else
			arr[i] -= shift;
	}

	phi->restoreArray((PetscReal*) arr);

	delete speed;
//	delete dist;

	return true;

}



#endif /* FMMWRAPPER_HPP_ */
