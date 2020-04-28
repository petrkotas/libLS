/*
 * Interface.h
 *
 *  Created on: Nov 7, 2012
 *      Author: petr
 */
#pragma once
#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "Grid.hpp"
#include <vector>
#include <string>
#include <petscdmadda.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <fstream>
#include <cstdio>

namespace lsm {
	#include "lsm_data_arrays.h"
}

template <typename type, int dim>
class Interface {
	// pointer to grid class that handles data distribution
	const Grid<type, dim>& gr;

	// local part of the data
	Vec localData;

	// global gathered data
	Vec globalData;

public:

	Interface(const Grid<type, dim>& _gr);
	virtual ~Interface();

	Vec getLocalData() const;
	void setLocalData(const std::vector<PetscReal>& data);

	PetscReal* getArray();
	void restoreArray(PetscReal* array);




	void writeData(char* file);

	void saveForMatlab(char* file);

	void saveForParaView(char* file);

	void dumpData();


	const Grid<type, dim>& getGrid() const;

};


//=========================================
//           Implementation
//=========================================





template <typename type, int dim>
Interface<type, dim>::Interface(const Grid<type, dim>& _gr) : gr(_gr){
	// create new data storage for interface evolution

	PetscErrorCode ierr;

	const DMDALocalInfo* info = gr.getLocalInfo();
	ierr                      = DMCreateLocalVector(info->da, &localData);


}

template <typename type, int dim>
Interface<type, dim>::~Interface() {
	// right now nothing
}


template <typename type, int dim>
Vec Interface<type, dim>::getLocalData() const {
	// return a copy of the internal data for whatever reason You need it

	return this->localData;

}

template <typename type, int dim>
void Interface<type, dim>::setLocalData(const std::vector<PetscReal>& data) {
	// set the local data from an arbitrary data field
	// probably better pass reference than copy of the data

	PetscErrorCode ierr;
	PetscReal*     arr;
	PetscInt       sz;

	ierr = VecGetSize(localData, &sz);

	if (data.size() != sz) {
		// probably a good time to throw exception but there is a problem,
		// this condition could hold on one node but will fail in other
		// question is the right way of handling this?
	}

	// otherwise everything went better than expected and let the copy begins

	ierr = VecGetArray(localData, &arr); // get array

	std::copy(data.begin(), data.end(), arr); // rewrite array

	ierr = VecRestoreArray(localData, &arr); // return array

}

template <typename type, int dim>
PetscReal* Interface<type, dim>::getArray() {
	// this method is used by solver only!!! it should not be used to access local data by any means
	// will have to look for better design

	PetscErrorCode ierr;
	PetscReal*     arr;

	ierr = VecGetArray(localData, &arr); // get array

	return arr;

}

template <typename type, int dim>
void Interface<type, dim>::restoreArray(PetscReal* array) {
	// this method shloud be used by solver only, it is higly volatile and could easily corrupt memory

	PetscErrorCode ierr;
	PetscInt       sz;

	ierr = VecRestoreArray(localData, &array);

}

//template <typename type, int dim> template <typename DataProvider>
//void Interface<type, dim>::initialize(Initializer<DataProvider, dim>* init) {
//	// loop thourg the data and initialize it
//
//	PetscReal* arr;
//	VecGetArray(localData, &arr);
//
//	int sz;
//	VecGetSize(localData, &sz);
//
//
//
//	init->updateArray(arr, gr); // call dataprovider and generate starting data
//
//	VecRestoreArray(localData, &arr);
//
////	VecView(localData, PETSC_VIEWER_STDOUT_SELF);
//
//}

template <typename type, int dim>
void Interface<type, dim>::dumpData() {
	const DMDALocalInfo* info = gr.getLocalInfo();	
	PetscErrorCode ierr;
	PetscInt       sz;
	PetscReal*     arr;

	ierr = VecGetSize(localData, &sz);
	ierr = VecGetArray(localData, &arr); // get array
	
	for (int i = 0; i < sz; ++sz) {
		std::cout << arr[i] << ", ";
	}
	std::cout << std::endl;
	
	ierr = VecRestoreArray(localData, &arr);


}

template <typename type, int dim>
void Interface<type, dim>::writeData(char* file) {

	int rank;
	const DMDALocalInfo* info = gr.getLocalInfo();


	long int length = gr.getM(0) * gr.getM(1) * gr.getM(2);


	Vec        glob, seq;
	VecScatter tolocalall, fromlocalall, sc;

	DMCreateGlobalVector(info->da, &glob);


	DMDAGlobalToNaturalAllCreate(info->da, &tolocalall);

	DMLocalToGlobalBegin(info->da, localData, INSERT_VALUES, glob);
	DMLocalToGlobalEnd(info->da, localData, INSERT_VALUES, glob);

	VecCreateSeq(PETSC_COMM_SELF, length, &seq);

	VecScatterBegin(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);


	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);



	if (rank == 0) {
		int d = dim;
		int dims[dim]  = {gr.getM(0), gr.getM(1), gr.getM(2)};
		double bl[dim] = {gr.getMin(0), gr.getMin(1), gr.getMin(2)};
		double tr[dim] = {gr.getMax(0), gr.getMax(1), gr.getMax(2)};
		double dx[dim] = {gr.getDx(0), gr.getDx(1), gr.getDx(2)};
				
		//lsm::Grid* toFile = lsm::createGridSetDxDyDz(dim, gr->getDx(0), gr->getDx(1), gr->getDx(2), bl, tr, lsm::VERY_HIGH);
		//lsm::writeGridToAsciiFile(toFile, "testingShit.ascii", NO_ZIP);

		FILE *fp;

		fp = fopen(file, "w");

		double * data;
		VecGetArray(seq, &data);

		fwrite(&d, sizeof(int), 1, fp);
		fwrite(dims, sizeof(int), 3, fp);
		fwrite(dx, sizeof(double), 3, fp);
		fwrite(bl, sizeof(double), 3, fp);
		fwrite(tr, sizeof(double), 3, fp);
		fwrite(data, sizeof(double), length, fp);

		fclose(fp);


		VecRestoreArray(seq, &data);

	}

}


template <typename type, int dim>
const Grid<type, dim>& Interface<type, dim>::getGrid() const {
	return gr;
}


template <typename type, int dim>
void Interface<type, dim>::saveForMatlab(char* file) {

	PetscErrorCode ierr;
	PetscViewer    view;
	DMDALocalInfo* info;

	info = gr->getLocalInfo();

	int size;

	ierr = DMCreateGlobalVector(info->da, &globalData);

	ierr = DMLocalToGlobalBegin(info->da, localData, INSERT_VALUES, globalData);

//	ierr = VecCreateSeq(PETSC_COMM_WORLD, dim, &size);
//	for (int i = 0; i < dim; ++i)
//		ierr = VecSetValue(size, i, gr->getM(i), INSERT_VALUES);

	ierr = DMLocalToGlobalEnd(info->da, localData, INSERT_VALUES, globalData);

	VecGetSize(globalData, &size);
	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "size: %d\n", size);


//	VecView(globalData, PETSC_VIEWER_STDOUT_WORLD);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &view);

	ierr = PetscObjectSetName((PetscObject)globalData, "phi");
	ierr = PetscObjectSetName((PetscObject)info->da, "da");
	ierr = DMView(info->da, view);
	ierr = VecView(globalData, view);

	ierr = PetscViewerDestroy(&view);

}

template <typename type, int dim>
void Interface<type, dim>::saveForParaView(char* file) {

	PetscErrorCode ierr;
	PetscViewer    view;
	DMDALocalInfo* info;

	info = gr->getLocalInfo();

	ierr = DMCreateGlobalVector(info->da, &globalData);

	ierr = DMLocalToGlobalBegin(info->da, localData, INSERT_VALUES, globalData);
	ierr = DMLocalToGlobalEnd(info->da, localData, INSERT_VALUES, globalData);

	ierr = PetscViewerCreate(PETSC_COMM_WORLD, &view);
	ierr = PetscViewerSetType(view, PETSCVIEWERASCII);
	ierr = PetscViewerFileSetName(view, file);
	ierr = PetscViewerSetFormat(view, PETSC_VIEWER_ASCII_VTK);
	ierr = PetscObjectSetName((PetscObject)globalData, "phi");
	ierr = DMView(info->da, view);
	ierr = VecView(globalData, view);
	ierr = PetscViewerDestroy(&view);

}
/*
template <typename type, int dim>
void Interface<type, dim>::saveToHDF5(char* file) {

	PetscErrorCode ierr;
	PetscViewer    view;
	DMDALocalInfo* info;

	info = gr->getLocalInfo();

	Vec size;

	ierr = DMCreateGlobalVector(info->da, &globalData);

	ierr = DMLocalToGlobalBegin(info->da, localData, INSERT_VALUES, globalData);
	ierr = DMLocalToGlobalEnd(info->da, localData, INSERT_VALUES, globalData);

	ierr = PetscViewerHDF5Open(PETSC_COMM_WORLD, file, FILE_MODE_WRITE, &view);
//
	ierr = PetscObjectSetName((PetscObject)globalData, "phi");
//	ierr = PetscObjectSetName((PetscObject)info->da, "da");
//	ierr = DMView(info->da, view);
	ierr = VecView(globalData, view);

	ierr = PetscViewerDestroy(&view);

}
*/

#endif /* INTERFACE_H_ */
