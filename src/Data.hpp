/*
 * Interface.h
 *
 *  Created on: Nov 7, 2012
 *      Author: petr
 */

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "Grid.hpp"
//#include "Initializer_bruteforce.hpp"
#include "Initializer_new.hpp"
//#include "Initializer.hpp" // initializer for the sphere
#include <vector>
#include <string>
#include <petscdmadda.h>
#include <petscdmda.h>
#include <petscdm.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <fstream>

template <typename type, int dim>
class Interface {

	// pointer to grid class that handles data distribution
	Grid<type, dim> *gr;

	// local part of the data
	Vec localData;

	// global gathered data
	Vec globalData;

public:

	Interface(Grid<type, dim> *gr);
	virtual ~Interface();

	std::vector<PetscReal> getLocalData();
	void setLocalData(const std::vector<PetscReal>& data);

	PetscReal* getArray();
	void restoreArray(PetscReal* array);

	template <typename DataProvider>
	void initialize(Initializer<DataProvider, dim>* init);

	void writeData(char* file);

	void saveForMatlab(char* file);

	void saveForParaView(char* file);

	void saveToHDF5(char* file);

	Grid<type, dim>* getGrid();

};


//=========================================
//           Implementation
//=========================================

template <typename type, int dim>
Interface<type, dim>::Interface(Grid<type, dim> *gr) {
	// create new data storage for interface evolution

	PetscErrorCode ierr;

	DMDALocalInfo* info = gr->getLocalInfo();
	this->gr            = gr;
	ierr                = DMCreateLocalVector(info->da, &localData);

}

template <typename type, int dim>
Interface<type, dim>::~Interface() {
	// right now nothing
}


template <typename type, int dim>
std::vector<PetscReal> Interface<type, dim>::getLocalData() {
	// return a copy of the internal data for whatever reason You need it

	PetscErrorCode         ierr;
	PetscScalar*           arr;
	PetscInt               sz;

	ierr = VecGetSize(localData, &sz);
	ierr = VecGetArray(localData, &arr);

	std::vector<PetscReal> retVec(arr, arr+sz);

//	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "sz = %d\n", sz);
//	PetscSynchronizedPrintf(PETSC_COMM_WORLD, "retVec.size() = %d\n", retVec.size());
//	//	for (int i = 0; i < retVec.size(); i++) {
////		PetscSynchronizedPrintf(PETSC_COMM_WORLD, "sz = %d\n", sz);
////	}
//	PetscSynchronizedFlush(PETSC_COMM_WORLD);

	ierr = VecRestoreArray(localData, &arr);

	return retVec;

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

template <typename type, int dim> template <typename DataProvider>
void Interface<type, dim>::initialize(Initializer<DataProvider, dim>* init) {
	// loop thourg the data and initialize it

	PetscReal* arr;
	VecGetArray(localData, &arr);

	int sz;
	VecGetSize(localData, &sz);

	init->updateArray(arr, gr); // call dataprovider and generate starting data

	VecRestoreArray(localData, &arr);

//	VecView(localData, PETSC_VIEWER_STDOUT_SELF);

}


template <typename type, int dim>
void Interface<type, dim>::writeData(char* file) {

	int rank;
	DMDALocalInfo* info = gr->getLocalInfo();

	Vec        glob, seq;
    VecScatter tolocalall, fromlocalall, sc;

    DMCreateGlobalVector(info->da, &glob);
    VecCreateSeq(PETSC_COMM_SELF, gr->getM(0)*gr->getM(1)*gr->getM(2), &seq);

    DMDAGlobalToNaturalAllCreate(info->da, &tolocalall);

    DMLocalToGlobalBegin(info->da, localData, INSERT_VALUES, glob);
    DMLocalToGlobalEnd(info->da, localData, INSERT_VALUES, glob);

//    VecScatterCreateToZero(glob, &sc, &seq);

//    VecView(glob, PETSC_VIEWER_STDOUT_WORLD);

	VecScatterBegin(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);


	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

	// ASSUMES: for all dimensions, ghostbox margin on two sides are equal and >= 1 !!!

	// Check that dx is equal in all dimensions
	if ((gr->getDx(0) != gr->getDx(1)) || (gr->getDx(1) != gr->getDx(2)) || (gr->getDx(0) != gr->getDx(2))) {
		// lol wut, UG needs the equidistant net for all dimensions
	}


	if (rank == 0) {
		// write the data to the file

		std::ofstream oFile;
		oFile.open(file);

		if (!oFile.is_open()) {
			std::cerr << "Problem opening file." << std::endl;
			// suxx, shit happend
		}

		long int numberOfCells = gr->getNumberOfCells();


		oFile << "{" << std::endl;
		oFile << "Gitterweite" << std::endl;
		oFile << gr->getDx(0) << std::endl;
		oFile << "BoundingBox" << std::endl;
		oFile << "nmaxx " << gr->getMax(0) << std::endl;
		oFile << "nminx " << gr->getMin(0) << std::endl;
		oFile << "nmaxy " << gr->getMax(1) << std::endl;
		oFile << "nminy " << gr->getMin(1) << std::endl;
		oFile << "nmaxz " << gr->getMax(2) << std::endl;
		oFile << "nminz " << gr->getMin(2) << std::endl;
		oFile << "Abstand " << numberOfCells << std::endl;

		PetscReal* arr;
		VecGetArray(seq, &arr);

		for (int z = 0; z < gr->getM(2); ++z) {
			for (int x = 0; x < gr->getM(0); ++x) {
				for (int y = 0; y < gr->getM(1); ++y) {
					long int idx = x + y*gr->getM(0) + z*gr->getM(0)*gr->getM(1);
					oFile << arr[idx] << std::endl;
				}
			}
		}

		VecRestoreArray(seq, &arr);

		oFile << "}";
		oFile.close();

	}

}


template <typename type, int dim>
Grid<type, dim>* Interface<type, dim>::getGrid() {
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


#endif /* INTERFACE_H_ */
