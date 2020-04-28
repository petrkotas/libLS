/*
 * LSMWritter.hpp
 *
 *  Created on: Dec 12, 2013
 *      Author: petr
 */
#pragma once
#ifndef LSMWRITTER_HPP_
#define LSMWRITTER_HPP_

#include <petscdmadda.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <string>

#include "Grid.hpp"
#include "Interface.hpp"

namespace lsm {
	#include "lsm_data_arrays.h"
}

class LSMWritter {
public:
	template <typename type, int dim>
	static void writeData(const Interface<type, dim>& data, char* file) {

		int rank;
		const Grid<type, dim>& gr = data.getGrid();
		const DMDALocalInfo* info = gr.getLocalInfo();


		long int length = gr.getM(0) * gr.getM(1) * gr.getM(2);


		Vec        glob, seq;
		VecScatter tolocalall, fromlocalall, sc;

		DMCreateGlobalVector(info->da, &glob);


		DMDAGlobalToNaturalAllCreate(info->da, &tolocalall);

		DMLocalToGlobalBegin(info->da, data.getLocalData(), INSERT_VALUES, glob);
		DMLocalToGlobalEnd(info->da, data.getLocalData(), INSERT_VALUES, glob);

		VecCreateSeq(PETSC_COMM_SELF, length, &seq);

		VecScatterBegin(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(tolocalall, glob, seq, INSERT_VALUES, SCATTER_FORWARD);


		MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

		std::string fname_g(file);
		std::string fname_d(file);

		fname_g += ".gr";
		fname_d += ".data";

		if (rank == 0) {
			int d = dim;
			int dims[dim]  = {gr.getM(0)  , gr.getM(1)  , gr.getM(2)};
			double bl[dim] = {gr.getMin(0), gr.getMin(1), gr.getMin(2)};
			double tr[dim] = {gr.getMax(0), gr.getMax(1), gr.getMax(2)};
			double dx[dim] = {gr.getDx(0) , gr.getDx(1) , gr.getDx(2)};

			std::cout << dx[0] << ", " << dx[1] << ", " << dx[2] << std::endl;

			bl[0] += 2*dx[0]; bl[1] += 2*dx[1]; bl[2] += 2*dx[2];
			tr[0] -= 2*dx[0]; tr[1] -= 2*dx[1]; tr[2] -= 2*dx[2];
			dims[0] -= 4; dims[1] -= 4; dims[2] -= 4;

			double * data;
			VecGetArray(seq, &data);

			lsm::Grid* toFile = lsm::createGridSetDx(dim, dx[0], bl, tr, lsm::LOW);
			lsm::writeGridToAsciiFile(toFile, &fname_g[0], NO_ZIP);
			lsm::writeDataArray(data, toFile, &fname_d[0], NO_ZIP);

			VecRestoreArray(seq, &data);

		}

	}

};

#endif /* LSMWRITTER_HPP_ */
