/*
 * LSMWritter.hpp
 *
 *  Created on: Dec 12, 2013
 *      Author: petr
 */
#pragma once
#ifndef BINWRITTER_HPP_
#define BINWRITTER_HPP_

#include <petscdmadda.h>
#include <petscdmda.h>
#include <petscsys.h>
#include <petscviewer.h>
#include <string>

#include "Grid.hpp"
#include "Interface.hpp"



class BINWritter {
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

		std::string fname_d(file);

		fname_d += ".data";

		if (rank == 0) {
			int d = dim;
			int dims[dim]  = {gr.getM(0)  , gr.getM(1)  , gr.getM(2)};
			double bl[dim] = {gr.getMin(0), gr.getMin(1), gr.getMin(2)};
			double tr[dim] = {gr.getMax(0), gr.getMax(1), gr.getMax(2)};
			double dx[dim] = {gr.getDx(0) , gr.getDx(1) , gr.getDx(2)};

			FILE *fp;

			fp = fopen(&fname_d[0], "w");

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

};

#endif /* BINWRITTER_HPP_ */
