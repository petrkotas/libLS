#pragma once
#ifndef _INITIALIZER_HPP_
#define _INITIALIZER_HPP_

#include <list>
#include <set>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include <petsctime.h>

#include "Grid.hpp"
#include "Interface.hpp"
#include "BoundingBox.hpp"
#include "TriangleElement.hpp"
#include "Geometry.hpp"

#include "tictoc.hpp"
#include "utility.h"



/// compares two signed distances
/// the signed distance and also prepedicular distance is compared
void sdfMin(cDistance* in, cDistance* inout, int len, MPI_Datatype* dptr) {
	int i = 0;
	for(i = 0; i < len; ++i) {
		double ax  = fabs(in[i].signedDistance), ay = fabs(inout[i].signedDistance);
		double eps = my_eps;
		if ( fabs(ax - ay) < eps ) {
			if (in[i].pDistance > inout[i].pDistance) {
				inout[i] = in[i];
			}
		} else if (ax > ay) {
			inout[i] = in[i];
		}
	}
}


class Initializer {
public:
	Initializer() {};

	// This is where the magic happens and the initializator puts the data in
	template <typename type, int dim>
	void operator() (Geometry& geom, Interface<type, dim>& interface, int groupID, int nGroups);


private:
	// some states here and utility functions needed to do the jolb

	MPI_Datatype mpi_cDist_type;
	MPI_Op mpi_sdfmin;

	void createTypes();
	template <typename type, int dim>
	void initAll(const Grid<type, dim>& gr, Geometry& geom, double*** data_ptr);
	template <typename type, int dim>
	void initBoundaries_nlb(const Grid<type, dim>& gr, Geometry& geom, double*** data_ptr, int boundaries);
	template <typename type, int dim>
	void initBoundaries_lb(const Grid<type, dim>& gr, Geometry& geom, double*** data_ptr, int boundaries, int groupID, int nGroups);

	int selectBoundaries(const Box<double, 3>& geometryAABB, const Box<double, 3>& procAABB);

	template <typename type, int dim>
	void putDataToBoundaries(const Grid<type, dim>& gr, double*** data_ptr, cDistance* distanceData, int cellNum, int boundaries);

	template <typename type, int dim>
	void putLocalDataInside(const Grid<type, dim>& gr, double*** data_ptr, const std::vector<int>& localTriangles, const Geometry& geom, double*** perpendicualDistance);

};

///
/// Main routine
///


void Initializer::createTypes() {
	// TODO this should be rewritten as create type placed in cDistance and create minFcn
	// define according data type for cDistance in MPI
	const int    nitems = 2;
	int          blocklengths[2] = {1, 1};
	MPI_Datatype types[2] = {MPI_DOUBLE, MPI_DOUBLE};
	MPI_Aint     offsets[2];

	offsets[0] = offsetof(cDistance, signedDistance);
	offsets[1] = offsetof(cDistance, pDistance);

	MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_cDist_type);
	MPI_Type_commit(&mpi_cDist_type);
	// end type definition

	// commit custom mpi reduce operation
	MPI_Op_create((MPI_User_function *) sdfMin, 1, &mpi_sdfmin);
	// end custom reduce operation commit
}

template <typename type, int dim>
void Initializer::initAll(const Grid<type, dim>& gr, Geometry& geom, double*** data_ptr) {

	DM cda;
	Vec gc;
	DMDACoor3d ***coors;

	DMGetCoordinatesLocal(gr.getDA() ,&gc);
	DMGetCoordinateDM(gr.getDA(), &cda);
	DMDAVecGetArray(cda, gc, &coors);


	int       x, y, z, m, n, p;
	int       cellCounter = 0, gh = 1;

	DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);
	double data_ref = 0.0, data_acc = 0.0;
	for (int i = x; i < x+m; ++i) {
		for (int j = y; j < y+n; ++j) {
			for (int k = z; k < z+p; ++k) {
				data_acc =
						geom.computeDistance(
							Eigen::Vector3d(coors[k][j][i].x,
											coors[k][j][i].y,
											coors[k][j][i].z),
							true
						);

				data_ref = geom.computeDistance(
							Eigen::Vector3d(coors[k][j][i].x,
											coors[k][j][i].y,
											coors[k][j][i].z),
							false
						);

			}
		}
	}

}

template <typename type, int dim>
void Initializer::initBoundaries_nlb(const Grid<type, dim>& gr, Geometry& geom, double*** data_ptr, int boundaries) {

	DM cda;
	Vec gc;
	DMDACoor3d ***coors;

	DMGetCoordinatesLocal(gr.getDA() ,&gc);
	DMGetCoordinateDM(gr.getDA(), &cda);
	DMDAVecGetArray(cda, gc, &coors);


	int       x, y, z, m, n, p;
	int       cellCounter = 0, gh = 1;

	DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);

	int tmpSum = 0;
	if ( boundaries & 1 ) {
		tmpSum+= gh*n*p;
		// this means lets initialize left side
		for (int i = x; i < x+gh; ++i) {
			for (int j = y; j < y+n; ++j) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 2 ) {
		tmpSum += gh*n*p;
		// this means lets initialize right side
		for (int i = x+m-1; i > x+m-1-gh; --i) {
			for (int j = y; j < y+n; ++j) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 4 ) {
		tmpSum += gh*m*p;
		// this means initialize the bottom side
		for (int j = y; j < y+gh; ++j) {
			for (int i = x; i < x+m; ++i) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 8 ) {
		tmpSum += gh*m*p;
		// this means initialize the top side
		for (int j = y+n-1; j > y+n-1-gh; --j) {
			for (int i = x; i < x+m; ++i) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 16 ) {
		tmpSum += gh*m*n;
		// this means initialize the front side
		for (int k = z; k < z+gh; ++k) {
			for (int i = x; i < x+m; ++i) {
				for (int j = y; j < y+n; ++j) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 32 ) {
		tmpSum += gh*m*n;
		// this means initialize the far side
		for (int k = z+p-1; k > z+p-1-gh; --k) {
			for (int i = x; i < x+m; ++i) {
				for (int j = y; j < y+n; ++j) {
					data_ptr[k][j][i] =
							geom.computeDistance(
								Eigen::Vector3d(coors[k][j][i].x,
												coors[k][j][i].y,
												coors[k][j][i].z),
								true
							);
					cellCounter++;
				}
			}
		}
	}

}

/* 
template <typename type, int dim>
void Initializer::initBoundaries_lb(const Grid<type, dim>& gr,
		Geometry& geom,
		double*** data_ptr,
		int boundaries,
		int groupID, int nGroups) {

	tictoc tick;
	int myWorldRank, nProcsInWorld;
	int myBegin, myEnd;   // local indices start, local indices stop
	int myBoundaries;
	MPI_Comm_rank(MPI_COMM_WORLD, &myWorldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcsInWorld);

	tick.tic();
	// prepare memory to store the list of all boundary cells
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > boundaryCellList;
	boundaryCellList.reserve( nProcsInWorld * 4 * gr.getM(0)*gr.getM(1) );

	// this is where the results go (final results)
	std::vector<cDistance> vDistance;
	// create the fucking list
	for (int p = 0; p < nProcsInWorld; ++p) {
		// loop through the processes and get the proc bounding box
		Box<double, 3> procBB    = gr.getNodeSpan(p);
		int boundariesToInit     = selectBoundaries(geom.aabb, procBB);

		std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > l =
				gr.getNodeBoundaryCells(p, boundariesToInit);

		// ok I need to remember what part of the boundary is mine, otherwise I am screwed
		if (myWorldRank == p) {
			myBegin      = boundaryCellList.size();
			myEnd        = myBegin + l.size();
			myBoundaries = boundariesToInit;
		}

		boundaryCellList.insert( boundaryCellList.end(), l.begin(), l.end() );
	}

	tick.toc(myWorldRank, MAKECELLLIST, "Create boundary list");

	// number of cells in the boundary list for the WHOLE FUCKING BOUNDARY
	int nBCells = boundaryCellList.size();
	std::cout << "ncells: " << nBCells << std::endl;
	tick.tic();
	// loadbalance the problem.
	// first create groups for the problem
	MPI_Comm groupCommunicator = MPI_COMM_WORLD;
	int nProcsInGroup;
	int myGroupRank;
	if ( nGroups > 1 ) {
		// but if the number of groups is larger than one we need to create the right communicator.
		// luckily the group assigment is already sorted in setup phase.
		MPI_Comm_split(MPI_COMM_WORLD, groupID, myWorldRank, &groupCommunicator);
	}
	MPI_Comm_rank(groupCommunicator, &myGroupRank);
	MPI_Comm_size(groupCommunicator, &nProcsInGroup);
	// now each group will split its workload of boundary cells and start doing stuff

	int numberOfCellsPerTask     = geom.getNumberOfCellsPerDomainLocal();
	assert( numberOfCellsPerTask >= 0 );

	int max_numberOfCellsPerTask = -std::numeric_limits<int>::max();
	size_t totalCellsToFillTheGaps  = 0;
	int exchangedNumberOfCellsPerTask[nProcsInGroup];

	MPI_Allgather(&numberOfCellsPerTask, 1, MPI_INT, exchangedNumberOfCellsPerTask, nProcsInGroup, MPI_INT, groupCommunicator);

	max_numberOfCellsPerTask = *std::max_element(exchangedNumberOfCellsPerTask, exchangedNumberOfCellsPerTask + nProcsInGroup);
	assert( max_numberOfCellsPerTask >= 0 ); // maximum should not be negative

	double diffs[nProcsInGroup];// diffs is the number of cells going to be added as balance for
	for (int p = 0 ; p < nProcsInGroup; ++p) {
		diffs[p] = max_numberOfCellsPerTask - exchangedNumberOfCellsPerTask[p];
		assert( diffs[p] >= 0 );
		totalCellsToFillTheGaps += size_t(max_numberOfCellsPerTask - exchangedNumberOfCellsPerTask[p]);
		 if (myWorldRank == 0)
			 std::cout << "DIFFS : " << (size_t)diffs[p] << "; MAX : " << max_numberOfCellsPerTask << "; EXCHANGED: " << exchangedNumberOfCellsPerTask[p] << std::endl;
	}
	if (myWorldRank == 0)
		 std::cout << "SUM: " << totalCellsToFillTheGaps << std::endl;
	size_t tmpSum = 0; double pSum = 0.0;
	
	if ( totalCellsToFillTheGaps > 0 ) {
		for (int p = 0 ; p < nProcsInGroup; ++p) {
			pSum += diffs[p]/double(totalCellsToFillTheGaps);
			tmpSum += size_t(diffs[p]);

			 if (myWorldRank == 0)
				 std::cout << "PROC: " << p << " DIFF: " << double(diffs[p])  << " / " << double(totalCellsToFillTheGaps) << " - " << diffs[p]/double(totalCellsToFillTheGaps) << " - " << pSum
				 << " - " << tmpSum << std::endl;
			
			diffs[p] = double(diffs[p]) / double(totalCellsToFillTheGaps);
		}
		 if (myWorldRank == 0)
			 std::cout << "TMPSUM: " << tmpSum << std::endl;
		assert( tmpSum == totalCellsToFillTheGaps ); // the percentage is not allowed to be greater then 100%
	}

	int theRest = nBCells - totalCellsToFillTheGaps;
	if ( theRest < 0 ) {
		// to make sure that assigment of negative shit is not happening
		theRest                 = 0;
		totalCellsToFillTheGaps = nBCells;
	}
	
	int start = 0, finish = 0;
	int send_offsets[nProcsInGroup];
	int send_lens[nProcsInGroup];
	int debugSum = 0;


	for (int p = 0; p < nProcsInGroup; ++p) {
		start  = finish;

		finish = (p == nProcsInGroup-1) ? nBCells : start + int( diffs[p]*double(totalCellsToFillTheGaps) ) + int( double(theRest)/double(nProcsInGroup) );
		 
		send_offsets[p] = start;
		send_lens[p]    = finish - start;
		debugSum += send_lens[p];
	}
	tick.toc(myWorldRank, SETUPGROUPS, "Setup groups");
	assert( debugSum <= nBCells );
	assert( send_lens[nProcsInGroup-1] + send_offsets[nProcsInGroup-1] == nBCells );	


	tick.tic();
	// std::cout << "send_offsets = " << send_offsets[myGroupRank] << ", send_lens = " << send_lens[myGroupRank] << ", nBCells = " << nBCells << std::endl;
	assert( send_lens[myGroupRank] >= 0 );
	assert( send_offsets[myGroupRank] >= 0 );
	assert( (send_offsets[myGroupRank] + send_lens[myGroupRank]) <= nBCells );
//	std::vector<cDistance> localDistances; localDistances.reserve( send_lens[myGroupRank] );
	cDistance localDistances[ send_lens[myGroupRank] ];
	// compute the distance within the group


	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >::const_iterator begin_range = boundaryCellList.begin() + send_offsets[myGroupRank];
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> >::const_iterator end_range = boundaryCellList.begin() + send_offsets[myGroupRank] + send_lens[myGroupRank];
	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > subList(begin_range, end_range);
	std::cout << "allocated for locality check" << std::endl;



	geom.computeDistance(subList, localDistances);
	tick.toc(myWorldRank, COMPUTEDIST, "Compute bval dist");

	boundaryCellList.clear();
	

	// put data together. That means call allgatherv and synchronize the shit.

	tick.tic();
	assert( nBCells > 0 );
	vDistance.resize( nBCells );
	MPI_Allgatherv(&localDistances[0], send_lens[myGroupRank], mpi_cDist_type,  &vDistance[0], send_lens, send_offsets, mpi_cDist_type, groupCommunicator);

	tick.toc(myWorldRank, ALLGATHER, "Allgatherv vector");
	// now the data is placed together in one nice vector
	std::cout << "dumping vDistance: " << std::endl;
	for (auto d: vDistance) {
		std::cout << d.signedDistance << ", ";
	} 
	std::cout << std::endl;

	// synchronize the groups by reduce operation. This is only needed when you have more then one group

	tick.tic();
	if ( nGroups > 1 ) {
//		mpi::all_reduce(world, &vDistance[0], nBCells, &vDistance[0], sdfMinimum<cDistance>());
		MPI_Allreduce(&vDistance[0], &vDistance[0], nBCells, mpi_cDist_type, mpi_sdfmin, MPI_COMM_WORLD);
	}
	tick.toc(myWorldRank, REDUCE, "Reduce between groups");

	tick.tic();
	// step THREE put into the grid
	assert( myBegin >= 0 );
	assert( (myEnd-myBegin) < nBCells );
	putDataToBoundaries(gr, data_ptr, &vDistance[myBegin], (myEnd-myBegin), myBoundaries);

	vDistance.clear();
	tick.toc(myWorldRank, PUTTOBOUNDARIES, "Place to boundaries");

}
*/

template <typename type, int dim>
void Initializer::operator ()(Geometry& geom, Interface<type, dim>& interface, int groupID, int nGroups) {


	std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > boundaryCells;    // xyz coordinated to compute the distance to

	int                        myBoundaries;     // local boundaries to init

	MPI_Group worldGroup;
	MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
	// perpendicular distance computed towards the triangle, it is used when we are no longer sure, which normal shall we use

	double*** perpendicualDistance;


	// initialize the data

	Vec       localData = interface.getLocalData();
	double*** data_ptr;
	int       x, y, z, m, n, p;
	int       vecSize;

	VecGetLocalSize(localData, &vecSize);
	assert( vecSize > 0 );

	int myWorldRank;
	int nProcsInWorld;

	MPI_Comm_rank(MPI_COMM_WORLD, &myWorldRank);
	MPI_Comm_size(MPI_COMM_WORLD, &nProcsInWorld);


	const Grid<type, dim>& gr = interface.getGrid();

	DMDAVecGetArray(gr.getDA(), localData, &data_ptr);
	DMDAGetArray(gr.getDA(), PETSC_TRUE, &perpendicualDistance);
	DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);

	assert( vecSize == m*n*p );

	for (int k = z, pk = 0; k < z+p; ++k, ++pk) {
		for (int j = y, pj = 0; j < y+n; ++j, ++pj) {
			for (int i = x, pi = 0; i < x+m; ++i, ++pi) {
				data_ptr[k][j][i]             =  std::numeric_limits<double>::max();
				perpendicualDistance[k][j][i] = -std::numeric_limits<double>::max();
			}
		}
	}

	// step FIVE compute triangle distance inside the narrowband

	putLocalDataInside(gr, data_ptr, geom.localElements, geom, perpendicualDistance);

	DMDARestoreArray(gr.getDA(), PETSC_TRUE, &perpendicualDistance);

	// =====================================================================================================

	if ( nProcsInWorld == 1 ) {
//		DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);
//		for (int k = z; k < z+p; ++k) {
//			for (int j = y; j < y+n; ++j) {
//		 		for (int i = x; i < x+m; ++i) {
//		 			if ( data_ptr[k][j][i] >= 10E+5 )
//		 				data_ptr[k][j][i] = -10;
//				}
//			}
//		}
		DMDAVecRestoreArray(gr.getDA(), localData, &data_ptr);
		return;
	}

	int boundariesToInit = selectBoundaries(geom.aabb, gr.getNodeSpan(myWorldRank));
	initBoundaries_nlb(gr, geom, data_ptr, boundariesToInit);
	// initAll(gr, geom, data_ptr);

//=======================================================================================
//	  DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);
//	  for (int k = z; k < z+p; ++k) {
//	  	for (int j = y; j < y+n; ++j) {
//	  		for (int i = x; i < x+m; ++i) {
//	  			if ( data_ptr[k][j][i] >= 10E+5 ) {
//	  				data_ptr[k][j][i] = -10;
//	  			}
//	  		}
//	  	}
//	  }
//=======================================================================================

	DMDAVecRestoreArray(gr.getDA(), localData, &data_ptr);

}


int Initializer::selectBoundaries(const Box<double, 3>& geometryAABB, const Box<double, 3>& procAABB) {
	// chooses which boundaries of the domain have to be initialized
	// ordering : [left, right, bottom, top, front, back]

	int selectedBoundaries_int = 0;

	for (int i = 0; i < 3; ++i) {
		//gr->getLocalGhostedMax(i) <= dataProvider->getMaxX(i)
		if ( procAABB.maxX(i) <= geometryAABB.maxX(i) ) {
			selectedBoundaries_int |= ( 1 << (i*2 + 1) );
		}
		// gr->getLocalGhostedMin(i) >= dataProvider->getMinX(i)
		if ( procAABB.minX(i) >= geometryAABB.minX(i) ) {
			selectedBoundaries_int |= ( 1 << (i*2) );
		}
	}

	return selectedBoundaries_int;

}


template <typename type, int dim>
void Initializer::putDataToBoundaries(const Grid<type, dim>& gr, double*** data_ptr, cDistance* distanceData, int cellNum, int boundaries) {

	int       x, y, z, m, n, p;
	int       cellCounter = 0, gh = 1;

	DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);
//	DMDAGetCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);
	int tmpSum = 0;
	if ( boundaries & 1 ) {
		tmpSum+= gh*n*p;
		// this means lets initialize left side
		for (int i = x; i < x+gh; ++i) {
			for (int j = y; j < y+n; ++j) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i+1] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 2 ) {
		tmpSum += gh*n*p;
		// this means lets initialize right side
		for (int i = x+m-1; i > x+m-1-gh; --i) {
			for (int j = y; j < y+n; ++j) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j][i-1] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 4 ) {
		tmpSum += gh*m*p;
		// this means initialize the bottom side
		for (int j = y; j < y+gh; ++j) {
			for (int i = x; i < x+m; ++i) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j+1][i] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 8 ) {
		tmpSum += gh*m*p;
		// this means initialize the top side
		for (int j = y+n-1; j > y+n-1-gh; --j) {
			for (int i = x; i < x+m; ++i) {
				for (int k = z; k < z+p; ++k) {
					data_ptr[k][j-1][i] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 16 ) {
		tmpSum += gh*m*n;
		// this means initialize the front side
		for (int k = z; k < z+gh; ++k) {
			for (int i = x; i < x+m; ++i) {
				for (int j = y; j < y+n; ++j) {
					data_ptr[k+1][j][i] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}
	if ( boundaries & 32 ) {
		tmpSum += gh*m*n;
		// this means initialize the far side
		for (int k = z+p-1; k > z+p-1-gh; --k) {
			for (int i = x; i < x+m; ++i) {
				for (int j = y; j < y+n; ++j) {
					data_ptr[k-1][j][i] = distanceData[cellCounter].signedDistance;
					cellCounter++;
				}
			}
		}
	}

	assert( cellNum == tmpSum );

}

template <typename type, int dim>
void Initializer::putLocalDataInside(const Grid<type, dim>& gr, double*** data_ptr, const std::vector<int>& localTriangles, const Geometry& geom, double*** perpendicualDistance) {
	int       x, y, z, m, n, p;
	double    narrowBand = gr.getDx(0) * 3;
	double    minDX      = Eigen::Array3d(gr.getDx(0), gr.getDx(1), gr.getDx(2)).minCoeff();
	double    eps        = my_eps;

	DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);

//	 std::cout << "number of triangles = " << localTriangles.size() << std::endl;
//	 for (auto i: localTriangles) {
//		 std::cout << i << ", ";
//	 }
//	 std::cout << std::endl;

	if ( localTriangles.size() == 0 )
		return;


	for (auto index: localTriangles) {
		TriangleElement<double> currentElement = geom.getElement(index);

//		std::cout << currentElement << std::endl;
//		std::cout << currentElement.getBoundingBox() << std::endl;

		Box<int, 3> ind = gr.getGlobalBoxIndices( currentElement.getBoundingBox() );

//		std::cout << ind << std::endl;


		auto bl = ind.bl();
		auto tr = ind.tr();

		for ( int k = bl[2]; k <= tr[2]; ++k ) {

			if ( k < z || k >= (z+p) )
				continue;

			for ( int j = bl[1]; j <= tr[1]; ++j ) {
				
				if ( j < y || j >= (y+n) )
					continue;

				for ( int i = bl[0]; i <= tr[0]; ++i ) {

					if ( i < x || i >= (x+m) ) 
						continue;

					Eigen::Vector3d        p       = gr.getCoord( Eigen::Vector3i(i,j,k) );

					SignedDistance<double> v       = TriangleElement<double>::computeDistance(currentElement, p);
					double                 pd      = TriangleElement<double>::computePerpendicularDistance(currentElement, p);


					if ( v.dist > narrowBand ) {
						continue;
					}

					if ( v.dist < eps ) {

						data_ptr[k][j][i]             = eps * v.sign;
						perpendicualDistance[k][j][i] = eps;

					}  else if ( fabs( v.dist - fabs(data_ptr[k][j][i]) ) < eps ) {

						// equal not funny
						if ( perpendicualDistance[k][j][i] <= pd ) {
							data_ptr[k][j][i]             = v.sign * v.dist;
							perpendicualDistance[k][j][i] = pd;
						}

					} else if ( v.dist < fabs(data_ptr[k][j][i]) ) {

						// less so ok than
						data_ptr[k][j][i]             = v.sign * v.dist;
						perpendicualDistance[k][j][i] = pd;
					}

				}
			}
		}

	}

}


#endif
