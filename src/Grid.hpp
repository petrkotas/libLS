/*
 * Grid.h
 *
 *  Created on: Nov 5, 2012
 *      Author: petr
 */
#pragma once
#ifndef GRID_H_
#define GRID_H_

#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <petscdmda.h>
#include <Eigen/Dense>
#include "BoundingBox.hpp"

/**
 * @brief GridData holds all information about grid distribution.
 * 
 * The information is hold for all nodes. This waste space, however all
 * other manipulation is faster and easier
 * 
 */
struct GridData {

    double minx[3]; ///< global ower left corner of the grid
    double maxx[3]; ///< global upper right corner of the grid

    double lminx[3]; ///< local lower left corner of the grid
    double lmaxx[3]; ///< local upper right corner of the grid
    
    double glminx[3]; ///< local ghosted lower left corner of the grid
    double glmaxx[3]; ///< local ghosted upper right corner of the grid
    
    
    int pm; ///< number of processors m
    int pn; ///< number of processors n
    int pp; ///< number of processors p
    
    
    int m; ///< number of gridcells m
    int n; ///< number of gridcells n
    int p; ///< number of gridcells p

    
    int ghosts; ///< number of ghost cells

    
    int lm; ///< local number of gridcells m
    int ln; ///< local number of gridcells n
    int lp; ///< local number of gridcells p
    
    int glm; ///< local number of ghosted gridcells m
    int gln; ///< local number of ghosted gridcells n
    int glp; ///< local number of ghosted gridcells p


    static void createType(MPI_Datatype* type) {
        const int    nitems = 2;
        int          blocklengths[2] = {18,
                                        13};

        MPI_Datatype types[2] = {MPI_DOUBLE,
                                 MPI_INT};

        MPI_Aint offsets[2], extent;

        offsets[0] = 0;

        MPI_Type_extent(MPI_DOUBLE, &extent);
        offsets[1] = blocklengths[0] * extent;

        MPI_Type_create_struct(nitems, blocklengths, offsets, types, type);
        MPI_Type_commit(type);
    }

    /**
     * @brief Return bounding box of the whole grid
     * @return grid bounding box
     */
    Box<double, 3> getBB() {
        return Box<double, 3>(minx, maxx);
    }
};





/**
 * @brief Grid class wrapping PETSc DMDA grid
 * @tparam type int or double
 * @tparam dim grid dimension
 */
template <typename type, PetscInt dim>
class Grid {

    MPI_Datatype mpi_grid_type;

    GridData* globalGridLayout; ///< Stores information about the global grid distribution

    // grid step size, with appropriate length
    type dx[dim]; ///< deprecated, todo remove

    // no shit Sherlock, yes it is stencil size
    type  stencilSize;

    // lower left corner of the grid
    type minX[dim]; ///< deprecated, todo remove

    // upper right corner of the grid
    type maxX[dim]; ///< deprecated, todo remove

    // store number of grid cells for each dimension - global data
    long int  numberOfGridCells[dim]; ///< deprecated, todo remove

    // store number of ghosted grid cells for each dimension - global data
    long int  numberOfGhostedGridCells[dim]; ///< deprecated, todo remove

    // store the local minimum of grid - ghosted one
    type localMinX[dim]; ///< deprecated, todo remove

    // store the local maximum of the grid - ghosted one
    type localMaxX[dim]; ///< deprecated, todo remove

    // store the local minimum of grid - ghosted one
    type localGhostedMinX[dim]; ///< deprecated, todo remove

    // store the local maximum of the grid - ghosted one
    type localGhostedMaxX[dim]; ///< deprecated, todo remove

    // holds the information about length domains
    type nodeSpan[dim]; ///< deprecated, todo remove

    // number of level sets on the grid
    int dof; ///< deprecated, todo remove

    // processor layout in each direction
    int pm, pn, pp; ///< deprecated, todo remove



    //========================
    // Petsc data structures
    //========================

    // grid associated with this object
    DM da; ///< internally stored DMDA grid

    // distribution and stuff
    DMDALocalInfo localInfo; ///< internally stored DMDA grid info, it is somehow redundant to GridData

    // duplicate info for storing type of the boundary (NONE, GHOSTED, MIRRORED, PERIODIC)
    DMDABoundaryType boundaryType;


    //==========================
    // Private helper function
    //==========================

    /** Initializes the grid class and prepare Petsc data stuctures
    *
    * \param *minX array containing lower point of the grid
    * \param *maxX array containing upper point of the grid
    * \param *M    array with number of grid cells in each direction
    * \param *NP   array with number of slices in each direction
    */
    void init(PetscReal* minX, PetscReal* maxX, PetscInt* M, PetscInt* NP);

    /**
     * \brief Friend ostream operator for debug print
     */
    template <typename T, int d>
    friend std::ostream& operator<< (std::ostream& os, const Grid<T, d>& gr);


public:

    GridData gr_layout;
    int processID, nProcs;


    /// Constructor for old interface, utilizes array pointers
    ///
    /// \param *minX array containing lower point of the grid
    /// \param *maxX array containing upper point of the grid
    /// \param *M    array with number of grid cells in each direction
    /// \param *NP   arrat with number of slices in each direction
    Grid(PetscReal* minX, PetscReal* maxX, PetscInt* M, PetscInt* NP);


    /// Constructor for new interface, utilizes Eigen classes
    ///
    /// \param minX array containing lower point of the grid
    /// \param maxX array containing upper point of the grid
    /// \param M    array with number of grid cells in each direction
    Grid(const Eigen::Array<type, dim, 1>& minX, const Eigen::Array<type, dim, 1>& maxX, const Eigen::Array<int, dim, 1>& M);


    /// Constructor used to create equidistant grids
    ///
    /// \param minX array containing lower point of the grid
    /// \param maxX array containing upper point of the grid
    /// \param M    number of grid cells in first dimension, the rest is computed accoringly to supply equidistant dx
    Grid(const Eigen::Array<type, dim, 1>& minX, const Eigen::Array<type, dim, 1>& maxX, int M);


    /// Constructor used to create equidistant grids
    ///
    /// \param gridSpan bounding box to set the grid span
    /// \param M number of grid cells in each direction
//  Grid(const Box<type, dim>& gridSpan, int M);

    /// Destructor
    ~Grid();

    const DMDALocalInfo* getLocalInfo() const;

    const DM getDA() const {return this->da;};


    /// returns minimum coordinate of the grid for asked dimension number
    PetscReal getMin(PetscInt dimNum) const;

    /// returns maximim coordinate of the grid for asked dimension number
    PetscReal getMax(PetscInt dimNum) const;

    /// return spacing between grid points for asked dimension number
    PetscReal getDx(PetscInt dimNum) const;

    /// return local minimum coordinate
    PetscReal getLocalMin(PetscInt dimNum) const;

    /// return local maximum coordinate
    PetscReal getLocalMax(PetscInt dimNum) const;

    /// returns local ghpsted minimum coordinate of the grid for asked dimension number
    PetscReal getLocalGhostedMin(PetscInt dimNum) const;

    /// returns local ghosted maximim coordinate of the grid for asked dimension number
    PetscReal getLocalGhostedMax(PetscInt dimNum)const;

    /// returns a ghosted bounding box for local region covered by this computational node
    Box<double, 3> getLocalGhostedRegion();


    PetscReal getStencilSize() const;

    /// return real world coordinate for given dimension and axis index
    PetscReal getCoordinate(PetscInt dimNum, PetscInt ind) const;

    /// return linear index from passed matrix indices (linear index is defined in row wise manner)
    long int getLinearIndex(PetscInt* ind) const;

    /// return linear index from passed matrix indices (linear index is defined in row wise manner)
    long int getLinearIndex(std::vector<PetscInt> ind) const;

    /// get coordinate from indices
    inline Eigen::Vector3d getCoord(const Eigen::Vector3i& ind) const;

    /// get coordinate from indices
    const std::vector<PetscReal> getLocalCoord(std::vector<PetscInt> ind) const;

    /// get local coordinate from indices
    Eigen::Vector3d getLocalCoord(PetscInt* ind, int procID) const;

    /// get local coordinates from indices
    const std::vector<PetscReal> getLocalCoord(PetscInt* ind) const;

    long int getNumberOfCells() const;

    long int getNumberOfLocalCells() const;

    PetscInt getM(PetscInt dimNum) const;

    PetscInt getLocalM(PetscInt dimNum) const;

    PetscReal getMaxDist() const;

    /// method that transform coordinate values to indices in the grid. Not sure if this is done in PETSc and
    Box<PetscInt, dim>& getPointSpan(const std::vector<PetscReal>& coord, Box<PetscInt, dim>& prealloc);

    /// trasfoms coordinate box into indices box
    Box<PetscInt, dim> getGlobalBoxIndices(const Box<type, dim>& bb) const;

    Box<PetscInt, 3> locateBoxIndices(const Box<PetscReal, 3>& bb);

    int getDataPosition(const type point[dim]);

    std::list<int> getProcessSpan(const Box<double, dim>& indBox) const;

    Box<double, dim> getNodeSpan(int indices[dim]) const;

    Box<double, dim> getNodeSpan() const;

    Box<double, dim> getNodeSpan(int procID) const;

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > getNodeBoundaryCells(int procID, int boundaries) const;
};





///==================================
///         Implementation
///==================================


template <typename type, int dim>
void Grid<type, dim>::init(PetscReal* minX, PetscReal* maxX, PetscInt* M, PetscInt* NP) {

    GridData::createType(&mpi_grid_type);

    PetscErrorCode ierr;

    this->boundaryType = DMDA_BOUNDARY_GHOSTED;

    // values setup for dim > N are ommited
    for (int i = 0; i < dim; ++i) {
        this->numberOfGridCells[i] = M[i];
        this->dx[i]                = (maxX[i] - minX[i])/ (M[i] - 1);
        this->minX[i]              = minX[i];
        this->maxX[i]              = maxX[i];
    }

    this->stencilSize = 1;

    //  this WILL alwasy be 1 because the data is stored in other class
    //  and changing this would corrupt it all. The petsc way is interleave the data
    //  this is not exactly what we are going to do.
    this->dof = 1;

    switch (dim) {
        case 1:
            break;
        case 2:
            break;
        case 3:
            // To make petsc work correctly, we need to convert number of gridcells to number of nodes => M[.] = M[.] + 1
            ierr = DMDACreate3d(PETSC_COMM_WORLD, boundaryType, boundaryType, boundaryType, DMDA_STENCIL_BOX,
                                    M[0], M[1], M[2], NP[0], NP[1], NP[2], dof, stencilSize,
                                    PETSC_NULL, PETSC_NULL, PETSC_NULL, &da);
            ierr = DMDASetUniformCoordinates(da, minX[0], maxX[0], minX[1], maxX[1], minX[2], maxX[2]);

            // very stupid way to get back the process layout
            ierr = DMDAGetInfo(da,
                               PETSC_NULL,
                               PETSC_NULL, PETSC_NULL, PETSC_NULL,
                               &pm, &pn, &pp,
                               PETSC_NULL,
                               PETSC_NULL,
                               PETSC_NULL, PETSC_NULL, PETSC_NULL,
                               PETSC_NULL);

            nodeSpan[0] = (maxX[0] - minX[0]) / (pm);
            nodeSpan[1] = (maxX[1] - minX[1]) / (pn);
            nodeSpan[2] = (maxX[2] - minX[2]) / (pp);

            break;
        default:
            break;
    }

    if (PetscUnlikely(ierr)) {
        // some error handling add lated here :)
    }

    ierr = DMDAGetLocalInfo(da, &localInfo);

    MPI_Comm_rank(PETSC_COMM_WORLD, &processID);
    MPI_Comm_size(PETSC_COMM_WORLD, &nProcs);




    // for convenience this part defines local chunks of grid, so user do not need to recompute
    // it every time he needs it
    switch (dim) {
        case 1:
            break;
        case 2:
            break;
        case 3:
            localMinX[0]        = minX[0] +  localInfo.xs * dx[0];
            localMaxX[0]        = minX[0] + (localInfo.xs + localInfo.xm - 1) * dx[0];
            localGhostedMinX[0] = minX[0] +  localInfo.gxs * dx[0];
            localGhostedMaxX[0] = minX[0] + (localInfo.gxs + localInfo.gxm - 1) * dx[0];

            localMinX[1]        = minX[1] +  localInfo.ys * dx[1];
            localMaxX[1]        = minX[1] + (localInfo.ys + localInfo.ym - 1) * dx[1];
            localGhostedMinX[1] = minX[1] +  localInfo.gys * dx[1];
            localGhostedMaxX[1] = minX[1] + (localInfo.gys + localInfo.gym - 1) * dx[1];

            localMinX[2]        = minX[2] +  localInfo.zs * dx[2];
            localMaxX[2]        = minX[2] + (localInfo.zs + localInfo.zm - 1) * dx[2];
            localGhostedMinX[2] = minX[2] +  localInfo.gzs * dx[2];
            localGhostedMaxX[2] = minX[2] + (localInfo.gzs + localInfo.gzm -1) * dx[2];

            numberOfGhostedGridCells[0] = localInfo.gxm;
            numberOfGhostedGridCells[1] = localInfo.gym;
            numberOfGhostedGridCells[2] = localInfo.gzm;

            break;
    }


    for (int i = 0; i < 3; ++i) {

        gr_layout.lminx[i] = localMinX[i];
        gr_layout.lmaxx[i] = localMaxX[i];
        gr_layout.glminx[i] = localGhostedMinX[i];
        gr_layout.glmaxx[i] = localGhostedMaxX[i];

    }

    gr_layout.glm = localInfo.gxm;
    gr_layout.gln = localInfo.gym;
    gr_layout.glp = localInfo.gzm;

    gr_layout.lm  = localInfo.mx;
    gr_layout.ln  = localInfo.my;
    gr_layout.lp  = localInfo.mz;

    gr_layout.m   = M[0];
    gr_layout.n   = M[1];
    gr_layout.p   = M[2];

    gr_layout.minx[0] = minX[0];
    gr_layout.minx[1] = minX[1];
    gr_layout.minx[2] = minX[2];

    gr_layout.maxx[0] = maxX[0];
    gr_layout.maxx[1] = maxX[1];
    gr_layout.maxx[2] = maxX[2];

    gr_layout.ghosts = 1;


    globalGridLayout = new GridData[nProcs]; // do not forget to delete this little shit


    MPI_Allgather(&gr_layout, 1, mpi_grid_type, globalGridLayout, 1, mpi_grid_type, MPI_COMM_WORLD);

    std::stringstream ss;
    
    ss << "====================" << std::endl;
    ss << "procid : " << processID << std::endl;
    ss << "min: " << gr_layout.glminx[0] << ", " << gr_layout.glminx[1] << ", " << gr_layout.glminx[2] << std::endl;
    ss << "max: " << gr_layout.glmaxx[0] << ", " << gr_layout.glmaxx[1] << ", " << gr_layout.glmaxx[2] << std::endl;
    ss << "====================" << std::endl;
     for (int i = 0; i < nProcs; ++i) {
        ss << "global proc id : " << i << std::endl;
        ss << "min: " << globalGridLayout[i].glminx[0] << ", " << globalGridLayout[i].glminx[1] << ", " << globalGridLayout[i].glminx[2] << std::endl;
        ss << "max: " << globalGridLayout[i].glmaxx[0] << ", " << globalGridLayout[i].glmaxx[1] << ", " << globalGridLayout[i].glmaxx[2] << std::endl;
     }
    ss << "====================" << std::endl;

    std::string s = ss.str();
//    PetscSynchronizedPrintf(PETSC_COMM_WORLD, s.c_str());
 //   PetscSynchronizedFlush(PETSC_COMM_WORLD);
    // create indexing for the data so, because we will need fast acces for inicialization computations




}

template <typename type, int dim>
Grid<type, dim>::Grid(const Eigen::Array<type, dim, 1>& minX, const Eigen::Array<type, dim, 1>& maxX, const Eigen::Array<int, dim, 1>& M) {

    int            m[3] = {M[0], M[1], M[2]},
                  np[3] = {0, 0, 0};
    type         min[3] = {minX[0], minX[1], minX[2]},
                 max[3] = {maxX[0], maxX[1], maxX[2]};

    init(min, max, m, np);

}

template <typename type, int dim>
Grid<type, dim>::Grid(const Eigen::Array<type, dim, 1>& minX, const Eigen::Array<type, dim, 1>& maxX, int M) {

    Eigen::Array<type, dim, 1> len   = maxX - minX;
    Eigen::Array<type, dim, 1> maxX_;   // new maximum size of the grid, possibly expanded a bit
    Eigen::Array<type, dim, 1> M_float; // stores array with accoring number of cells

    M_float[0] = M;

    type dx = len[0] / (M - 1); // size of one cell on dim = 0

//  std::cout << dx << std::endl;

    M_float[1] = std::ceil( len[1] / dx + 1);
    M_float[2] = std::ceil( len[2] / dx + 1);

    Eigen::Array<type, dim, 1> newLen = dx * (M_float - 1);
    Eigen::Array<type, dim, 1> diff   = newLen - len;

    maxX_ = maxX + diff;

    Eigen::Array<int, dim, 1> M_int = M_float.template cast<int>();

//  std::cout << M_int << std::endl;

    int  m[3]   = {M_int[0], M_int[1], M_int[2]},
         np[3]  = {0, 0, 0};
    type min[3] = {minX[0], minX[1], minX[2]},
         max[3] = {maxX_[0], maxX_[1], maxX_[2]};

    init(min, max, m, np);

}

template <typename type, int dim>
Grid<type, dim>::Grid(PetscReal* minX, PetscReal* maxX, PetscInt* M, PetscInt* NP) {
    // this constructor version has one huge problem. There is now way of checking if user passed
    // the inicialized values or just some random stuff.

    init(minX, maxX, M, NP);

}



template <typename type, int dim>
Grid<type, dim>::~Grid() {
    // definitely is nice to clean
//  std::cout << "call destructor in grid" << std::endl;
    PetscBool finalized;
    PetscFinalized(&finalized);

    delete [] globalGridLayout;

    if (!finalized) {
        DMDestroy(&da);
    } // else I am not really sure what happen, will the object stay hanging?
    else {
        std::cout << "already finalized" << std::endl;
    }

}



// this should be removed !!!!!!!!!!!!!
template <typename type, int dim>
const DMDALocalInfo* Grid<type, dim>::getLocalInfo() const{

    return &localInfo;

}



//=====================//
//      Getters        //
//=====================//

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getMin(PetscInt dimNum) const {
    return minX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getMax(PetscInt dimNum) const {
    return maxX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getDx(PetscInt dimNum) const {
    return dx[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getLocalMin(PetscInt dimNum) const {
    return localMinX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getLocalMax(PetscInt dimNum) const {
    return localMaxX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getLocalGhostedMin(PetscInt dimNum) const {
    return localGhostedMinX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getLocalGhostedMax(PetscInt dimNum) const {
    return localGhostedMaxX[dimNum];
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getStencilSize() const {
    return stencilSize;
}

template <typename type, int dim>
inline PetscReal Grid<type, dim>::getCoordinate(PetscInt dimNum, PetscInt ind) const {
    return localMinX[dimNum] + dx[dimNum]*ind;
}

template <typename type, int dim>
inline long int Grid<type, dim>::getLinearIndex(PetscInt* ind) const {
    if (dim == 2) {
        return ind[1]*localInfo.gxm + ind[0];
    } else if (dim == 3) {
        return ind[2]*localInfo.gxm*localInfo.gym + ind[1]*localInfo.gxm + ind[0];
    }
    return -1;
}

template <typename type, int dim>
inline long int Grid<type, dim>::getLinearIndex(std::vector<PetscInt> ind) const {
    if (dim == 2) {
        return ind[1]*localInfo.gxm + ind[0];
    } else if (dim == 3) {
        return ind[2]*localInfo.gxm*localInfo.gym + ind[1]*localInfo.gxm + ind[0];
    }
    return -1;
}

template <typename type, int dim>
inline Eigen::Vector3d Grid<type, dim>::getCoord(const Eigen::Vector3i& ind) const {
    Eigen::Vector3d x_lo(minX), x_dx(dx);

    return x_lo + x_dx.cwiseProduct( ind.cast<double>() );
}


template <typename type, int dim>
inline const std::vector<PetscReal> Grid<type, dim>::getLocalCoord(PetscInt* ind) const {
    std::vector<PetscReal> coord;
    for (int i = 0; i < dim; ++i) {
        coord.push_back( localGhostedMinX[i] + ind[i]*dx[i] );
    }

    return coord;
}

template <typename type, int dim>
inline Eigen::Vector3d Grid<type, dim>::getLocalCoord(PetscInt* ind, int procID) const {
    Eigen::Vector3d coord;
    for (int i = 0; i < dim; ++i) {
        coord[i] = globalGridLayout[procID].glminx[i] + ind[i]*dx[i];
    }

    return coord;
}

template <typename type, int dim>
inline const std::vector<PetscReal> Grid<type, dim>::getLocalCoord(std::vector<PetscInt> ind) const {
    std::vector<PetscReal> coord;
    for (int i = 0; i < ind.size(); ++i) {
        coord.push_back( localGhostedMinX[i] + ind[i]*dx[i] );
    }

    return coord;
}

template <typename type, int dim>
inline long int Grid<type, dim>::getNumberOfCells() const {
    long int numberOfCells = 1;
    for (int i = 0; i < dim; ++i) {
        numberOfCells *= numberOfGridCells[i];
    }
    return numberOfCells;
}

template <typename type, int dim>
inline long int Grid<type, dim>::getNumberOfLocalCells() const {
    long int numberOfCells = 1;
    if (dim >= 1) {
        numberOfCells *= localInfo.gxm;
    }
    if (dim >= 2) {
        numberOfCells *= localInfo.gym;
    }
    if (dim == 3) {
        numberOfCells *= localInfo.gzm;
    }

    return numberOfCells;
}


template <typename type, int dim>
inline PetscInt Grid<type, dim>::getM(PetscInt dimNum) const {
    return numberOfGridCells[dimNum];
}


template <typename type, int dim>
inline PetscInt Grid<type, dim>::getLocalM(PetscInt dimNum) const {
    switch (dimNum) {
        case 0:
            return localInfo.gxm;
            break;
        case 1:
            return localInfo.gym;
            break;
        case 2:
            return localInfo.gzm;
            break;
        default:
            return -1;
            break;
    }
}


template <typename type, int dim>
inline Box<double, 3> Grid<type, dim>::getLocalGhostedRegion() {

    return Box<double, 3>(localGhostedMinX, localGhostedMaxX);

}


template <typename type, int dim>
inline PetscReal Grid<type, dim>::getMaxDist() const {
//  return sqrt( (localGhostedMaxX[0] - localGhostedMinX[0])*(localGhostedMaxX[0] - localGhostedMinX[0]) +
//               (localGhostedMaxX[1] - localGhostedMinX[1])*(localGhostedMaxX[1] - localGhostedMinX[1]) +
//               (localGhostedMaxX[2] - localGhostedMinX[2])*(localGhostedMaxX[2] - localGhostedMinX[2])   );
    return     ( (maxX[0] - minX[0])*(maxX[0] - minX[0]) +
                 (maxX[1] - minX[1])*(maxX[1] - minX[1]) +
                 (maxX[2] - minX[2])*(maxX[2] - minX[2])   );

}



template <typename type, int dim>
Box<PetscInt, dim>& Grid<type, dim>::getPointSpan(const std::vector<PetscReal>& coord, Box<PetscInt, dim>& prealloc) {
// method that transform coordinate values to indices in the grid. Not sure if this is done in PETSc and
// I have no time to go through the doucmentation.
// TODO: resolve this in future :)

    int       bandSize = 3;
    int       indRow[dim];

    for (PetscInt i = 0; i < dim; ++i) {

        PetscInt ind = round( (coord[i] - localGhostedMinX[i]) / dx[i] );
        indRow[i] = ind;

    }


    if (dim == 1) {
// pass
    } if (dim == 2) {
// pass
    } if (dim == 3) {

        int k_low  = (indRow[2] - bandSize) <  0         ? 0                       : (indRow[2] - bandSize);
        int k_high = (indRow[2] + bandSize) >= localInfo.gzm ? (localInfo.gzm - 1) : (indRow[2] + bandSize);

        int j_low  = (indRow[1] - bandSize) <  0         ? 0                       : (indRow[1] - bandSize);
        int j_high = (indRow[1] + bandSize) >= localInfo.gym ? (localInfo.gym - 1) : (indRow[1] + bandSize);

        int i_low  = (indRow[0] - bandSize) <  0         ? 0                       : (indRow[0] - bandSize);
        int i_high = (indRow[0] + bandSize) >= localInfo.gxm ? (localInfo.gxm - 1) : (indRow[0] + bandSize);


        int minx[3] = {i_low, j_low, k_low};
        int maxx[3] = {i_high, j_high, k_high};

        prealloc.updateFromMinMax(minx, maxx);

    }

    return prealloc;

}

template <typename type, int dim>
Box<PetscInt, dim> Grid<type, dim>::getGlobalBoxIndices(const Box<type, dim>& bb) const {
    
    PetscInt padding = 2;
    PetscInt gMin[dim], gMax[dim];

    for (int i = 0; i < dim; ++i) {

        gMin[i] = floor( (bb.minX(i) - this->minX[i]) / dx[i] ) - padding;
        gMax[i] = ceil(  (bb.maxX(i) - this->minX[i]) / dx[i] ) + padding;

        gMin[i] = (gMin[i] < 0)                     ? 0                        : gMin[i];
        gMax[i] = (gMax[i] >= numberOfGridCells[i]) ? (numberOfGridCells[i]-1) : gMax[i];

    }

//  std::cout << "gmin : " << gMin[0] << ", " << gMin[1] << ", " << gMin[2] << std::endl;
//  std::cout << "gmax : " << gMax[0] << ", " << gMax[1] << ", " << gMax[2] << std::endl;

    return Box<PetscInt, dim>(gMin, gMax);

}

template <typename type, int dim>
Box<PetscInt, 3> Grid<type, dim>::locateBoxIndices(const Box<PetscReal, 3>& bb) {
    // method that transform coordinate values to indices in the grid. Not sure if this is done in PETSc and
    // I have no time to go through the documentation.
    // TODO: resolve this in future :)

    PetscInt padding = 2;
    PetscInt minx[dim], maxx[dim];

    for (int i = 0; i < dim; ++i) {

        minx[i] = floor( (bb.minX(i) - localGhostedMinX[i]) / dx[i] ) - padding;
        maxx[i] = ceil(  (bb.maxX(i) - localGhostedMinX[i]) / dx[i] ) + padding;

        minx[i] = (minx[i] < 0)                            ? 0                               : minx[i];
        maxx[i] = (maxx[i] >= numberOfGhostedGridCells[i]) ? (numberOfGhostedGridCells[i]-1) : maxx[i];


    }

    return Box<PetscInt, dim>(minx, maxx);

}


// Returns physical position of the data in question
template <typename type, int dim>
int Grid<type, dim>::getDataPosition(const type point[dim]) {
    Eigen::Array<type, dim, 1> x_lo(minX), p(point), span(nodeSpan), x_hi(maxX);
    Eigen::Array<int, dim, 1> procs;

    procs << 1, pm, pm*pn;

    auto pos_cor = ((p - x_lo) / span).template cast<int>();

    return (pos_cor*procs).sum();

}

template <typename type, int dim>
std::list<int> Grid<type, dim>::getProcessSpan(const Box<double, dim>& indBox) const {
    Eigen::Array<type, dim, 1> x_lo(minX), span(nodeSpan), x_hi(maxX);
    Eigen::Array<int, dim, 1> procs;
    std::list<int> processes;
    double eps = 1E-10; // this is purely local epsilon need to be machine eps

    procs << 1, pm, pm*pn;

    auto pos_cor_bl = ((indBox.bl() - x_lo - eps) / span).template cast<int>();
    auto pos_cor_tr = ((indBox.tr() - x_lo - eps) / span).template cast<int>();

    for (int k = pos_cor_bl[2]; k <= pos_cor_tr[2]; ++k) {
        for (int j = pos_cor_bl[1]; j <= pos_cor_tr[1]; ++j) {
            for (int i = pos_cor_bl[0]; i <= pos_cor_tr[0]; ++i) {
                processes.push_back( (Eigen::Array3i(i,j,k)*procs).sum() );
            }
        }
    }

    return processes;

}


template <typename type, int dim>
Box<double, dim> Grid<type, dim>::getNodeSpan(int indices[dim]) const {

    Eigen::Array<type, dim, 1> x_lo(minX), span(nodeSpan);
    Eigen::Array<int, dim, 1> inds(indices);


    Eigen::Array<double, dim, 1> c = x_lo + (0.5*span)*(inds+1).template cast<double>();
    Eigen::Array<double, dim, 1> e = 0.5*span;


    return Box<double, 3>( c.data(), e.data() );

}

template <typename type, int dim>
Box<double, dim> Grid<type, dim>::getNodeSpan() const {

    return Box<double, 3>(globalGridLayout[processID].glminx, globalGridLayout[processID].glmaxx);

}

template <typename type, int dim>
Box<double, dim> Grid<type, dim>::getNodeSpan(int procID) const {
    std::stringstream ss;
    ss << "ProcID = " << procID << std::endl;
    ss << globalGridLayout[procID].glminx[0] << ", " << globalGridLayout[procID].glminx[1] << ", " << globalGridLayout[procID].glminx[2] << std::endl;
    ss << globalGridLayout[procID].glmaxx[0] << ", " << globalGridLayout[procID].glmaxx[1] << ", " << globalGridLayout[procID].glmaxx[2] << std::endl;
    std::string s = ss.str();
    PetscSynchronizedPrintf(PETSC_COMM_WORLD, s.c_str());
    PetscSynchronizedFlush(PETSC_COMM_WORLD);
    return Box<double, 3>(globalGridLayout[procID].glminx, globalGridLayout[procID].glmaxx);

}

template <typename type, int dim>
std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > Grid<type, dim>::getNodeBoundaryCells(int procID, int boundaries) const {

    // std::cout << "proc: " << procID << " , boundaries : " << boundaries << std::endl;

    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > cellList;
    long int reserveSize = (globalGridLayout[procID].glm + globalGridLayout[procID].gln + globalGridLayout[procID].glp) / 3;
    reserveSize          = reserveSize*reserveSize*4;

    cellList.reserve( reserveSize );

    if ( boundaries & 1 ) {
        // this means lets initialize left side
        for (int i = 0; i < globalGridLayout[procID].ghosts; ++i) {
            for (int j = 0; j < globalGridLayout[procID].gln; ++j) {
                for (int k = 0; k < globalGridLayout[procID].glp; ++k) {

                    PetscInt inds[dim]    = {i, j, k};
                    Eigen::Vector3d coord  = getLocalCoord(inds, procID);
                    cellList.push_back(coord);

                }
            }
        }
    }
    if ( boundaries & 2 ) {
        // this means lets initialize right side
        for (int i = 0; i < globalGridLayout[procID].ghosts; ++i) {
            for (int j = 0; j < globalGridLayout[procID].gln; ++j) {
                for (int k = 0; k < globalGridLayout[procID].glp; ++k) {

                    PetscInt inds[dim]    = {( globalGridLayout[procID].glm - 1) - i, j, k};
                    Eigen::Vector3d coord  = getLocalCoord(inds, procID);
                    cellList.push_back(coord);

                }
            }
        }
    }
    if ( boundaries & 4 ) {
        // this means initialize the bottom side
        for (int j = 0; j < globalGridLayout[procID].ghosts; ++j) {
            for (int i = 0; i < globalGridLayout[procID].glm; ++i) {
                for (int k = 0; k < globalGridLayout[procID].glp; ++k) {

                    PetscInt inds[dim]           = {i, j, k};
                    Eigen::Vector3d coord = getLocalCoord(inds, procID);
                    cellList.push_back(coord);

                }
            }
        }
    }
    if ( boundaries & 8 ) {
        // this means initialize the top side
        for (int j = 0; j < globalGridLayout[procID].ghosts; ++j) {
            for (int i = 0; i < globalGridLayout[procID].glm; ++i) {
                for (int k = 0; k < globalGridLayout[procID].glp; ++k) {

                    PetscInt inds[dim]           = {i, (globalGridLayout[procID].gln-1) - j, k};
                    Eigen::Vector3d coord = getLocalCoord(inds, procID);
                    cellList.push_back(coord);

                }
            }
        }
    }
    if ( boundaries & 16 ) {
        // this means initialize the front side
        for (int k = 0; k < globalGridLayout[procID].ghosts; ++k) {
            for (int i = 0; i < globalGridLayout[procID].glm; ++i) {
                for (int j = 0; j < globalGridLayout[procID].gln; ++j) {

                    PetscInt inds[dim]   = {i, j, k};
                    Eigen::Vector3d coord = getLocalCoord(inds, procID);
                    cellList.push_back(coord);



                }
            }
        }
    }
    if ( boundaries & 32 ) {
        // this means initialize the far side
        for (int k = 0; k < globalGridLayout[procID].ghosts; ++k) {
            for (int i = 0; i < globalGridLayout[procID].glm; ++i) {
                for (int j = 0; j < globalGridLayout[procID].gln; ++j) {

                    PetscInt inds[dim]           = {i, j, (globalGridLayout[procID].glp-1) - k};
                    Eigen::Vector3d coord = getLocalCoord(inds, procID);
                    cellList.push_back(coord);

                }
            }
        }
    }

    return cellList;

}
















///////////////////////////////////////
/// Ostream operator for debug print
///////////////////////////////////////




template <typename T, int d>
std::ostream& operator<< (std::ostream& os, const Grid<T, d>& gr) {

    os << "Grid info for process: " << gr.processID << std::endl;
    os << "-------------------------------------------------------------------------------------" << std::endl;
    if (gr.processID == 0) {
        os << "Global grid size: " << "X: " << gr.minX[0] << " -- " << gr.maxX[0] << std::endl <<
              "                  " << "Y: " << gr.minX[1] << " -- " << gr.maxX[1] << std::endl <<
              "                  " << "Z: " << gr.minX[2] << " -- " << gr.maxX[2] << std::endl <<
              "Global number of grid points: " << "[X, Y, Z] : [" << gr.numberOfGridCells[0] << ", " << gr.numberOfGridCells[1] << ", " << gr.numberOfGridCells[2] << "]" << std::endl <<
              "Resolution for the grid is: " << "[dx, dy, dz] : [" << gr.dx[0] << ", " << gr.dx[1] << ", " << gr.dx[2] << "]" << std::endl;
        os << "-------------------------------------------------------------------------------------" << std::endl;
    }
    os <<  "Local grid size: " << "X: " << gr.localMinX[0] << " -- " << gr.localMaxX[0] << std::endl <<
           "                 " << "Y: " << gr.localMinX[1] << " -- " << gr.localMaxX[1] << std::endl <<
           "                 " << "Z: " << gr.localMinX[2] << " -- " << gr.localMaxX[2] << std::endl <<
           "Local ghosted grid size: " << "X: " << gr.localGhostedMinX[0] << " -- " << gr.localGhostedMaxX[0] << std::endl <<
           "                         " << "Y: " << gr.localGhostedMinX[1] << " -- " << gr.localGhostedMaxX[1] << std::endl <<
           "                         " << "Z: " << gr.localGhostedMinX[2] << " -- " << gr.localGhostedMaxX[2] << std::endl <<
           "Local grid numbering: " << "X: " << gr.localInfo.xs << " -- " << gr.localInfo.xs + gr.localInfo.xm - 1 << std::endl <<
           "                      " << "Y: " << gr.localInfo.ys << " -- " << gr.localInfo.ys + gr.localInfo.ym - 1 << std::endl <<
           "                      " << "Z: " << gr.localInfo.zs << " -- " << gr.localInfo.zs + gr.localInfo.zm - 1 << std::endl <<
           "Local ghosted grid numbering: " << "X: " << gr.localInfo.gxs << " -- " << gr.localInfo.gxs + gr.localInfo.gxm - 1 << std::endl <<
           "                              " << "Y: " << gr.localInfo.gys << " -- " << gr.localInfo.gys + gr.localInfo.gym - 1 << std::endl <<
           "                              " << "Z: " << gr.localInfo.gzs << " -- " << gr.localInfo.gzs + gr.localInfo.gzm - 1 << std::endl;
    os << "-------------------------------------------------------------------------------------" << std::endl;

    return os;

}







#endif /* GRID_H_ */
