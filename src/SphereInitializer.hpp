/*
 * SphereInitializer.hpp
 *
 *  Created on: Nov 8, 2012
 *      Author: petr
 */
#pragma once
#ifndef SPHEREINITIALIZER_HPP_
#define SPHEREINITIALIZER_HPP_

#include <petscdmadda.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Interface.hpp"


class SphereInitializer{

    // where the sphere is located
    Eigen::Vector3d center;

    // how big is the sphere
    PetscReal              radius;

    size_t                 d;

    // geometry is stores here
    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > geometry;
 
    // select which boundaries need to be initialized
    int selectBoundaries(const Box<double, 3>& procAABB);

    // defines the distance to the sphere
    PetscReal sphereEquation(const Eigen::Vector3d& point);

    // initialize the geometry (in this case creates the sphere)
    void initializeGeometry();


public:

    SphereInitializer(const Eigen::Vector3d& center, PetscReal radius);


    // This is where the magic happens and the initializator puts the data in
    template <typename type, int dim>
    void operator() (Interface<type, dim>& interface, bool initAll);



};


///===========================================
///              Implementation
///===========================================

SphereInitializer::SphereInitializer(const Eigen::Vector3d& center, PetscReal radius) {

    this->center = center;
    this->radius = radius;
    this->d = center.size();
    initializeGeometry();

}

template <typename type, int dim>
void SphereInitializer::operator ()(Interface<type, dim>& interface, bool initAll) {


    std::vector<Eigen::Vector3d, Eigen::aligned_allocator<Eigen::Vector3d> > boundaryCells;    // xyz coordinated to compute the distance to

    int                        myBoundaries;     // local boundaries to init

    MPI_Group worldGroup;
    MPI_Comm_group(MPI_COMM_WORLD, &worldGroup);
    // perpendicular distance computed towards the triangle, it is used when we are no longer sure, which normal shall we use

    // initialize the data

    Vec           localData = interface.getLocalData();
    double***     data_ptr;
    int           x, y, z, m, n, p;
    int           vecSize;
    DM            cda;
    Vec           gc;
    DMDACoor3d*** coors;

    VecSet(localData, std::numeric_limits<double>::max());
//    VecSet(localData, -0.5);
    VecGetLocalSize(localData, &vecSize);
    assert( vecSize > 0 );

    int myWorldRank;
    int nProcsInWorld;

    MPI_Comm_rank(MPI_COMM_WORLD, &myWorldRank);
    MPI_Comm_size(MPI_COMM_WORLD, &nProcsInWorld);


    const Grid<type, dim>& gr = interface.getGrid();

    DMDAVecGetArray(gr.getDA(), localData, &data_ptr);
    
    DMDAGetGhostCorners(gr.getDA(), &x, &y, &z, &m, &n, &p);

    DMGetCoordinatesLocal(gr.getDA() ,&gc);
    DMGetCoordinateDM(gr.getDA(), &cda);
    DMDAVecGetArray(cda, gc, &coors);


    if (initAll) {
        // initialize the whole domain
        
        for (int i = x; i < x+m; ++i) {
            for (int j = y; j < y+n; ++j) {
                for (int k = z; k < z+p; ++k) {
                    data_ptr[k][j][i] = 
                        sphereEquation(
                            Eigen::Vector3d(coors[k][j][i].x,
                                            coors[k][j][i].y,
                                            coors[k][j][i].z)                            
                        );

                }
            }
        }

    } else {
        // insert sphere boundaries
    	double    narrowBand = gr.getDx(0) * 3;
        for (int i = x; i < x+m; ++i) {
            for (int j = y; j < y+n; ++j) {
                for (int k = z; k < z+p; ++k) {
                    double dist = 
                        sphereEquation(
                            Eigen::Vector3d(coors[k][j][i].x,
                                            coors[k][j][i].y,
                                            coors[k][j][i].z)                            
                        );
                    if (fabs(dist) > narrowBand) {
                        continue;
                    }
                    data_ptr[k][j][i] = dist;
                }
            }
        }

        // boundaries
        int gh = 1;
        int boundaries = 63;//selectBoundaries( gr.getNodeSpan(myWorldRank) );

        if ( boundaries & 1 ) {
            // this means lets initialize left side
            for (int i = x; i < x+gh; ++i) {
                for (int j = y; j < y+n; ++j) {
                    for (int k = z; k < z+p; ++k) {
                        data_ptr[k][j][i+1] = sphereEquation(
                                                Eigen::Vector3d(coors[k][j][i+1].x,
                                                                coors[k][j][i+1].y,
                                                                coors[k][j][i+1].z)
                                                );
                    }
                }
            }
        }

        if ( boundaries & 2 ) {
            // this means lets initialize right side
            for (int i = x+m-1; i > x+m-1-gh; --i) {
                for (int j = y; j < y+n; ++j) {
                    for (int k = z; k < z+p; ++k) {
                        data_ptr[k][j][i-1] = sphereEquation(
                                                Eigen::Vector3d(coors[k][j][i-1].x,
                                                                coors[k][j][i-1].y,
                                                                coors[k][j][i-1].z)
                                                );
                    }
                }
            }
        }

        if ( boundaries & 4 ) {
            // this means initialize the bottom side
            for (int j = y; j < y+gh; ++j) {
                for (int i = x; i < x+m; ++i) {
                    for (int k = z; k < z+p; ++k) {
                        data_ptr[k][j+1][i] = sphereEquation(
                                                Eigen::Vector3d(coors[k][j+1][i].x,
                                                                coors[k][j+1][i].y,
                                                                coors[k][j+1][i].z)
                                                );
                    }
                }
            }
        }

        if ( boundaries & 8 ) {
            // this means initialize the top side
            for (int j = y+n-1; j > y+n-1-gh; --j) {
                for (int i = x; i < x+m; ++i) {
                    for (int k = z; k < z+p; ++k) {
                        data_ptr[k][j-1][i] = sphereEquation(
                                                Eigen::Vector3d(coors[k][j-1][i].x,
                                                                coors[k][j-1][i].y,
                                                                coors[k][j-1][i].z)
                                                );
                    }
                }
            }
        }

        if ( boundaries & 16 ) {
            // this means initialize the front side
            for (int k = z; k < z+gh; ++k) {
                for (int i = x; i < x+m; ++i) {
                    for (int j = y; j < y+n; ++j) {
                        data_ptr[k+1][j][i] = sphereEquation(
                                                Eigen::Vector3d(coors[k+1][j][i].x,
                                                                coors[k+1][j][i].y,
                                                                coors[k+1][j][i].z)
                                                );
                    }
                }
            }
        }

        if ( boundaries & 32 ) {
            // this means initialize the far side
            for (int k = z+p-1; k > z+p-1-gh; --k) {
                for (int i = x; i < x+m; ++i) {
                    for (int j = y; j < y+n; ++j) {
                        data_ptr[k-1][j][i] = sphereEquation(
                                                Eigen::Vector3d(coors[k-1][j][i].x,
                                                                coors[k-1][j][i].y,
                                                                coors[k-1][j][i].z)
                                                );
                    }
                }
            }
        }

    }

    DMDAVecRestoreArray(cda, gc, &coors);
    DMDAVecRestoreArray(gr.getDA(), localData, &data_ptr);

}

PetscReal SphereInitializer::sphereEquation(const Eigen::Vector3d& point) {
    // this function really speaks for itself. It is just sphere equation

    PetscReal dist = 0;

    
    dist = sqrt( (point[0] - center[0])*(point[0] - center[0]) +
                 (point[1] - center[1])*(point[1] - center[1]) +
                 (point[2] - center[2])*(point[2] - center[2])) - radius;
        

    return dist;

}


int SphereInitializer::selectBoundaries(const Box<double, 3>& procAABB) {
    // chooses which boundaries of the domain have to be initialized
    // ordering : [left, right, bottom, top, front, back]

    int selectedBoundaries_int = 0;

    for (int i = 0; i < 3; ++i) {
        //gr->getLocalGhostedMax(i) <= dataProvider->getMaxX(i)
        if ( procAABB.maxX(i) <= (center[i]+radius) ) {
            selectedBoundaries_int |= ( 1 << (i*2 + 1) );
        }
        // gr->getLocalGhostedMin(i) >= dataProvider->getMinX(i)
        if ( procAABB.minX(i) >= (center[i]-radius) ) {
            selectedBoundaries_int |= ( 1 << (i*2) );
        }
    }

    return selectedBoundaries_int;

}


void SphereInitializer::initializeGeometry() {
    // creates the sphere points in Nd

    int    numberOfPoints = 300;


    // x(t) = center(1) - radius*cos(angle_step_1 * t)*sin(angle_step_2 * t);
    // y(t) = center(2) - radius*sin(angle_step_1 * t)*cos(angle_step_2 * t);
    // z(t) =             radius*                      cos(angle_step_2 * t);

    double angle_step_1 = (2*3.1415926535897) / numberOfPoints;
    double angle_step_2 = (  3.1415926535897) / numberOfPoints;

    geometry.reserve(numberOfPoints*numberOfPoints);

    for (int i = 0; i < numberOfPoints; ++i) {

        for (int j = 0; j < numberOfPoints; ++j) {
                
            double x = center[0] + radius*cos(angle_step_1 * i)*sin(angle_step_2*j);
            double y = center[1] + radius*sin(angle_step_1 * i)*sin(angle_step_2*j);
            double z = center[2] + radius*cos(angle_step_2*j);


            geometry.push_back( Eigen::Vector3d(x, y, z) );
            
        }

    }
    
}













#endif /* SPHEREINITIALIZER_HPP_ */
