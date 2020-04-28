static char help[] = "main file launching the whole computation - kind of setup script";

#include <petscdmda.h>
#include <petsctime.h>
#include <petsclog.h>

#include <vector>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <Eigen/Dense>
#include <cassert>

#include "Geometry.hpp"
#include "Grid.hpp"
#include "Initializer.hpp"
#include "Interface.hpp"
#include "FmmWrapper.hpp"
#include "TriangleElement.hpp"
#include "tictoc.hpp"
#include "BINWritter.hpp"



#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv) {

	PetscErrorCode ierr;
	ierr = PetscInitialize(&argc, &argv, (char*)0, help); CHKERRQ(ierr);


    int init_event, loadData_event, solve_event;

    bool write_out = true;

    Eigen::Array3i M;            // number of nodes in each dimension

    

    char           fname[120];
    char           oname[120];
    PetscBool      flg;
    PetscLogDouble t0, t1, t2, t3, t4, t5, fmm_time, init_time, load_time;
    int rank, size;


    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    ///////////////////
    /// GET OPTIONS ///
    ///////////////////

    PetscOptionsGetString(PETSC_NULL, "-f", fname, 120, &flg);
    if (!flg) {
        // no input filename
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "C'mon gimme something!\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscFinalize();
        return 1;
    }

    PetscOptionsGetString(PETSC_NULL, "-o", oname, 120, &flg);
    if (!flg) {
        // no output filename default will be used
        write_out = false;
    }

    int tmpM, tmpN, tmpP;
    int M[3], NP[3];

    PetscOptionsGetInt(PETSC_NULL,"-m", &M[0], &flg);
    if (!flg) {
        // no input size
        PetscSynchronizedPrintf(PETSC_COMM_WORLD, "At least one dimension must be given!\n");
        PetscSynchronizedFlush(PETSC_COMM_WORLD);
        ierr = PetscFinalize();
        return 1;
    }

    PetscOptionsGetInt(PETSC_NULL,"-n", &M[1], &flg);
    if (!flg) {
        M[1] = M[0];
        // no worry, everything is ok
    }

    PetscOptionsGetInt(PETSC_NULL,"-p", &M[2], &flg);
    if (!flg) {
        M[2] = M[0];
        // no worry, everything is ok
    }

    PetscOptionsGetInt(PETSC_NULL,"-nm", &NP[0], &flg);
    if (!flg) {
        NP[0] = PETSC_DECIDE;
        // no worry, everything is ok
    } 

    PetscOptionsGetInt(PETSC_NULL,"-nn", &NP[1], &flg);
    if (!flg) {
        NP[1] = PETSC_DECIDE;
        // no worry, everything is ok
    } 

    PetscOptionsGetInt(PETSC_NULL,"-np", &NP[2], &flg);
    if (!flg) {
        NP[2] = PETSC_DECIDE;
        // no worry, everything is ok
    }   


    int numberOfGroups = 1;
    PetscOptionsGetInt(PETSC_NULL,"-g", &numberOfGroups, &flg);
    if (!flg) {
        // no worry, everything is ok
        numberOfGroups = 1;
    }

    double growCoef = 1.01;
    PetscOptionsGetReal(PETSC_NULL,"-grow", &growCoef, &flg);
    if (!flg) {
        // no worry, everything is ok
        growCoef = 1.01;
    }

    int initAll = 0;
	PetscOptionsGetInt(PETSC_NULL,"-all", &initAll, &flg);
	if (!flg) {
		// no worry, everything is ok
		initAll = 0;
	}

    ////////////////////////////
    /// END SETUP PARAMETERS ///
    ////////////////////////////

    int procsPerGroup  = ceil( double(size) / double(numberOfGroups) );
    int groupID        = rank / procsPerGroup;

    Eigen::Array3d min, max;

    min << -0.5, -0.5, -0.5;
    max <<  0.5,  0.5,  0.5;


    PetscLogEventRegister("loadData", 0, &loadData_event);
    PetscLogEventBegin(loadData_event, 0, 0, 0, 0);

    Geometry geom = Geometry::fromFile(fname, growCoef);
    geom.initSearchAccelerator();

    PetscLogEventEnd(loadData_event, 0, 0, 0, 0);


    // SphereInitializer init( Eigen::Vector3d(0,0,0), 0.25 );
    Initializer init;
    Grid<double, 3>  gr(geom.aabb.bl(), geom.aabb.tr(), M, NP);
    // Grid<double, 3>  gr(min, max, tmpM);


    geom.preSortElements(gr.getNodeSpan(), gr.getDx(0), numberOfGroups, groupID);
    // std::cout << "presort done" << std::endl;


    Interface<double, 3>* interface = new Interface<double, 3>(gr);
    

    PetscLogEventRegister("initData", 0, &init_event);
	PetscLogEventBegin(init_event, 0, 0, 0, 0);

    init(geom, *interface, groupID, numberOfGroups);
    // init(*interface, initAll);

    PetscLogEventEnd(init_event, 0, 0, 0, 0);


    double shift = gr.getMaxDist();

    PetscLogEventRegister("solve", 0, &solve_event);
	PetscLogEventBegin(solve_event, 0, 0, 0, 0);

	if (!initAll) {
		std::cout << "FMM solver" << std::endl;
		FmmWrapper::solveEikonalEquation<3>(gr, interface, shift);
	}

	PetscLogEventEnd(solve_event, 0, 0, 0, 0);

    if (write_out) {
        // std::cout << "Writting out" << std::endl;
        BINWritter::writeData<double,3>(*interface, oname);
        // std::cout << "write" << std::endl;
    }


//  delete interface;
//  delete gr;


    ierr = PetscFinalize();
    return 0;

}
