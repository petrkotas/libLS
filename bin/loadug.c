/*
 * loadug.c
 *
 *  Created on: Apr 23, 2013
 *      Author: petr
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {

	if (nrhs != 1) {
		mexErrMsgIdAndTxt("LoadUG:load:nrhs", "Filename required.");
	}

	if (nlhs != 2) {
		mexErrMsgIdAndTxt("LoadUG:filename:nlhs", "Output is required.");
	}

	if (mxIsChar(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("LoadUG:load:inputNoString", "Input must be a string.");
	}

	if (mxGetM(prhs[0]) != 1) {
		mexErrMsgIdAndTxt("LoadUG:load:inputNoVector", "Input must be a row vector.");
	}

	/* Load the data and convert it for the Matlab */

	char *fname = mxArrayToString(prhs[0]);

	double*  phi;
	double   bb[6];
	int      mnp[3];
	double   dx;
	long int nGridPoints;

	/**
	 * Load the fucking data
	 */

	const int MAX_LINE_LENGTH = 30;
	char *pEnd;
	char line[MAX_LINE_LENGTH];
	FILE *file;

	file = fopen(fname, "r");


	if( file == NULL ) {
		mexErrMsgIdAndTxt("LoadUG:load:openFileError", "Can't open file.");
	}

    printf("File opened\n");
    
	/* The rest of the code is readig directly a distance function from the file */
	if( fgets(line, MAX_LINE_LENGTH, file) == NULL ) {
		mexErrMsgIdAndTxt("LoadUG:load:readLineError", "File can't be read.");
	}


	while( strncmp((const char*) line, "Gitterweite", 11) != 0 ) {
		if( fgets(line,MAX_LINE_LENGTH,file) == NULL ) {
			fclose(file);
			mexErrMsgIdAndTxt("LoadUG:load:noGitterweiteFound", "Error, parameter Gitterweite/h not found.");
		}
	}


	if( fgets(line, MAX_LINE_LENGTH, file) != NULL ) {
		dx = strtod(line, &pEnd);
        printf("dx = %f\n", dx);
	} else {
		fclose(file);
		mexErrMsgIdAndTxt("LoadUG:load:noGridLength", "Error, grid length was not read correctly.");
	}


	while( strncmp((const char*) line, "BoundingBox", 11) != 0 ) {
		if( fgets(line, MAX_LINE_LENGTH, file) == NULL ) {
			fclose(file);
			mexErrMsgIdAndTxt("LoadUG:load:noBoundingBox", "Error, BoundingBox not found.");
		}
	}


	int j = 0;
	while( strncmp((const char*) line, "Abstand", 7) != 0 ) {
		if( fgets(line, MAX_LINE_LENGTH, file) != NULL ) {
			bb[j] = atof(line+5);
            printf("bb[%i] = %f\n", j, bb[j]);
			j++;
		} else {
			fclose(file);
			mexErrMsgIdAndTxt("LoadUG:load:noBoundingBox", "Error, BoundingBox not found.");
		}
	}


	nGridPoints = atoi(line+8);
    printf("nGridPoints = %li\n", nGridPoints);
/*
	phi = (double*) mxMalloc( ( (nGridPoints) ) * sizeof(double) );

	if( phi == NULL ) {
		mexErrMsgIdAndTxt("LoadUG:load:memoryError", "Error, Memory error in distance.");
	}
*/
    
	long int i = 0;
	/*
    while( (fgets(line, MAX_LINE_LENGTH, file) != NULL) &&  (strncmp((const char*) line, "}", 1) != 0) ) {
		phi[i] = strtod(line, &pEnd);
		if ( i >= (nGridPoints) ) {
			mexErrMsgIdAndTxt("LoadUG:load:wrongNumberOfEntries", "Error, Too much grid points.");
		}
		i++;
	}
    */

	for( i = 0; i < 3; i++ ) {
		mnp[i] = (int) ceil( (bb[2*i] - bb[2*i+1]) / dx) + 1;
        printf("mnp[%i] = %i\n", i, mnp[i]);
	}

	fclose(file);



	/* create outphut matrix */
	plhs[0] = mxCreateNumericArray(3, mnp, mxDOUBLE_CLASS, mxREAL);

	/* populate the real part of the created array */
	double* start_of_pr;
	size_t bytes_to_copy;
/*
	start_of_pr = (double*) mxGetData(plhs[0]);
    bytes_to_copy = nGridPoints * mxGetElementSize(plhs[0]);
    memcpy(start_of_pr, phi, bytes_to_copy);
*/
/*
    int bbmnp[2] = {mnp[0], 3};
    plhs[1] = mxCreateNumericArray(2, bbmnp, mxDOUBLE_CLASS, mxREAL);

    double* data_ptr = (double*) mxGetData(plhs[1]);
    /* compute axis for dim 1 
    for (i = 0; i < mnp[0]; ++i, ++data_ptr) {
		*data_ptr = bb[1] + i*dx;
	}
    /* compute axis for dim 2 
    for (i = 0; i < mnp[1]; ++i, ++data_ptr) {
		*data_ptr = bb[3] + i*dx;
	}
    /* compute axis for dim 3 
    for (i = 0; i < mnp[2]; ++i, ++data_ptr) {
		*data_ptr = bb[5] + i*dx;
	}
*/
}



