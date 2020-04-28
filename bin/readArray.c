/*
 * readArray.c
 *
 *  Created on: Apr 23, 2013
 *      Author: petr
 */

#include <stdio.h>
#include <stdlib.h>


int readArray(int*, int**, double**, double**, double**, double**, char*);


int main (int argc, char* argv[]) {

	printf("This tests the data exchange file.\n");

	int dim;
	int *M;
	double *dx, *x_lo, *x_hi, *data;

	FILE* fp;

	if ( readArray( &dim, &M, &dx, &x_lo, &x_hi, &data, argv[1]) ) {

		printf("number of dimension is %d\n", dim);
		printf("size of the grid is %d, %d, %d\n", M[0], M[1], M[2]);

		printf("bottom left is %f, %f, %f\n", x_lo[0], x_lo[1], x_lo[2]);
		printf("upper right is %f, %f, %f\n", x_hi[0], x_hi[1], x_hi[2]);

		int i;
		for (i = 0; i < M[0]*M[1]*M[2]; ++i) {
			printf("%f, ", data[i]);
		}
		printf("\n");
	} 

	free(M);
	free(dx);
	free(x_lo);
	free(x_hi);
	free(data);

	return 0;

}

int readArray(int* dim, int** M, double** dx, double** x_lo, double** x_hi, double** data, char* filename) {

	printf("filename = %s\n", filename);

	FILE* fp; 
	fp = fopen(filename, "r");


	if (fp == NULL) {
		return -1;
	}


	fread(dim, sizeof(int), 1, fp);          // read number of dimensions
	
	*M = malloc(sizeof(int) * *dim);
	fread(*M, sizeof(int), *dim, fp);        // read actual size of the grid
	
	*dx = malloc(sizeof(double) * *dim);
	fread(*dx, sizeof(double), *dim, fp);    // read resolution of the grid
	
	*x_lo = malloc(sizeof(double) * *dim);
	fread(*x_lo, sizeof(double), *dim, fp);  // read bottom left corner of the grid
	
	*x_hi = malloc(sizeof(double) * *dim);
	fread(*x_hi, sizeof(double), *dim, fp);  // read top right corner of the grid

	int nCells = 1;
	int i;
	for(i = 0; i < *dim; ++i) {
		nCells *= *(*M+i);
	}


	*data = malloc( nCells * sizeof(double) );

	fread(*data, sizeof(double), nCells, fp); // read the distance values 

	fclose(fp);

	return 1;

}


