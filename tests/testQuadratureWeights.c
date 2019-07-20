/*
 * testQuadratureWeights.c
 *
 *  Created on: Jan 16, 2019
 *      Author: nathanvaughn
 */


#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



/* initializing array for Clenshaw-Curtis weights */



int main()
{
	printf("Entering main() in tests...\n\n");

	// Set the global variables
	torder = 4;
	torderlim = torder+1;


	// Generate the Global Variable unscaledQuadratureWeights.  Nodes use this to compute their scaled quadrature weights
	double CCtheta, b;
	double *unscaledQuadratureWeights;
	int i,j;

	make_vector(unscaledQuadratureWeights, torderlim);

	for ( i = 0; i < torderlim; i++ ){
		  CCtheta = ( double ) ( i ) * M_PI / ( double ) ( torderlim - 1 );

		  unscaledQuadratureWeights[i] = 1.0;

		  for ( j = 1; j <= ( torderlim - 1 ) / 2; j++ )
		  {
			if ( 2 * j == ( torderlim - 1 ) ){
			  b = 1.0;
			}

			else{
			  b = 2.0;
			}

			unscaledQuadratureWeights[i] = unscaledQuadratureWeights[i] - b * cos ( 2.0 * ( double ) ( j ) * CCtheta )/ ( double ) ( 4 * j * j - 1 );
		  }
		}

	unscaledQuadratureWeights[0] = unscaledQuadratureWeights[0] / ( double ) ( torderlim - 1 );
	for ( i = 1; i < torderlim - 1; i++ ){
		unscaledQuadratureWeights[i] = 2.0 * unscaledQuadratureWeights[i] / ( double ) ( torderlim - 1 );
	}
	unscaledQuadratureWeights[torderlim-1] = unscaledQuadratureWeights[torderlim-1] / ( double ) ( torderlim - 1 );


	struct tnode *testNode = NULL;
	testNode = malloc(sizeof(struct tnode));



	// Setup the test node
	double high=1.0;
	double low=-1.0;
	double mid = (high+low)/2.0;

	testNode->x_max=high;
	testNode->y_max=high;
	testNode->z_max=high;

	testNode->x_min=low;
	testNode->y_min=low;
	testNode->z_min=low;

	testNode->x_mid=mid;
	testNode->y_mid=mid;
	testNode->z_mid=mid;

	testNode->exist_ms=0;

	pc_comp_weights(testNode)


	for (i=0; i<torderlim; i++){
		printf("x, y, z weights: %2.5f, %2.5f, %2.5f \n", testNode->wx[i], testNode->wy[i], testNode->wz[i] );
	}

	return 0;
}
