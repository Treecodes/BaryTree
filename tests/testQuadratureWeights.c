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



int main()
{
	struct tnode *testNode = NULL;
	testNode = malloc(sizeof(struct tnode));

	// Set the global variables
	torder = 4;
	torderlim = torder+1;

	// Set bounds for the cluster
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


	for (int i=0; i<torderlim; i++){
		printf("x, y, z weights: %2.5f, %2.5f, %2.5f \n", testNode->wx[i], testNode->wy[i], testNode->wz[i] );
	}

	return 0;
}
