/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"

#include <omp.h>


void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level)
{
//	printf("Entering pc_create_tree_n0.\n");

    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild, idx;
    
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];
    
    
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }
                        

    (*p) = malloc(sizeof(struct tnode));


    /* increment number of nodes */
    (*p)->node_index = numnodes;
    numnodes++;

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;
    
    (*p)->x_min = minval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(sources->z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(sources->z + ibeg - 1, (*p)->numpar);
    


    /*compute aspect ratio*/
    xl = (*p)->x_max - (*p)->x_min;
    yl = (*p)->y_max - (*p)->y_min;
    zl = (*p)->z_max - (*p)->z_min;
        
    lmax = max3(xl, yl, zl);
    t1 = lmax;
    t2 = min3(xl, yl, zl);


    if (t2 != 0.0)
        (*p)->aspect = t1/t2;
    else
        (*p)->aspect = 0.0;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    t1 = (*p)->x_max - (*p)->x_mid;
    t2 = (*p)->y_max - (*p)->y_mid;
    t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    /*set particle limits, tree level of node, and nullify child pointers*/
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    if (maxlevel < level) maxlevel = level;

    (*p)->num_children = 0;
    for (i = 0; i < 8; i++)
        (*p)->child[i] = NULL;

    
    if ((*p)->numpar > maxparnode) {

    /*
     * set IND array to 0, and then call PARTITION_8 routine.
     * IND array holds indices of the eight new subregions.
     * Also, setup XYZMMS array in the case that SHRINK = 1.
     */
        xyzmms[0][0] = (*p)->x_min;
        xyzmms[1][0] = (*p)->x_max;
        xyzmms[2][0] = (*p)->y_min;
        xyzmms[3][0] = (*p)->y_max;
        xyzmms[4][0] = (*p)->z_min;
        xyzmms[5][0] = (*p)->z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        x_mid = (*p)->x_mid;
        y_mid = (*p)->y_mid;
        z_mid = (*p)->z_mid;

        pc_partition_8(sources->x, sources->y, sources->z, sources->q, sources->w,
                       xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;
        if (level==0){  // level 0, spawn openMP threads

//#pragma omp parallel
//        	{
//		#pragma omp for private (j, idx, loclev,lxyzmm)
        for (i = 0; i < numposchild; i++) {
//			#pragma omp barrier
            if (ind[i][0] <= ind[i][1]) {

//				#pragma omp critical
//            	{
                (*p)->num_children = (*p)->num_children + 1;
                idx = (*p)->num_children - 1;
                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
//            	}
                pc_create_tree_n0(&((*p)->child[idx]),
                                  sources, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev);

            }
        }
//        	} // end omp parallel
        }else{ // not level 0, don't spawn openMP threads
        	for (i = 0; i < numposchild; i++) {
        	            if (ind[i][0] <= ind[i][1]) {

        	                (*p)->num_children = (*p)->num_children + 1;

        	                for (j = 0; j < 6; j++)
        	                    lxyzmm[j] = xyzmms[j][i];

        	                pc_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
        	                                  sources, ind[i][0], ind[i][1],
        	                                  maxparnode, lxyzmm, loclev);
        	            }
        	        }
        }
        
    } else {
        
        if (level < minlevel) minlevel = level;
        if (minpars > (*p)->numpar) minpars = (*p)->numpar;
        if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
        /* increment number of leaves */
        numleaves++;
    }

//    printf("Exiting pc_create_tree_n0.\n");
    return;

} /* END of function create_tree_n0 */



void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array)
{
//	printf("Entering pc_create_tree_array.\n");
    int i;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    tree_array->x_mid[p->node_index] = p->x_mid;
    tree_array->y_mid[p->node_index] = p->y_mid;
    tree_array->z_mid[p->node_index] = p->z_mid;

    tree_array->x_min[p->node_index] = p->x_min;
	tree_array->y_min[p->node_index] = p->y_min;
	tree_array->z_min[p->node_index] = p->z_min;

	tree_array->x_max[p->node_index] = p->x_max;
	tree_array->y_max[p->node_index] = p->y_max;
	tree_array->z_max[p->node_index] = p->z_max;

    tree_array->ibeg[p->node_index] = p->ibeg;
    tree_array->iend[p->node_index] = p->iend;

    for (i = 0; i < p->num_children; i++) {
        pc_create_tree_array(p->child[i], tree_array);
    }

    return;

} /* END of function create_tree_n0 */



void pc_partition_8(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    /* local variables */
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {

        pc_partition(x, y, z, q, w, orderarr, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];

        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;

    }

    if (yl >= critlen) {

        for (i = 0; i < *numposchild; i++) {
            pc_partition(y, x, z, q, w, orderarr, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);
                        
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    if (zl >= critlen) {

        for (i = 0; i < *numposchild; i++) {
            pc_partition(z, x, y, q, w, orderarr, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);
                        
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function partition_8 */



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order, int numDevices, int numThreads, struct tnode_array * tree_array){

	int pointsPerCluster = (order+1)*(order+1)*(order+1);
	int numInterpPoints = numnodes * pointsPerCluster;
	make_vector(clusters->x, numInterpPoints);
	make_vector(clusters->y, numInterpPoints);
	make_vector(clusters->z, numInterpPoints);
	make_vector(clusters->q, numInterpPoints);
	make_vector(clusters->w, numInterpPoints);  // will be used in singularity subtraction
//	memset(clusters->q,0,numInterpPoints*sizeof(double));
	clusters->num=numInterpPoints;

	for (int i=0; i< numInterpPoints; i++){
		clusters->q[i]=0.0;
		clusters->w[i]=0.0;
	}


#pragma omp parallel num_threads(numThreads)
//#pragma omp parallel num_threads(1)
	{
		if (omp_get_thread_num()<numDevices){
			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
		}

		int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
		if (this_thread==0){printf("numDevices: %i\n", numDevices);}
		if (this_thread==0){printf("num_threads: %i\n", num_threads);}
		printf("this_thread: %i\n", this_thread);

		double *tempQ, *tempX, *tempY, *tempZ;
		make_vector(tempX,clusters->num);
		make_vector(tempY,clusters->num);
		make_vector(tempZ,clusters->num);
		make_vector(tempQ,clusters->num);
		for (int i = 0; i < clusters->num; i++)
		{
			tempX[i] = 0.0;
			tempY[i] = 0.0;
			tempZ[i] = 0.0;
			tempQ[i] = 0.0;
		}

		double *xS = sources->x;
		double *yS = sources->y;
		double *zS = sources->z;
		double *qS = sources->q;
		double *wS = sources->w;

		double *xC = clusters->x;
		double *yC = clusters->y;
		double *zC = clusters->z;
		double *qC = clusters->q;

		int clusterNum = clusters->num;
		int sourceNum = sources->num;

#pragma acc data copyin(tt[0:torderlim], \
		xS[0:sourceNum], yS[0:sourceNum], zS[0:sourceNum], qS[0:sourceNum], wS[0:sourceNum]) \
		copy(tempX[0:clusterNum], tempY[0:clusterNum], tempZ[0:clusterNum], tempQ[0:clusterNum])
		{
//	addNodeToArray(troot, sources, clusters, order, numInterpPoints, pointsPerCluster);

			#pragma omp for schedule(guided)
			for (int i=1;i<numnodes; i++){  // start from i=1, as we do not need to compute moments for root, and this is realtively expensive.
				pc_comp_ms_modifiedF(tree_array, i, xS, yS, zS, qS, wS, \
							tempX,tempY,tempZ,tempQ);
				}
			#pragma acc wait
			} // end ACC DATA REGION

		int counter=0;
		for (int j = 0; j < clusters->num; j++)
		{

			if (tempQ[j]!=0.0){
				clusters->x[j] = tempX[j];
				clusters->y[j] = tempY[j];
				clusters->z[j] = tempZ[j];
				clusters->q[j] += tempQ[j];
			}



		} // end j loop
//		#pragma omp barrier



		free_vector(tempX);
		free_vector(tempY);
		free_vector(tempZ);
		free_vector(tempQ);

		} // end OMP PARALLEL REGION

//	printf("outside omp parallel region: %f, %f\n\n", clusters->q[0], clusters->q[213599]);
//	double tempSum = sum(clusters->q, clusters->num);
//	printf("\n\n\nSum of cluster q: %f\n\n\n", tempSum);

//#pragma acc data copyin(clusters->q[0:clusters->num])



/////// REDO WITH ONE THREAD THEN COMPARE
//#pragma omp parallel num_threads(1)
////#pragma omp parallel num_threads(1)
//	{
//		if (omp_get_thread_num()<numDevices){
//			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
//		}
//
//		int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
//		if (this_thread==0){printf("numDevices: %i\n", numDevices);}
//		if (this_thread==0){printf("num_threads: %i\n", num_threads);}
//		printf("this_thread: %i\n", this_thread);
//
//		double *tempQ2;
//		make_vector(tempQ2,clusters->num);
//		for (int i = 0; i < clusters->num; i++)
//			tempQ2[i] = 0.0;
//
//		double *xS = sources->x;
//		double *yS = sources->y;
//		double *zS = sources->z;
//		double *qS = sources->q;
//		double *wS = sources->w;
//
//		double *xC = clusters->x;
//		double *yC = clusters->y;
//		double *zC = clusters->z;
//		double *qC = clusters->q;
//
//		int clusterNum = clusters->num;
//		int sourceNum = sources->num;
//
//#pragma acc data copyin(tt[0:torderlim], \
//		xS[0:sourceNum], yS[0:sourceNum], zS[0:sourceNum], qS[0:sourceNum], wS[0:sourceNum]) \
//		copy(xC[0:clusterNum], yC[0:clusterNum], zC[0:clusterNum], tempQ2[0:clusterNum])
//			{
////	addNodeToArray(troot, sources, clusters, order, numInterpPoints, pointsPerCluster);
//
//			#pragma omp for schedule(guided)
//			for (int i=1;i<numnodes; i++){  // start from i=1, as we do not need to compute moments for root, and this is realtively expensive.
//				pc_comp_ms_modifiedF(tree_array, i, xS, yS, zS, qS, wS, \
//							xC,yC,zC,tempQ2);
//				}
//			#pragma acc wait
//			} // end ACC DATA REGION
////			#pragma acc wait
////		#pragma omp critical
////			{
//		int counter=0;
//		printf("clusters->num = %i\n",clusters->num);
//
//		printf("Comparing tempQ2 and clusters->q elementwise.\n");
//		for (int j = 0; j < clusters->num; j++)
//		{
//			if (abs(tempQ2[j] - clusters->q[j])>1e-15){
//				printf("WARNING: tempQ2 and clusters->q not agreeing at index %i\n",j);
//			}
//
//
//
//		} // end j loop
//		printf("done/n");
////			} // end omp critical
//		#pragma omp barrier
////			for (int j = 0; j < clusters->num; j++)
////					{
////						if (clusters->q[j]==0.0)
////						{
////							printf("clusters->q[%i] = 0.0\n", j);
////						}
////					}
//
//
//		free_vector(tempQ2);
//
//		} // end OMP PARALLEL REGION

	return;
}

void addNodeToArray(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;


	if (1==1){ // don't compute moments for clusters that won't get used

//		pc_comp_ms_modifiedF(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
//				clusters->x,clusters->y,clusters->z,clusters->q);

		printf("Commented out old pc_comp_ms_modifiedF.\n");

		p->exist_ms = 1;


	}

	for (i = 0; i < p->num_children; i++) {
		addNodeToArray(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
	}

	return;
}

void addNodeToArray_nonRecursive(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;



//	pc_comp_ms_modifiedF(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
//			clusters->x,clusters->y,clusters->z,clusters->q);

	p->exist_ms = 1;



//	for (i = 0; i < p->num_children; i++) {
//		addNodeToArray(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
//	}

	return;
}





void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets, struct particles *clusters,
                 double *tpeng, double *EnP, int numDevices, int numThreads)
{
//	printf("Entered pc_treecoode.\n");
    /* local variables */
    int i, j;

    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;
    


#pragma omp parallel num_threads(numThreads)
	{
    	if (omp_get_thread_num()<numDevices){
    		acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
    	}

        int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
		if (this_thread==0){printf("numDevices: %i\n", numDevices);}
		if (this_thread==0){printf("num_threads: %i\n", num_threads);}
		printf("this_thread: %i\n", this_thread);

		double *EnP2;
		make_vector(EnP2,targets->num);
		for (i = 0; i < targets->num; i++)
			EnP2[i] = 0.0;




#pragma acc data copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num]) \
		copy(EnP2[0:targets->num])
    {


	#pragma omp for private(j)
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, EnP2,
				clusters->x, clusters->y, clusters->z, clusters->q);
        }
    }



    } // end acc data region
    for (int k = 0; k < targets->num; k++){
    	if (EnP2[k] != 0.0)
			EnP[k] += EnP2[k];
		}

    free_vector(EnP2);
	} // end omp parallel region

    printf("Exited the main comp_pc call.\n");
    *tpeng = sum(EnP, targets->num);

    return;

} /* END of function pc_treecode */



void compute_pc(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM )
{
//	printf("Entering compute_cp.  Batch start: %d.  Batch end: %d.\n", batch_ind[0]-1, batch_ind[1]);
    /* local variables */
    double dist;
    double tx, ty, tz;
    int i, j, k, ii, kk;
    double dxt,dyt,dzt,tempPotential;
    double temp_i[torderlim], temp_j[torderlim], temp_k[torderlim];

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);

    int smallEnoughLeaf=0;
    if (p->numpar < torderlim*torderlim*torderlim){
    	smallEnoughLeaf=1;
    }else{
    	smallEnoughLeaf=0;
    }

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (smallEnoughLeaf==0) ) {
//	if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) ) {


	int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
	int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
	int clusterStart = numberOfInterpolationPoints*p->node_index;



	double xi,yi,zi ;
	int batchStart = batch_ind[0] - 1;
//	double pointvals[4];
# pragma acc kernels present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
    {
	#pragma acc loop independent
	for (i = 0; i < numberOfTargets; i++){
		tempPotential = 0.0;
		xi = xT[ batchStart + i];
		yi = yT[ batchStart + i];
		zi = zT[ batchStart + i];
		for (j = 0; j < numberOfInterpolationPoints; j++){

			// Compute x, y, and z distances between target i and interpolation point j
			dxt = xi - clusterX[clusterStart + j];
			dyt = yi - clusterY[clusterStart + j];
			dzt = zi - clusterZ[clusterStart + j];


//			tempPotential += localMoments[j] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
			tempPotential += clusterM[clusterStart + j] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

						}

		EnP[batchStart + i] += tempPotential;
	}
    }



    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if ( (p->num_children == 0) | (smallEnoughLeaf==1) ) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
                compute_pc(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP,
							clusterX, clusterY, clusterZ, clusterM);
            }
        }
    }

    return;

} /* END of function compute_pc */



void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
                     double *xS, double *yS,  double *zS,  double *qS,  double *wS,
					 double *xT,  double *yT,  double *zT,  double *qT,  double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;

    int batch_start=batch_ibeg - 1;
    int batch_end = batch_iend;

    int source_start=ibeg - 1;
    int source_end=iend;

    double d_peng, r;
    int streamID = rand() % 2;
# pragma acc kernels async(streamID) present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
//# pragma acc kernels present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
    {
	#pragma acc loop independent
    for (ii = batch_start; ii < batch_end; ii++) {
        d_peng = 0.0;
        for (i = source_start; i < source_end; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            r = sqrt(tx*tx + ty*ty + tz*tz);
            if (r > DBL_MIN) {
            	d_peng += qS[i] * wS[i] / r;
            }
        }
		#pragma acc atomic
        EnP[ii] += d_peng;
    }
    }
    return;

} /* END function pc_comp_direct */



void pc_comp_ms_modifiedF(struct tnode_array * tree_array, int idx, double *xS, double *yS, double *zS, double *qS, double *wS,
		double *clusterX, double *clusterY, double *clusterZ, double *clusterQ){

	int i,j,k;
	int pointsPerCluster = torderlim*torderlim*torderlim;
	int pointsInNode = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
	int startingIndexInClusters = idx * pointsPerCluster;
	int startingIndexInSources = tree_array->ibeg[idx]-1;

	double x0, x1, y0, y1, z0, z1;  // bounding box


	double weights[torderlim];
	double dj[torderlim];
	double *modifiedF;
	make_vector(modifiedF,pointsInNode);

	double nodeX[torderlim], nodeY[torderlim], nodeZ[torderlim];

	int *exactIndX, *exactIndY, *exactIndZ;
	make_vector(exactIndX, pointsInNode);
	make_vector(exactIndY, pointsInNode);
	make_vector(exactIndZ, pointsInNode);



	x0 = tree_array->x_min[idx];  // 1e-15 fails for large meshes, mysteriously.
	x1 = tree_array->x_max[idx];
	y0 = tree_array->y_min[idx];
	y1 = tree_array->y_max[idx];
	z0 = tree_array->z_min[idx];
	z1 = tree_array->z_max[idx];

	// Make and zero-out arrays to store denominator sums
	double sumX, sumY, sumZ;
	double sx,sy,sz,cx,cy,cz,denominator,w;

	int streamID = rand() % 4;
//	async(streamID)
#pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ,tt) \
	create(modifiedF[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
			nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],weights[0:torderlim],dj[0:torderlim])
	{


	#pragma acc loop independent
	for (j=0;j<pointsInNode;j++){
		modifiedF[j] = qS[startingIndexInSources+j]*wS[startingIndexInSources+j];
		exactIndX[j] = -1;
		exactIndY[j] = -1;
		exactIndZ[j] = -1;
	}

	//  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
	#pragma acc loop independent
	for (i = 0; i < torderlim; i++) {
		nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
		nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
		nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

	}

	// Compute weights
//	dj[0] = 0.5;
//	dj[torder] = 0.5;
	#pragma acc loop independent
	for (j = 0; j < torder+1; j++){
		dj[j] = 1.0;
		if (j==0) dj[j] = 0.5;
		if (j==torder) dj[j]=0.5;
	}

	#pragma acc loop independent
	for (j = 0; j < torderlim; j++) {
		weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
	}


	// Compute modified f values



	#pragma acc loop independent
	for (i=0; i<pointsInNode;i++){ // loop through the source points

		sumX=0.0;
		sumY=0.0;
		sumZ=0.0;

		sx = xS[startingIndexInSources+i];
		sy = yS[startingIndexInSources+i];
		sz = zS[startingIndexInSources+i];


		#pragma acc loop independent
		for (j=0;j<torderlim;j++){  // loop through the degree

			cx = sx-nodeX[j];
			cy = sy-nodeY[j];
			cz = sz-nodeZ[j];

			if (fabs(cx)<DBL_MIN) exactIndX[i]=j;
			if (fabs(cy)<DBL_MIN) exactIndY[i]=j;
			if (fabs(cz)<DBL_MIN) exactIndZ[i]=j;

			// Increment the sums
			w = weights[j];
			sumX += w / (cx);
			sumY += w / (cy);
			sumZ += w / (cz);

		}

//		denominator = sumX*sumY*sumZ;

		denominator = 1.0;
		if (exactIndX[i]==-1) denominator *= sumX;
		if (exactIndY[i]==-1) denominator *= sumY;
		if (exactIndZ[i]==-1) denominator *= sumZ;

		modifiedF[i] /= denominator;

	}


	// Compute moments for each interpolation point
	double numerator, xn, yn, zn, temp;
	int k1, k2, k3, kk;
	double w1,w2,w3;

	#pragma acc loop independent
	for (j=0;j<pointsPerCluster;j++){ // loop over interpolation points, set (cx,cy,cz) for this point
		// compute k1, k2, k3 from j
		k1 = j%torderlim;
		kk = (j-k1)/torderlim;
		k2 = kk%torderlim;
		kk = kk - k2;
		k3 = kk / torderlim;

		cz = nodeZ[k3];
		w3 = weights[k3];

		cy = nodeY[k2];
		w2 = weights[k2];

		cx = nodeX[k1];
		w1 = weights[k1];


		// Fill cluster X, Y, and Z arrays
		clusterX[startingIndexInClusters + j] = cx;
		clusterY[startingIndexInClusters + j] = cy;
		clusterZ[startingIndexInClusters + j] = cz;


		// Increment cluster Q array
		temp = 0.0;
		#pragma acc loop independent
		for (i=0;i<pointsInNode; i++){  // loop over source points
			sx = xS[startingIndexInSources+i];
			sy = yS[startingIndexInSources+i];
			sz = zS[startingIndexInSources+i];

			numerator=1.0;


			// If exactInd[i] == -1, then no issues.
			// If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
			if (exactIndX[i]==-1){
				numerator *=  w1 / (sx - cx);
			}else{
//				printf("ExactIndX != -1\n");
				if (exactIndX[i]!=k1) numerator *= 0;
			}

			if (exactIndY[i]==-1){
				numerator *=  w2 / (sy - cy);
			}else{
				if (exactIndY[i]!=k2) numerator *= 0;
			}

			if (exactIndZ[i]==-1){
				numerator *=  w3 / (sz - cz);
			}else{
				if (exactIndZ[i]!=k3) numerator *= 0;
			}





			temp += numerator*modifiedF[i];


		}
//		#pragma acc atomic
		clusterQ[startingIndexInClusters + j] += temp;

	}


	}

	free_vector(modifiedF);
	free_vector(exactIndX);
	free_vector(exactIndY);
	free_vector(exactIndZ);


	return;
}




void pc_make_interaction_list(struct tnode *p, struct batch *batches,
                              int *tree_inter_list, int *direct_inter_list)
{
    /* local variables */
    int i;
    int tree_index_counter;
    int direct_index_counter;

    for (i = 0; i < batches->num * numnodes; i++) {
        tree_inter_list[i] = -1;
    }

    for (i = 0; i < batches->num * numleaves; i++) {
        direct_inter_list[i] = -1;
    }
    
    for (i = 0; i < batches->num; i++) {
        tree_index_counter = 0;
        direct_index_counter = 0;
        
        pc_compute_interaction_list(p,
                batches->index[i], batches->center[i], batches->radius[i],
                &(tree_inter_list[i*numnodes]), &(direct_inter_list[i*numleaves]),
                &tree_index_counter, &direct_index_counter);

        batches->index[i][2] =tree_index_counter;
        batches->index[i][3] =direct_index_counter;

    }

    return;

} /* END of function pc_treecode */




void pc_compute_interaction_list(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list,
                int *tree_index_counter, int *direct_index_counter)
{
    /* local variables */
    double tx, ty, tz, dist;
    int i;

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)
       && (torder*torder*torder < p->numpar)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */

        batch_tree_list[*tree_index_counter] = p->node_index;
        (*tree_index_counter)++;

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            batch_direct_list[*direct_index_counter] = p->node_index;
            (*direct_index_counter)++;

        } else {
            for (i = 0; i < p->num_children; i++) {
                pc_compute_interaction_list(p->child[i], batch_ind, batch_mid, batch_rad,
                           batch_tree_list, batch_direct_list,
                           tree_index_counter, direct_index_counter);
            }
        }
    }

    return;

} /* END of function compute_pc */

void pc_interaction_list_treecode(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double *EnP, int numDevices, int numThreads)
{
	    int i, j;

	    for (i = 0; i < targets->num; i++)
	        EnP[i] = 0.0;


	#pragma omp parallel num_threads(numThreads)
		{
	    	if (omp_get_thread_num()<numDevices){
	    		acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
	    	}

	        int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
			if (this_thread==0){printf("numDevices: %i\n", numDevices);}
			if (this_thread==0){printf("num_threads: %i\n", num_threads);}
			printf("this_thread: %i\n", this_thread);

			double *EnP2, *EnP3;
			make_vector(EnP2,targets->num);
			make_vector(EnP3,targets->num);
			for (i = 0; i < targets->num; i++)
				EnP3[i] = 0.0;
				EnP2[i] = 0.0;


			double *xS = sources->x;
			double *yS = sources->y;
			double *zS = sources->z;
			double *qS = sources->q;
			double *wS = sources->w;

			double *xT = targets->x;
			double *yT = targets->y;
			double *zT = targets->z;
			double *qT = targets->q;

			double *xC = clusters->x;
			double *yC = clusters->y;
			double *zC = clusters->z;
			double *qC = clusters->q;

//			printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//			printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

			int * ibegs = tree_array->ibeg;
			int * iends = tree_array->iend;

	#pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], qS[0:sources->num], wS[0:sources->num], \
			xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
			xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num], tree_inter_list[0:numnodes*batches->num], \
			direct_inter_list[0:batches->num * numleaves], ibegs[0:numnodes], iends[0:numnodes]) \
			copy(EnP3[0:targets->num], EnP2[0:targets->num])
	    {


		int batch_ibeg, batch_iend, node_index;
		double dist;
		double tx, ty, tz;
		int i, j, k, ii, jj;
		double dxt,dyt,dzt,tempPotential;
		double temp_i[torderlim], temp_j[torderlim], temp_k[torderlim];

		int source_start;
		int source_end;

		double d_peng, r;
		double xi,yi,zi;

		int numberOfTargets;
		int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
		int clusterStart, batchStart;

		int numberOfClusterApproximations, numberOfDirectSums;

		int streamID;
		#pragma omp for private(j,ii,jj,batch_ibeg,batch_iend,numberOfClusterApproximations,numberOfDirectSums,numberOfTargets,batchStart,node_index,clusterStart,streamID)
	    for (i = 0; i < batches->num; i++) {
	    	batch_ibeg = batches->index[i][0];
			batch_iend = batches->index[i][1];
			numberOfClusterApproximations = batches->index[i][2];
			numberOfDirectSums = batches->index[i][3];


			numberOfTargets = batch_iend - batch_ibeg + 1;
			batchStart =  batch_ibeg - 1;

			for (j = 0; j < numberOfClusterApproximations; j++) {
				node_index = tree_inter_list[i * numnodes + j];
				clusterStart = numberOfInterpolationPoints*node_index;



				streamID = j%3;
				#pragma acc kernels async(streamID) //present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
				{
				#pragma acc loop independent
				for (ii = 0; ii < numberOfTargets; ii++){
					tempPotential = 0.0;
					xi = xT[ batchStart + ii];
					yi = yT[ batchStart + ii];
					zi = zT[ batchStart + ii];

					for (jj = 0; jj < numberOfInterpolationPoints; jj++){

						// Compute x, y, and z distances between target i and interpolation point j
						dxt = xi - xC[clusterStart + jj];
						dyt = yi - yC[clusterStart + jj];
						dzt = zi - zC[clusterStart + jj];
						tempPotential += qC[clusterStart + jj] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

									}
					#pragma acc atomic
					EnP3[batchStart + ii] += tempPotential;
				}
					} // end kernel
	        } // end for loop over cluster approximations

			for (j = 0; j < numberOfDirectSums; j++) {

				node_index = direct_inter_list[i * numleaves + j];

				source_start=ibegs[node_index]-1;
				source_end=iends[node_index];


				streamID = j%3;
		    # pragma acc kernels async(streamID)
		    {
			#pragma acc loop independent
		    for (ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
		        d_peng = 0.0;
		        for (jj = source_start; jj < source_end; jj++) {
		            tx = xS[jj] - xT[ii];
		            ty = yS[jj] - yT[ii];
		            tz = zS[jj] - zT[ii];
		            r = sqrt(tx*tx + ty*ty + tz*tz);
		            if (r > DBL_MIN) {
		            	d_peng += qS[jj] * wS[jj] / r;
		            }
		        }
				#pragma acc atomic
		        EnP2[ii] += d_peng;
		    }
		    } // end kernel
			} // end loop over number of leaves
//#pragma acc wait

	    } // end loop over target batches


		#pragma acc wait
		#pragma omp barrier
	    } // end acc data region


	    for (int k = 0; k < targets->num; k++){
	    	if (EnP2[k] != 0.0)
				#pragma omp critical
	    	{
				EnP[k] += EnP2[k];
				EnP[k] += EnP3[k];
	    	} // end omp critical
			}

	    free_vector(EnP2);
	    free_vector(EnP3);
		} // end omp parallel region

	    printf("Exited the main comp_pc call.\n");
	    *tpeng = sum(EnP, targets->num);

	    return;

	} /* END of function pc_treecode */
