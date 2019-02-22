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

#include "mkl.h"
#include <omp.h>


void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level)
{

    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild;
    
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
    
    //(*p)->x_min = xyzmm[0];
    //(*p)->x_max = xyzmm[1];
    //(*p)->y_min = xyzmm[2];
    //(*p)->y_max = xyzmm[3];
    //(*p)->z_min = xyzmm[4];
    //(*p)->z_max = xyzmm[5];
    

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
        
    } else {
        
        if (level < minlevel) minlevel = level;
        if (minpars > (*p)->numpar) minpars = (*p)->numpar;
        if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
        /* increment number of leaves */
        numleaves++;
    }

    return;

} /* END of function create_tree_n0 */



void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array)
{
    int i;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    tree_array->x_mid[p->node_index] = p->x_mid;
    tree_array->y_mid[p->node_index] = p->y_mid;
    tree_array->z_mid[p->node_index] = p->z_mid;

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



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order){

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

#pragma acc data copyin(tt[0:torderlim], \
		sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num]) \
		copy(clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num])

	addNodeToArray(troot, sources, clusters, order, numInterpPoints, pointsPerCluster);

	return;
}

void addNodeToArray(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;




//	double * testingQ;
//	make_vector(testingQ,numInterpPoints);
//	printf("number of interpolation points: %i\n\n", numInterpPoints);

	if (torderlim*torderlim*torderlim < p->numpar){ // don't compute moments for clusters that won't get used
//	if (torderlim*torderlim*torderlim < 1e10){ // don't compute moments for clusters that won't get used

//		make_vector(p->tx, torderlim);
//		make_vector(p->ty, torderlim);
//		make_vector(p->tz, torderlim);
//		pc_comp_ms(p, sources->x, sources->y, sources->z, sources->q, sources->w, clusters->q);
//		int k1,k2,k3;
//		int kk = -1;
//		for (k3 = 0; k3 < torderlim; k3++) {
//			for (k2 = 0; k2 < torderlim; k2++) {
//				for (k1 = 0; k1 < torderlim; k1++) {
//					kk++;
//					clusters->x[kk+startingIndex] = p->tx[k1];
//					clusters->y[kk+startingIndex] = p->ty[k2];
//					clusters->z[kk+startingIndex] = p->tz[k3];
//				}
//			}
//		}
//		if (p->node_index>72){
		pc_comp_ms_modifiedF(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
				clusters->x,clusters->y,clusters->z,clusters->q);
//		}

//		for (i=startingIndex;i<startingIndex+torderlim*torderlim*torderlim;i++) testingQ[i] = clusters->q[i];
//		for (i=startingIndex;i<startingIndex+torderlim*torderlim*torderlim; i++) clusters->q[i]=0.0;

//		pc_comp_ms_gpu(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
//				clusters->x,clusters->y,clusters->z,clusters->q);



//		for (i=startingIndex;i<startingIndex+torderlim*torderlim*torderlim;i++) testingQ[i] = clusters->q[i];

//		double maxDiff=0.0;
//		double relDiff=0.0;
//		int errindex;
//		for (i=startingIndex;i<startingIndex+torderlim*torderlim*torderlim;i++){
//			if (fabs(testingQ[i]-clusters->q[i])>maxDiff){
//				maxDiff=fabs(testingQ[i]-clusters->q[i]);
//				relDiff=fabs(testingQ[i]-clusters->q[i])/fabs(clusters->q[i]);
//				errindex = i-startingIndex;
//			}
//		}
////		printf("Starting index: %i\n", startingIndex);
//		printf("MaxDiff %e relDiff %e for node %i at index %i\n", maxDiff, relDiff, p->node_index, errindex);

		p->exist_ms = 1;


		// fill in arrays, starting at startingIndex


	}

	for (i = 0; i < p->num_children; i++) {
		addNodeToArray(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
	}

	return;
}





void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets, struct particles *clusters,
                 double *tpeng, double *EnP)
{
//	printf("Entered pc_treecoode.\n");
    /* local variables */
    int i, j;

    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;
    
#pragma acc data copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num]) \
		copy(EnP[0:targets->num])
    {
//#pragma acc data copyin(targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], EnP[0:targets->num]) \
//		copyout(EnP[0:targets->num]) \
//		present(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
//				clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num])
//    {

    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, EnP,
				clusters->x, clusters->y, clusters->z, clusters->q);
        }
    }
}

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



    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (torderlim*torderlim*torderlim < p->numpar) ) {


	int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
	int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
	int clusterStart = numberOfInterpolationPoints*p->node_index;

//	double clusterX[numberOfInterpolationPoints], clusterY[numberOfInterpolationPoints], clusterZ[numberOfInterpolationPoints], localMoments[numberOfInterpolationPoints];
//	double clusterXYZM[4*numberOfInterpolationPoints];


//	// Fill some local cluster arrays from the cluster itself
//	int k1,k2,k3;
//	kk = -1;
//	for (k3 = 0; k3 < torderlim; k3++) {
//		for (k2 = 0; k2 < torderlim; k2++) {
//			for (k1 = 0; k1 < torderlim; k1++) {
//				kk++;
//				localMoments[kk] = p->ms[kk];
//				clusterX[kk] = p->tx[k1];
//				clusterY[kk] = p->ty[k2];
//				clusterZ[kk] = p->tz[k3];
//			}
//		}
//	}





	double xi,yi,zi;
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


        if (p->num_children == 0) {
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






/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
//void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
//                    __restrict__ double *xS,__restrict__ double *yS, __restrict__ double *zS, __restrict__ double *qS, __restrict__ double *wS,
//					__restrict__ double *xT, __restrict__ double *yT, __restrict__ double *zT, __restrict__ double *qT, __restrict__ double *EnP)
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
# pragma acc kernels present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
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
        EnP[ii] += d_peng;
    }
    }
    return;

} /* END function pc_comp_direct */








/*
 * cp_comp_ms computes the moments for node p needed in the Taylor approximation
 */
//void pc_comp_ms(struct tnode *p, double __restrict__ *xS, double __restrict__ *yS, double __restrict__ *zS, double __restrict__ *qS, double __restrict__ *wS, double __restrict__ *clusterQ)
void pc_comp_ms(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS, double *clusterQ)
{

	int pointsPerCluster = torderlim*torderlim*torderlim;
	int startingIndex = p->node_index * pointsPerCluster;

//	printf("Entering pc_comp_ms.\n");
    int i, j, k1, k2, k3, kk;
    int a1exactind, a2exactind, a3exactind;
//    double dx, dy, dz, tx, ty, tz, qloc;
    double x0, x1, y0, y1, z0, z1;
    double sumA1, sumA2, sumA3;
    double xx, yy, zz, qq, ww;
//    double *xibeg, *yibeg, *zibeg, *qibeg, *wibeg;
    int xibeg, yibeg, zibeg, qibeg, wibeg;
    
    double w1i[torderlim], w2j[torderlim], w3k[torderlim], dj[torderlim];
    double Dd; //, **a1i, **a2j, **a3k;
    double *a1i, *a2j, *a3k;
    double *node_x, *node_y, *node_z;
    

//    xibeg = &(x[p->ibeg-1]);
//    yibeg = &(y[p->ibeg-1]);
//    zibeg = &(z[p->ibeg-1]);
//    qibeg = &(q[p->ibeg-1]);
//    wibeg = &(w[p->ibeg-1]);

    xibeg = p->ibeg-1;
	yibeg = p->ibeg-1;
	zibeg = p->ibeg-1;
	qibeg = p->ibeg-1;
	wibeg = p->ibeg-1;
    

    x0 = p->x_min;
    x1 = p->x_max;
    y0 = p->y_min;
    y1 = p->y_max;
    z0 = p->z_min;
    z1 = p->z_max;

//	x0 = p->x_min-1e-6*(p->x_max-p->x_min);
//	x1 = p->x_max+1e-6*(p->x_max-p->x_min);
//	y0 = p->y_min-1e-6*(p->y_max-p->y_min);
//	y1 = p->y_max+1e-6*(p->y_max-p->y_min);
//	z0 = p->z_min-1e-6*(p->z_max-p->z_min);
//	z1 = p->z_max+1e-6*(p->z_max-p->z_min);
    
    for (i = 0; i < torderlim; i++) {
    	p->tx[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        p->ty[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        p->tz[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

    }
    
//    make_matrix(a1i, torderlim, p->numpar);
//    make_matrix(a2j, torderlim, p->numpar);
//    make_matrix(a3k, torderlim, p->numpar);

	int pointsInNode = p->numpar;

//    make_vector(a1i, torderlim *pointsInNode);
//    make_vector(a2j, torderlim * pointsInNode);
//	make_vector(a3k, torderlim * pointsInNode);

	make_vector(a1i, torderlim );
	make_vector(a2j, torderlim );
	make_vector(a3k, torderlim );

//    make_vector(Dd, pointsInNode);
    make_vector(node_x, torderlim);
    make_vector(node_y, torderlim);
    make_vector(node_z, torderlim);
    
    dj[0] = 0.5;
    dj[torder] = 0.5;
    for (j = 1; j < torder; j++)
        dj[j] = 1.0;
    
    for (j = 0; j < torderlim; j++) {
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
        w2j[j] = w1i[j];
        w3k[j] = w1i[j];
        a1i[j]=0.0;
        a2j[j]=0.0;
        a3k[j]=0.0;
    }
    
    sumA1 = 0.0;
    sumA2 = 0.0;
    sumA3 = 0.0;
    

    for (i=0;i<torderlim;i++){
    	node_x[i] = p->tx[i];
    	node_y[i] = p->ty[i];
    	node_z[i] = p->tz[i];
    }

//
//    double * qMatrix;
//    make_vector(qMatrix, pointsInNode*torderlim*torderlim*torderlim );
//    for (i=0; i< pointsInNode*torderlim*torderlim*torderlim; i++){
//    	qMatrix[i]=0.0;
//    }

//#pragma acc region copy(clusterQ[startingIndex:startingIndex+torderlim*torderlim*torderlim])  \
//    create(qMatrix[0:pointsInNode*torderlim*torderlim*torderlim], a1i[0:torderlim], a2j[0:torderlim], a3k[0:torderlim]) present(xS,yS,zS,qS,wS)
//    {
//
//   	#pragma acc loop independent
//    for (i=0; i< pointsInNode*torderlim*torderlim*torderlim; i++){
//    	qMatrix[i]=0.0;
//    }

//   	#pragma acc loop independent
    for (i = 0; i < pointsInNode; i++) {
//        xx = xibeg[i];
//        yy = yibeg[i];
//        zz = zibeg[i];
//        qq = qibeg[i];
//        ww = wibeg[i];
    	xx = xS[xibeg+i];
		yy = yS[yibeg+i];
		zz = zS[zibeg+i];
		qq = qS[qibeg+i];
		ww = wS[wibeg+i];
        
        a1exactind = -1;
        a2exactind = -1;
        a3exactind = -1;
//		#pragma acc loop independent
        for (j = 0; j < torderlim; j++) {
//            a1i[i*torderlim + j] = w1i[j] / (xx - node_x[j]);
//            a2j[i*torderlim + j] = w2j[j] / (yy - node_y[j]);
//            a3k[i*torderlim + j] = w3k[j] / (zz - node_z[j]);

        	a1i[j] = w1i[j] / (xx - node_x[j]);
			a2j[j] = w2j[j] / (yy - node_y[j]);
			a3k[j] = w3k[j] / (zz - node_z[j]);
            


//            sumA1 += a1i[i*torderlim + j];
//			sumA2 += a2j[i*torderlim + j];
//			sumA3 += a3k[i*torderlim + j];

			sumA1 += a1i[j];
			sumA2 += a2j[j];
			sumA3 += a3k[j];
            
            if (fabs(xx - node_x[j]) < DBL_MIN) a1exactind = j;
            if (fabs(yy - node_y[j]) < DBL_MIN) a2exactind = j;
            if (fabs(zz - node_z[j]) < DBL_MIN) a3exactind = j;
        }
        
        if (a1exactind > -1) {
            sumA1 = 1.0;
            for (j = 0; j < torderlim; j++) a1i[j] = 0.0;
            a1i[a1exactind] = 1.0;
        }

        if (a2exactind > -1) {
            sumA2 = 1.0;
            for (j = 0; j < torderlim; j++) a2j[j] = 0.0;
            a2j[a2exactind] = 1.0;
        }

        if (a3exactind > -1) {
            sumA3 = 1.0;
            for (j = 0; j < torderlim; j++) a3k[j] = 0.0;
            a3k[a3exactind] = 1.0;
        }


        Dd = 1.0 / (sumA1 * sumA2 * sumA3);

		#pragma acc loop independent
        for (k3 = 0; k3 < torderlim; k3++) {
			#pragma acc loop independent
			for (k2 = 0; k2 < torderlim; k2++) {
				#pragma acc loop independent
				for (k1 = 0; k1 < torderlim; k1++) {
					j=k3*torderlim*torderlim + k2*torderlim + k1;
//					#pragma acc atomic update
					clusterQ[startingIndex + k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1] * a2j[k2] * a3k[k3] * Dd * qq * ww ;
//					qMatrix[i*torderlim*torderlim*torderlim + j] += a1i[k1] * a2j[k2] * a3k[k3] * Dd * qq * ww ;
//					p->ms[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1] * a2j[k2] * a3k[k3] * Dd * qibeg[i] * wibeg[i] ;
				}
			}
        }
        sumA1 = 0.0;
        sumA2 = 0.0;
        sumA3 = 0.0;
    }
//    } // end acc data region

//    printf("Filled qMatrix.\n");


//    double tempSum;
//	#pragma acc loop independent
//    for (j = 0;j<torderlim*torderlim*torderlim;j++){
//		tempSum=0.0;
//		#pragma acc loop independent
//	    for (i = 0; i < pointsInNode; i++) {
//	    	tempSum += qMatrix[i*torderlim*torderlim*torderlim + j];
//    	}
//	    clusterQ[startingIndex + j] = tempSum;
//    }
//    } // end acc data region

//    }
    
//    kk = -1;

////#pragma acc parallel loop
//    for (k3 = 0; k3 < torderlim; k3++) {
//        for (k2 = 0; k2 < torderlim; k2++) {
//            for (k1 = 0; k1 < torderlim; k1++) {
////                kk++;
//                for (i = 0; i < p->numpar; i++) {
////                	d_ms[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1][i] * a2j[k2][i] * a3k[k3][i] * Dd[i] * qibeg[i] * wibeg[i] ;  // in this case, multiple the function value by the quadrature weight.  Revisit for sing. subt.
//                	p->ms[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[i*torderlim + k1] * a2j[i*torderlim + k2] * a3k[i*torderlim + k3] * Dd[i] * qibeg[i] * wibeg[i] ;  // in this case, multiple the function value by the quadrature weight.  Revisit for sing. subt.
////                	d_ms2[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1][i] * a2j[k2][i] * a3k[k3][i] * Dd[i] *  wibeg[i] ;  // Anterpolate only the quadrature weights w_i for ms2, instead of w_i*f_i
//                    p->ms2[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[i*torderlim + k1] * a2j[i*torderlim + k2] * a3k[i*torderlim + k3] * Dd[i] *  wibeg[i] ;  // Anterpolate only the quadrature weights w_i for ms2, instead of w_i*f_i
//
//                    //if (p->ms[kk] != p->ms[kk])
//                    //    printf("%d, %d, %d, %d: %f, %f, %f, %f, %f, %f\n", k1, k2, k3, i,
//                    //    p->ms[kk], Dd[i], q[p->ibeg-1+i], a1i[k1][i], a2j[k2][i], a3k[k3][i]);
//                }
//            }
//        }
//    }
    

    free_vector(a1i);
    free_vector(a2j);
    free_vector(a3k);
//    free_vector(qMatrix);
//    free_vector(Dd);
    
    return;
    
} /* END function cp_comp_ms */


//void pc_comp_ms_gpu(struct tnode *p, double __restrict__ *xS, double __restrict__ *yS, double __restrict__ *zS, double __restrict__ *qS, double __restrict__ *wS,
//		double __restrict__ *clusterX, double __restrict__ *clusterY, double __restrict__ *clusterZ, double __restrict__ *clusterQ)
void pc_comp_ms_gpu(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
		double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{

	int pointsPerCluster = torderlim*torderlim*torderlim;
	int pointsInNode = p->numpar;
	int startingIndexInClusters = p->node_index * pointsPerCluster;
    int startingIndexInSources = p->ibeg-1;




//	printf("Entering pc_comp_ms.\n");
    int i, j, k;

    double x0, x1, y0, y1, z0, z1;  // bounding box
    double sumA1, sumA2, sumA3;
    double xx, yy, zz, qq, ww;

//    double w1i[torderlim], w2j[torderlim], w3k[torderlim], dj[torderlim];
    double *w1i, *w2j, *w3k, *dj;
    double *w3d;
    make_vector(w1i,torderlim);
    make_vector(w2j,torderlim);
    make_vector(w3k,torderlim);
    make_vector(dj,torderlim);
    make_vector(w3d,pointsPerCluster);
    double denominator, numerator, ak;



    // Set the bounding box.
    x0 = p->x_min-1e-6*(p->x_max-p->x_min);
    x1 = p->x_max+1e-6*(p->x_max-p->x_min);
    y0 = p->y_min-1e-6*(p->y_max-p->y_min);
    y1 = p->y_max+1e-6*(p->y_max-p->y_min);
    z0 = p->z_min-1e-6*(p->z_max-p->z_min);
    z1 = p->z_max+1e-6*(p->z_max-p->z_min);


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
    double nodeX[torderlim], nodeY[torderlim], nodeZ[torderlim];
//    double *nodeX, *nodeY, *nodeZ;
    for (i = 0; i < torderlim; i++) {
    	nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
    	nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
    	nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
    dj[0] = 0.5;
	dj[torder] = 0.5;
	for (j = 1; j < torder; j++){
		dj[j] = 1.0;
	}
    for (j = 0; j < torderlim; j++) {
            w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
            w2j[j] = w1i[j];
            w3k[j] = w1i[j];
	}

    // Fill in the global arrays of cluster point data.  This cluster only fills a portion.   Also fill local array of weights (outer product of 1d weights)
    int kk = -1;
    for (k=0;k<torderlim;k++){
    	for (j=0;j<torderlim;j++){
    		for (i=0;i<torderlim;i++){
    			kk++;
				clusterX[kk+startingIndexInClusters] = nodeX[i];
				clusterY[kk+startingIndexInClusters] = nodeY[j];
				clusterZ[kk+startingIndexInClusters] = nodeZ[k];
				w3d[kk] = w3k[k]*w2j[j]*w1i[i];
    		}
    	}
    }


    double cx, cy, cz, cw, px, py, pz, pw, pq;  // coordinates of cluster interpolation point and particle
#pragma acc kernels present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ)
	#pragma acc loop independent
    for (i = 0; i < pointsPerCluster; i++) {

		#pragma acc loop independent
        for (j = 0; j < pointsInNode; j++) {


        	px = xS[startingIndexInSources + j];
        	py = yS[startingIndexInSources + j];
        	pz = zS[startingIndexInSources + j];
        	pq = qS[startingIndexInSources + j];
        	pw = wS[startingIndexInSources + j];



			// Compute the denominator
        	sumA1 = 0.0;
			sumA2 = 0.0;
			sumA3 = 0.0;
			#pragma acc loop independent
			for (k=0; k<torderlim;k++){
				cx = nodeX[k];
				cy = nodeY[k];
				cz = nodeZ[k];
				sumA1 += w1i[k] / (px - cx);
				sumA2 += w2j[k] / (py - cy);
				sumA3 += w3k[k] / (pz - cz);
			}

			denominator = (sumA1 * sumA2 * sumA3);


			// Compute the numerator
			cx = clusterX[startingIndexInClusters+i];
			cy = clusterY[startingIndexInClusters+i];
			cz = clusterZ[startingIndexInClusters+i];
			cw = w3d[i]; // product of the three 1-dimensional weights


			numerator = cw / (px - cx) / (py - cy) / (pz - cz);

			// Construct weight
			ak = numerator / denominator;


			// Increment modified weight by adding contribution from j^th point
			clusterQ[startingIndexInClusters + i] += ak * pq * pw;
        }

    }

    free_vector(w1i);
    free_vector(w2j);
    free_vector(w3k);
    free_vector(dj);
    free_vector(w3d);


    return;

} /* END function cp_comp_ms_gpu */



void pc_comp_ms_modifiedF(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
		double *clusterX, double *clusterY, double *clusterZ, double *clusterQ){

	int i,j,k;
	int pointsPerCluster = torderlim*torderlim*torderlim;
	int pointsInNode = p->numpar;
	int startingIndexInClusters = p->node_index * pointsPerCluster;
	int startingIndexInSources = p->ibeg-1;

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







	// Set the bounding box.


//	x0 = p->x_min-1e-13*(p->x_max-p->x_min);  // 1e-15 fails for large meshes, mysteriously.
//	x1 = p->x_max+2e-13*(p->x_max-p->x_min);
//	y0 = p->y_min-1e-13*(p->y_max-p->y_min);
//	y1 = p->y_max+2e-13*(p->y_max-p->y_min);
//	z0 = p->z_min-1e-13*(p->z_max-p->z_min);
//	z1 = p->z_max+2e-13*(p->z_max-p->z_min);

	x0 = p->x_min;  // 1e-15 fails for large meshes, mysteriously.
	x1 = p->x_max;
	y0 = p->y_min;
	y1 = p->y_max;
	z0 = p->z_min;
	z1 = p->z_max;

	// Make and zero-out arrays to store denominator sums
	double sumX, sumY, sumZ;


#pragma acc kernels present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ,tt) \
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
	dj[0] = 0.5;
	dj[torder] = 0.5;
	#pragma acc loop independent
	for (j = 1; j < torder; j++){
		dj[j] = 1.0;
	}
	#pragma acc loop independent
	for (j = 0; j < torderlim; j++) {
		weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
	}


	// Compute modified f values
	double sx,sy,sz,cx,cy,cz,denominator,w;



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
                              int **tree_inter_list, int **direct_inter_list)
{
    /* local variables */
    int i, j;
    int tree_index_counter;
    int direct_index_counter;

    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < numnodes; j++) {
            tree_inter_list[i][j] = -1;
        }
        for (j = 0; j < numleaves; j++) {
            direct_inter_list[i][j] = -1;
        }
    }
    
    for (i = 0; i < batches->num; i++) {
        tree_index_counter = 0;
        direct_index_counter = 0;
        
        pc_compute_interaction_list(p,
                batches->index[i], batches->center[i], batches->radius[i],
                tree_inter_list[i], direct_inter_list[i],
                &tree_index_counter, &direct_index_counter);
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

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
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
