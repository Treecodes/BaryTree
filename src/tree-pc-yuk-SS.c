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
#include "tools.h"

#include "partition.h"
#include "tree.h"
#include "mkl.h"



void pc_treecode_yuk_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappa, double *tpeng, double *EnP)
{
    /* local variables */
    int i, j;
    
    for (i = 0; i < targets->num; i++){
        EnP[i] = 4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
    }


#pragma acc data copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num]) \
		copy(EnP[0:targets->num])
    {
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc_yuk_SS(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, kappa, EnP,
				clusters->x, clusters->y, clusters->z, clusters->q, clusters->w);
        }
    }
    }
    
//    for (i=0;i<targets->num;i++){
//    	EnP[i] *= sources->w[i];
//    }
    *tpeng = sum(EnP, targets->num);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_yuk_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM, double * clusterM2 )
{
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

    int smallEnoughLeaf;
	if (torderlim*torderlim*torderlim < p->numpar){
		smallEnoughLeaf=0;
	}else{
		smallEnoughLeaf=1;
	}
	if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (smallEnoughLeaf==0)  ) {


		int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
		int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
		int clusterStart = numberOfInterpolationPoints*p->node_index;



		double xi,yi,zi,qi,r;
		int batchStart = batch_ind[0] - 1;
# pragma acc kernels present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM, clusterM2)
		{
		#pragma acc loop independent
		for (i = 0; i < numberOfTargets; i++){
			tempPotential = 0.0;
			xi = xT[ batchStart + i];
			yi = yT[ batchStart + i];
			zi = zT[ batchStart + i];
			qi = qT[ batchStart + i];
			#pragma acc loop independent
			for (j = 0; j < numberOfInterpolationPoints; j++){

				// Compute x, y, and z distances between target i and interpolation point j
				dxt = xi - clusterX[clusterStart + j];
				dyt = yi - clusterY[clusterStart + j];
				dzt = zi - clusterZ[clusterStart + j];

				r = sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

//				tempPotential += (clusterM[clusterStart + j]-qi ) * exp(-kappa*r) / r;
				tempPotential += (clusterM[clusterStart + j]-qi*clusterM2[clusterStart + j] ) * exp(-kappa*r) / r;
	//			tempPotential += clusterM[clusterStart + j]* exp(-kappa*r) / r - qi*clusterM2[clusterStart + j] * exp(-kappa*r) / r;

							}

			EnP[batchStart + i] += tempPotential;
		}
		}

        
    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if ( (p->num_children == 0)|(smallEnoughLeaf==1) ) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct_yuk_SS(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, kappa, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
            	compute_pc_yuk_SS(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, kappa, EnP,
							clusterX, clusterY, clusterZ, clusterM, clusterM2);
            }
        }
    }

    return;

} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct_yuk_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;
    double d_peng, r;

# pragma acc kernels present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
    {
	#pragma acc loop independent
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
        d_peng = 0.0;
		#pragma acc loop independent
        for (i = ibeg - 1; i < iend; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            r = sqrt(tx*tx + ty*ty + tz*tz);
            if (r > DBL_MIN) {
            	d_peng += (qS[i] - qT[ii]) * wS[i] * exp(-kappa*r) / r;
            }
        }
        EnP[ii] += d_peng;
    }
    }

    return;

} /* END function pc_comp_direct_yuk */


