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

void fill_in_cluster_data_SS(struct particles *clusters, struct particles *sources, struct tnode *troot, int order){

	int pointsPerCluster = (order+1)*(order+1)*(order+1);
	int numInterpPoints = numnodes * pointsPerCluster;
	make_vector(clusters->x, numInterpPoints);
	make_vector(clusters->y, numInterpPoints);
	make_vector(clusters->z, numInterpPoints);
	make_vector(clusters->q, numInterpPoints);
	make_vector(clusters->w, numInterpPoints);  // will be used in singularity subtraction
	clusters->num=numInterpPoints;

	for (int i=0; i< numInterpPoints; i++){
		clusters->q[i]=0.0;
		clusters->w[i]=0.0;
	}

#pragma acc data copyin(tt[0:torderlim], \
		sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num]) \
		copy(clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num])



	addNodeToArray_SS(troot, sources, clusters, order, numInterpPoints, pointsPerCluster);

	return;
}

void addNodeToArray_SS(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;


	pc_comp_ms_modifiedF_SS(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
					clusters->x,clusters->y,clusters->z,clusters->q,clusters->w);


	for (i = 0; i < p->num_children; i++) {
		addNodeToArray_SS(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
	}

	return;
}




void pc_comp_ms_modifiedF_SS(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
		double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW){

	int i,j,k;
	int pointsPerCluster = torderlim*torderlim*torderlim;
	int pointsInNode = p->numpar;
	int startingIndexInClusters = p->node_index * pointsPerCluster;
	int startingIndexInSources = p->ibeg-1;

	double x0, x1, y0, y1, z0, z1;  // bounding box

	double weights[torderlim];
	double dj[torderlim];
	double *modifiedF, *modifiedF2;
	make_vector(modifiedF,pointsInNode);
	make_vector(modifiedF2,pointsInNode);

	double nodeX[torderlim], nodeY[torderlim], nodeZ[torderlim];

	int *exactIndX, *exactIndY, *exactIndZ;
	make_vector(exactIndX, pointsInNode);
	make_vector(exactIndY, pointsInNode);
	make_vector(exactIndZ, pointsInNode);


	// Set the bounding box.
//	x0 = p->x_min-1e-13*(p->x_max-p->x_min);
//	x1 = p->x_max+2e-13*(p->x_max-p->x_min);
//	y0 = p->y_min-1e-13*(p->y_max-p->y_min);
//	y1 = p->y_max+2e-13*(p->y_max-p->y_min);
//	z0 = p->z_min-1e-13*(p->z_max-p->z_min);
//	z1 = p->z_max+2e-13*(p->z_max-p->z_min);

	x0 = p->x_min;
	x1 = p->x_max;
	y0 = p->y_min;
	y1 = p->y_max;
	z0 = p->z_min;
	z1 = p->z_max;

	// Make and zero-out arrays to store denominator sums
	double sumX, sumY, sumZ;


#pragma acc kernels present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ, clusterW, tt) \
	create(modifiedF[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
			modifiedF2[0:pointsInNode],nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],weights[0:torderlim],dj[0:torderlim])
	{

	#pragma acc loop independent
	for (j=0;j<pointsInNode;j++){
		modifiedF[j] = qS[startingIndexInSources+j]*wS[startingIndexInSources+j];
		modifiedF2[j] = wS[startingIndexInSources+j];
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

			cx = sx - nodeX[j];
			cy = sy - nodeY[j];
			cz = sz - nodeZ[j];

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
		modifiedF2[i] /= denominator;

	}


	// Compute moments for each interpolation point
	double numerator, xn, yn, zn, temp, temp2;
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
		temp2 = 0.0;
		#pragma acc loop independent
		for (i=0;i<pointsInNode; i++){  // loop over source points
			sx = xS[startingIndexInSources+i];
			sy = yS[startingIndexInSources+i];
			sz = zS[startingIndexInSources+i];

			numerator=1.0;

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


//			numerator *=  w1 / (sx - cx);
//			numerator *=  w2 / (sy - cy);
//			numerator *=  w3 / (sz - cz);



			temp += numerator*modifiedF[i];
			temp2 += numerator*modifiedF2[i];


		}

		clusterQ[startingIndexInClusters + j] += temp;
		clusterW[startingIndexInClusters + j] += temp2;

	}


	}

	free_vector(modifiedF);
	free_vector(modifiedF2);
	free_vector(exactIndX);
	free_vector(exactIndY);
	free_vector(exactIndZ);


	return;
}

void pc_treecode_coulomb_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappa, double *tpeng, double *EnP, int numDevices, int numThreads)
{
    /* local variables */
    int i, j;
    double kappaSq=kappa*kappa;
    for (i = 0; i < targets->num; i++)
        EnP[i] = 2.0*M_PI*kappaSq*targets->q[i]; // change this to whatever it should be
    
#pragma omp parallel num_threads(numThreads)
	{
    	if (omp_get_thread_num()<numDevices){
    		acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
    	}
        int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
		if (this_thread==0){printf("numDevices: %i\n", numDevices);}
		if (this_thread==0){printf("num_threads: %i\n", num_threads);}

		double *EnP2;
		make_vector(EnP2,targets->num);
		for (i = 0; i < targets->num; i++)
			EnP2[i] = 0.0;

#pragma acc data copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num]) \
		copy(EnP2[0:targets->num])
    {

	#pragma omp for private(j) schedule(guided)
	for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc_coulomb_SS(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, kappaSq, EnP2,
				clusters->x, clusters->y, clusters->z, clusters->q, clusters->w);
        }
    }

	#pragma acc wait
    } // end acc data region

    for (int k = 0; k < targets->num; k++){
    	if (EnP2[k] != 0.0)
    		EnP[k] += EnP2[k];
    	}
    free_vector(EnP2);
	} // end omp parallel region

    *tpeng = sum(EnP, targets->num);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_coulomb_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP,
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


	int smallEnoughLeaf=0;
	if (p->numpar < torderlim*torderlim*torderlim){
		smallEnoughLeaf=1;
	}else{
		smallEnoughLeaf=0;
	}
	if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (smallEnoughLeaf==0)  ) {

	int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
	int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
	int clusterStart = numberOfInterpolationPoints*p->node_index;



	double xi,yi,zi,qi,r;
	int batchStart = batch_ind[0] - 1;

	int streamID = rand() % 4;
	# pragma acc kernels async(streamID) present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM, clusterM2)
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

			tempPotential += (clusterM[clusterStart + j]-qi*clusterM2[clusterStart + j]* exp(-r*r/kappaSq) )  / r;

						}
		#pragma acc atomic
		EnP[batchStart + i] += tempPotential;
	}
	}



//		EnP[batch_ind[0] - 1 + i] +=  ( Weights1[j] - Weights2[j]*qT[ batch_ind[0] - 1 + i]*exp(-r*r/kappaSq) ) / r
    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if ( (p->num_children == 0)|(smallEnoughLeaf==1) ) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct_coulomb_SS(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, kappaSq, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
            	compute_pc_coulomb_SS(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, kappaSq, EnP,
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
void pc_comp_direct_coulomb_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;
    double d_peng, r;

    int streamID = rand() % 4;
	# pragma acc kernels async(streamID) present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
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
            	d_peng += wS[i]* ( qS[i] - qT[ii]* exp(-r*r/kappaSq) )  / r;
            }
        }
		#pragma acc atomic
        EnP[ii] += d_peng;
    }
    }

    return;

} /* END function pc_comp_direct_yuk */


//d_peng += wS[i]* ( qS[i] - qT[ii]* exp(-r*r/kappaSq) )  / r;


void pc_interaction_list_treecode_Coulomb_SS(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng,double kappaSq, double *EnP, int numDevices, int numThreads)
{
	    int i, j;

	    for (i = 0; i < targets->num; i++)
//	        EnP[i] = 0.0;
	    	EnP[i] = 2.0*M_PI*kappaSq*targets->q[i];

	    printf("Using interaction lists!\n");

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
			double *wC = clusters->w;

//			printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//			printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

			int * ibegs = tree_array->ibeg;
			int * iends = tree_array->iend;

	#pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], qS[0:sources->num], wS[0:sources->num], \
			xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
			xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num], wC[0:clusters->num], tree_inter_list[0:numnodes*batches->num], \
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
		double xi,yi,zi,qi;

		int numberOfTargets;
		int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
		int clusterStart, batchStart;

		int numberOfClusterApproximations, numberOfDirectSums;

		int streamID;
		#pragma omp for private(j,ii,jj,batch_ibeg,batch_iend,numberOfClusterApproximations,numberOfDirectSums,numberOfTargets,batchStart,node_index,clusterStart,streamID) schedule(guided)
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
					qi = qT[ batchStart + ii];

					for (jj = 0; jj < numberOfInterpolationPoints; jj++){

						// Compute x, y, and z distances between target i and interpolation point j
						dxt = xi - xC[clusterStart + jj];
						dyt = yi - yC[clusterStart + jj];
						dzt = zi - zC[clusterStart + jj];
						r = sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
//						tempPotential += qC[clusterStart + jj] *exp(-kappa*r)/ r;
						tempPotential += (qC[clusterStart + jj]-qi*wC[clusterStart + jj]* exp(-r*r/kappaSq) )  / r;
//						tempPotential += (clusterM[clusterStart + j]-qi*clusterM2[clusterStart + j]* exp(-r*r/kappaSq) )  / r;
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
		            r = sqrt(tx*tx + ty*ty +  tz*tz);
		            if (r > DBL_MIN) {
//		            	d_peng += qS[jj] * wS[jj] *exp(-kappa*r)/ r;
		            	d_peng += wS[jj]*( qS[jj] - qT[ii]* exp(-r*r/kappaSq) )/r;
//		            	d_peng += wS[i]* ( qS[i] - qT[ii]* exp(-r*r/kappaSq) )  / r;

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
	    } // end acc data region


	    for (int k = 0; k < targets->num; k++){
	    	if (EnP2[k] != 0.0) EnP[k] += EnP2[k];
	    	if (EnP2[k] != 0.0)EnP[k] += EnP3[k];
			}

	    free_vector(EnP2);
	    free_vector(EnP3);
		} // end omp parallel region

	    printf("Exited the main comp_pc call.\n");
	    *tpeng = sum(EnP, targets->num);

	    return;

	} /* END of function pc_treecode */


