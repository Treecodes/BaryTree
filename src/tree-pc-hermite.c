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


void fill_in_cluster_data_hermite(struct particles *clusters, struct particles *sources, struct tnode *troot, int order){

	int pointsPerCluster = (order+1)*(order+1)*(order+1);
	int numInterpPoints = numnodes * pointsPerCluster;
	make_vector(clusters->x, numInterpPoints);
	make_vector(clusters->y, numInterpPoints);
	make_vector(clusters->z, numInterpPoints);
	make_vector(clusters->q, numInterpPoints);
	make_vector(clusters->qx, numInterpPoints);
	make_vector(clusters->qy, numInterpPoints);
	make_vector(clusters->qz, numInterpPoints);
	make_vector(clusters->qxy, numInterpPoints);
	make_vector(clusters->qyz, numInterpPoints);
	make_vector(clusters->qxz, numInterpPoints);
	make_vector(clusters->qxyz, numInterpPoints);
	make_vector(clusters->w, numInterpPoints);  // will be used in singularity subtraction
//	memset(clusters->q,0,numInterpPoints*sizeof(double));
	clusters->num=numInterpPoints;

	for (int i=0; i< numInterpPoints; i++){
		clusters->w[i]=0.0;
	}
	for (int i=0; i< numInterpPoints; i++){
		clusters->q[i]=0.0;
		clusters->qx[i]=0.0;
		clusters->qy[i]=0.0;
		clusters->qz[i]=0.0;
		clusters->qxy[i]=0.0;
		clusters->qyz[i]=0.0;
		clusters->qxz[i]=0.0;
		clusters->qxyz[i]=0.0;
	}

#pragma acc data copy(tt[0:torderlim],ww[0:torderlim], \
		sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num], \
		clusters->qx[0:clusters->num],clusters->qy[0:clusters->num],clusters->qz[0:clusters->num],clusters->qxy[0:clusters->num],clusters->qyz[0:clusters->num], \
		clusters->qxz[0:clusters->num], clusters->qxyz[0:clusters->num])

	addNodeToArray_hermite(troot, sources, clusters, order, numInterpPoints, pointsPerCluster);

	return;
}




void pc_comp_ms_modifiedF_hermite(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
		double *clusterX, double *clusterY, double *clusterZ,
		double *clusterQ, double *clusterQx, double *clusterQy, double *clusterQz, double *clusterQxy,
		double *clusterQyz, double *clusterQxz, double *clusterQxyz){

	int i,j,k;
	int pointsPerCluster = torderlim*torderlim*torderlim;
	int pointsInNode = p->numpar;
	int startingIndexInClusters = p->node_index * pointsPerCluster;
	int startingIndexInSources = p->ibeg-1;

	double x0, x1, y0, y1, z0, z1;  // bounding box


//	double weights[torderlim];
	double dj[torderlim],wx[torderlim],wy[torderlim],wz[torderlim];
	double *modifiedF;
	make_vector(modifiedF,pointsInNode);

	double nodeX[torderlim], nodeY[torderlim], nodeZ[torderlim];

	int *exactIndX, *exactIndY, *exactIndZ;
	make_vector(exactIndX, pointsInNode);
	make_vector(exactIndY, pointsInNode);
	make_vector(exactIndZ, pointsInNode);

	// Set the bounding box.

	x0 = p->x_min;
	x1 = p->x_max;
	y0 = p->y_min;
	y1 = p->y_max;
	z0 = p->z_min;
	z1 = p->z_max;

	// Make and zero-out arrays to store denominator sums
	double sumX, sumY, sumZ;

	int streamID = rand() % 2;
#pragma acc kernels present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ,tt,ww, \
		clusterQ,clusterQx,clusterQy,clusterQz,clusterQxy,clusterQyz,clusterQxz,clusterQxyz) \
	create(modifiedF[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
			nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],dj[0:torderlim], \
			wx[0:torderlim],wy[0:torderlim],wz[0:torderlim])
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

	#pragma acc loop independent
	for (j = 0; j < torderlim;j++){
		dj[j] = 1.0;
		wx[j] = -4.0 * ww[j] / (x1 - x0);
		wy[j] = -4.0 * ww[j] / (y1 - y0);
		wz[j] = -4.0 * ww[j] / (z1 - z0);
	}
	dj[0] = 0.25;
	dj[torder] = 0.25;



	// Compute modified f values
	double sx,sy,sz,dx,dy,dz,denominator,w;



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

			dx = sx-nodeX[j];
			dy = sy-nodeY[j];
			dz = sz-nodeZ[j];

			if (fabs(dx)<DBL_MIN) exactIndX[i]=j;
			if (fabs(dy)<DBL_MIN) exactIndY[i]=j;
			if (fabs(dz)<DBL_MIN) exactIndZ[i]=j;

			// Increment the sums
//			w = weights[j];
			sumX += dj[j] / (dx*dx) + wx[j]/dx;
			sumY += dj[j] / (dy*dy) + wy[j]/dy;
			sumZ += dj[j] / (dz*dz) + wz[j]/dz;

		}

//		denominator = sumX*sumY*sumZ;

		denominator = 1.0;
		if (exactIndX[i]==-1) denominator *= sumX;
		if (exactIndY[i]==-1) denominator *= sumY;
		if (exactIndZ[i]==-1) denominator *= sumZ;

		modifiedF[i] /= denominator;
//		printf("%1.2e \n", modifiedF[i]);
	}


	// Compute moments for each interpolation point
	double numerator0,numerator1,numerator2,numerator3,numerator4,numerator5,numerator6,numerator7,  xn, yn, zn;
	double Ax,Ay,Az,Bx,By,Bz;
	double temp0,temp1,temp2,temp3,temp4,temp5,temp6,temp7;
	int k1, k2, k3, kk,sourcePointIndex,interpolationPointIndex;
	double cx,cy,cz;

	#pragma acc loop independent
	for (j=0;j<pointsPerCluster;j++){ // loop over interpolation points, set (dx,dy,dz) for this point
		// compute k1, k2, k3 from j
		k1 = j%torderlim;
		kk = (j-k1)/torderlim;
		k2 = kk%torderlim;
		kk = kk - k2;
		k3 = kk / torderlim;

		cz = nodeZ[k3];

		cy = nodeY[k2];

		cx = nodeX[k1];


		// Fill cluster X, Y, and Z arrays
		interpolationPointIndex = startingIndexInClusters + j;
		clusterX[interpolationPointIndex] = cx;
		clusterY[interpolationPointIndex] = cy;
		clusterZ[interpolationPointIndex] = cz;


		// Increment cluster Q array
		temp0 = 0.0; temp1=0.0; temp2 = 0.0; temp3=0.0;
		temp4 = 0.0; temp5=0.0; temp6 = 0.0; temp7=0.0;
		#pragma acc loop independent
		for (i=0;i<pointsInNode; i++){  // loop over source points
			sourcePointIndex = startingIndexInSources+i;
			dx = xS[sourcePointIndex]-cx;
			dy = yS[sourcePointIndex]-cy;
			dz = zS[sourcePointIndex]-cz;

			numerator0=1.0;numerator1=1.0;numerator2=1.0;numerator3=1.0;numerator4=1.0;numerator5=1.0;numerator6=1.0;numerator7=1.0;

			Ax = ( dj[k1]/(dx*dx) + wx[k1]/dx );
			Ay = ( dj[k2]/(dy*dy) + wy[k2]/dy );
			Az = ( dj[k3]/(dz*dz) + wz[k3]/dz );
			Bx = ( dj[k1]/(dx) );
			By = ( dj[k2]/(dy) );
			Bz = ( dj[k3]/(dz) );



			if (exactIndX[i]==-1){
				numerator0 *=  Ax; 					// Aaa

				numerator1 *=  Bx; 					// Baa
				numerator2 *=  Ax; 					// Aba
				numerator3 *=  Ax; 					// Aab

				numerator4 *=  Bx; 					// Bba
				numerator5 *=  Ax; 					// Abb
				numerator6 *=  Bx; 					// Bab

				numerator7 *=  Bx; 					// Bbb

			}else{
				if (exactIndX[i]!=k1){ numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0; numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
				}else{
					numerator1 *= 0; numerator4 *= 0;numerator6 *= 0; numerator7 *= 0;
				}
			}

			if (exactIndY[i]==-1){
				numerator0 *=  Ay;					// aAa

				numerator1 *=  Ay; 					// bAa
				numerator2 *=  By; 					// aBa
				numerator3 *=  Ay; 					// aAb

				numerator4 *=  By; 					// bBa
				numerator5 *=  By; 					// aBb
				numerator6 *=  Ay; 					// bAb

				numerator7 *=  By; 					// bBb

			}else{
				if (exactIndY[i]!=k2) { numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0; numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
				}else{
					numerator2 *= 0; numerator4 *= 0;numerator5 *= 0; numerator7 *= 0;
				}

			}

			if (exactIndZ[i]==-1){
				numerator0 *=  Az;					// aaA

				numerator1 *=  Az;					// baA
				numerator2 *=  Az;					// abA
				numerator3 *=  Bz;					// aaB

				numerator4 *=  Az; 					// bbA
				numerator5 *=  Bz; 					// abB
				numerator6 *=  Bz; 					// baB

				numerator7 *=  Bz; 					// bbB

			}else{
				if (exactIndZ[i]!=k3) { numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0; numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
				}else{
					numerator3 *= 0; numerator5 *= 0;numerator6 *= 0; numerator7 *= 0;
				}
			}





			temp0 += numerator0*modifiedF[i];
			temp1 += numerator1*modifiedF[i];
			temp2 += numerator2*modifiedF[i];
			temp3 += numerator3*modifiedF[i];
			temp4 += numerator4*modifiedF[i];
			temp5 += numerator5*modifiedF[i];
			temp6 += numerator6*modifiedF[i];
			temp7 += numerator7*modifiedF[i];


		}

		clusterQ[interpolationPointIndex] += temp0;
		clusterQx[interpolationPointIndex] += temp1;
		clusterQy[interpolationPointIndex] += temp2;
		clusterQz[interpolationPointIndex] += temp3;
		clusterQxy[interpolationPointIndex] += temp4;
		clusterQyz[interpolationPointIndex] += temp5;
		clusterQxz[interpolationPointIndex] += temp6;
		clusterQxyz[interpolationPointIndex] += temp7;

//		printf("\n\n %1.2e\n", clusterQ[8*(startingIndexInClusters + j) + 0]);
//		printf("%1.2e\n", clusterQ[8*(startingIndexInClusters + j) + 1]);
//		printf("%1.2e\n", clusterQ[8*(startingIndexInClusters + j) + 2]);
//		printf("%1.2e\n", clusterQ[8*(startingIndexInClusters + j) + 3]);
//		printf("%1.2e\n", clusterQ[8*(startingIndexInClusters + j) + 5]);
//		printf("%1.2e\n\n", clusterQ[8*(startingIndexInClusters + j) + 7]);

	}


	}

	free_vector(modifiedF);
	free_vector(exactIndX);
	free_vector(exactIndY);
	free_vector(exactIndZ);


	return;
}


void addNodeToArray_hermite(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;


//	if (torderlim*torderlim*torderlim < p->numpar){

	pc_comp_ms_modifiedF_hermite(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
			clusters->x,clusters->y,clusters->z,clusters->q,clusters->qx,clusters->qy,clusters->qz,clusters->qxy, \
			clusters->qyz, clusters->qxz, clusters->qxyz);

	p->exist_ms = 1;



	for (i = 0; i < p->num_children; i++) {
		addNodeToArray_hermite(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
	}

	return;
}




void pc_interaction_list_treecode_hermite_coulomb(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double *EnP)
{
	    int i, j;

	    for (i = 0; i < targets->num; i++)
	        EnP[i] = 0.0;




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
			double *qxC = clusters->qx;
			double *qyC = clusters->qy;
			double *qzC = clusters->qz;
			double *qxyC = clusters->qxy;
			double *qyzC = clusters->qyz;
			double *qxzC = clusters->qxz;
			double *qxyzC = clusters->qxyz;

//			printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//			printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

			int * ibegs = tree_array->ibeg;
			int * iends = tree_array->iend;

	#pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], qS[0:sources->num], wS[0:sources->num], \
			xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
			xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num],tree_inter_list[0:numnodes*batches->num], \
			qxC[0:clusters->num],qyC[0:clusters->num],qzC[0:clusters->num],qxyC[0:clusters->num], \
			qyzC[0:clusters->num],qxzC[0:clusters->num],qxyzC[0:clusters->num],qC[0:clusters->num], \
			direct_inter_list[0:batches->num * numleaves], ibegs[0:numnodes], iends[0:numnodes]) \
			copy(EnP3[0:targets->num], EnP2[0:targets->num])

//#pragma acc data copyin(targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
//		sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
//		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], \
//		clusters->qx[0:clusters->num],clusters->qy[0:clusters->num],clusters->qz[0:clusters->num],clusters->qxy[0:clusters->num],clusters->qyz[0:clusters->num], \
//		clusters->qxz[0:clusters->num], clusters->qxyz[0:clusters->num]) \
//		copy(EnP3[0:targets->num], EnP2[0:targets->num])
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
		double rinv,r3inv,r5inv,r7inv;

		int numberOfTargets;
		int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
		int clusterStart, batchStart, sourceIdx;

		int numberOfClusterApproximations, numberOfDirectSums;

		int streamID;
//		#pragma omp for private(j,ii,jj,sourceIdx,batch_ibeg,batch_iend,numberOfClusterApproximations,numberOfDirectSums,numberOfTargets,batchStart,node_index,clusterStart,streamID,rinv,r3inv,r5inv,r7inv)
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
						sourceIdx = clusterStart + jj;
						// Compute x, y, and z distances between target i and interpolation point j
						dxt = xi - xC[sourceIdx];
						dyt = yi - yC[sourceIdx];
						dzt = zi - zC[sourceIdx];
//						tempPotential += qC[clusterStart + jj] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

						rinv = 1 / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
						r3inv = rinv*rinv*rinv;
						r5inv = r3inv*rinv*rinv;
						r7inv = r5inv*rinv*rinv;



						tempPotential +=       rinv  * ( qC[sourceIdx])
														+      r3inv * ( qxC[sourceIdx]*dxt +  qyC[sourceIdx]*dyt +  qzC[sourceIdx]*dzt )
														+ 3 *  r5inv * ( qxyC[sourceIdx]*dxt*dyt +  qyzC[sourceIdx]*dyt*dzt +  qxzC[sourceIdx]*dxt*dzt )
														+ 15 * r7inv *   qxyzC[sourceIdx]*dxt*dyt*dzt
														;



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
	    } // end acc data region


	    for (int k = 0; k < targets->num; k++){
	    	if (EnP2[k] != 0.0)
//			#pragma omp critical
//	    	{
				EnP[k] += EnP2[k];
				EnP[k] += EnP3[k];
//	    	} // end omp critical
			}

	    free_vector(EnP2);
	    free_vector(EnP3);
//		} // end omp parallel region

	    printf("Exited the main comp_pc call.\n");
	    *tpeng = sum(EnP, targets->num);

	    return;

	} /* END of function pc_treecode */

