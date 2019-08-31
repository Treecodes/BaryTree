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


void pc_interaction_list_treecode_dcf(struct tnode_array *tree_array,
                                      struct particles *clusters, struct batch *batches,
                                      int *tree_inter_list, int *direct_inter_list,
                                      struct particles *sources, struct particles *targets,
                                      double *tpeng, double eta, double *EnP,
                                      int numDevices, int numThreads)
{
	    int i, j;

	    for (i = 0; i < targets->num; i++)
	        EnP[i] = 0.0;

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
						r = sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
						tempPotential += qC[clusterStart + jj] * erf(r / eta) / r;

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
		            	d_peng += qS[jj] * wS[jj] * erf(r / eta) / r;
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

