/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void pc_interaction_list_treecode_yuk_SS(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double kappa, double *EnP)
{
        int i, j;
        int rank; int numProcs;	int ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

//        for (i = 0; i < targets->num; i++)
//            EnP[i] = 0.0;

        for (i = 0; i < targets->num; i++){
			EnP[i] = 4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
		}



            double *EnP2, *EnP3;
            make_vector(EnP2,targets->num);
            make_vector(EnP3,targets->num);

            for (i = 0; i < targets->num; i++) {
                EnP3[i] = 0.0;
                EnP2[i] = 0.0;
            }

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

//          printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//          printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

            int * ibegs = tree_array->ibeg;
            int * iends = tree_array->iend;

            int * clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], \
                            qS[0:sources->num], wS[0:sources->num], \
                            xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
                            xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num],wC[0:clusters->num], \
                            tree_inter_list[0:tree_numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) copy(EnP3[0:targets->num], EnP2[0:targets->num])
#endif
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

        for (i = 0; i < batches->num; i++) {
            batch_ibeg = batches->index[i][0];
            batch_iend = batches->index[i][1];
            numberOfClusterApproximations = batches->index[i][2];
            numberOfDirectSums = batches->index[i][3];

//            printf("Rank %i, batch %i, number of cluster approximations: %i\n", rank, i, numberOfClusterApproximations);
//            printf("Rank %i, batch %i, number of direct interactions: %i\n", rank, i, numberOfDirectSums);

            numberOfTargets = batch_iend - batch_ibeg + 1;
            batchStart =  batch_ibeg - 1;

            for (j = 0; j < numberOfClusterApproximations; j++) {
                node_index = tree_inter_list[i * tree_numnodes + j];
//                clusterStart = numberOfInterpolationPoints*node_index;
                clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                streamID = j%3;
#ifdef OPENACC_ENABLED
                #pragma acc kernels async(streamID) //present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
                {
                #pragma acc loop independent
#endif
                for (ii = 0; ii < numberOfTargets; ii++) {
                    tempPotential = 0.0;
                    xi = xT[ batchStart + ii];
                    yi = yT[ batchStart + ii];
                    zi = zT[ batchStart + ii];
                    qi = qT[ batchStart + ii];;

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
                        // Compute x, y, and z distances between target i and interpolation point j
                        dxt = xi - xC[clusterStart + jj];
                        dyt = yi - yC[clusterStart + jj];
                        dzt = zi - zC[clusterStart + jj];
                        r=sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
//                        tempPotential += qC[clusterStart + jj]*exp(-kappa*r) / r;
                        tempPotential += ( qC[clusterStart + jj]-qi*wC[clusterStart + jj] )*exp(-kappa*r) / r;
//                        (qS[i] - qT[ii])

                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
#endif
                    EnP3[batchStart + ii] += tempPotential;
                }
#ifdef OPENACC_ENABLED
                } // end kernel
#endif
            } // end for loop over cluster approximations

            for (j = 0; j < numberOfDirectSums; j++) {

                node_index = direct_inter_list[i * tree_numnodes + j];

                source_start=ibegs[node_index]-1;
                source_end=iends[node_index];

                streamID = j%3;
#ifdef OPENACC_ENABLED
                # pragma acc kernels async(streamID)
                {
                #pragma acc loop independent
#endif
                for (ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                    d_peng = 0.0;

                    for (jj = source_start; jj < source_end; jj++) {
                        tx = xS[jj] - xT[ii];
                        ty = yS[jj] - yT[ii];
                        tz = zS[jj] - zT[ii];
                        r = sqrt(tx*tx + ty*ty + tz*tz);

                        if (r > DBL_MIN) {
                            d_peng += ( qS[jj] - qT[ii] ) * wS[jj]*exp(-kappa*r) / r;
                        }
                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
#endif
                    EnP2[ii] += d_peng;
                }
#ifdef OPENACC_ENABLED
                } // end kernel
#endif
            } // end loop over number of leaves
        } // end loop over target batches
#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
        } // end acc data region

        double totalDueToApprox = 0.0, totalDueToDirect = 0.0;
        totalDueToApprox = sum(EnP3, targets->num);
        totalDueToDirect = sum(EnP2, targets->num);
//        printf("Potential due to approximations: %f\n",totalDueToApprox);
//        printf("Potential due to direct: %f\n",totalDueToDirect);
        for (int k = 0; k < targets->num; k++) {
            if (EnP2[k] != 0.0)
                EnP[k] += EnP2[k];
                EnP[k] += EnP3[k];
            }

            free_vector(EnP2);
            free_vector(EnP3);

        *tpeng = sum(EnP, targets->num);

        return;

    } /* END of function pc_treecode */

