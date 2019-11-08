/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "kernels/kernels.h"
#include "tree.h"




void pc_interaction_list_treecode(struct tnode_array *tree_array, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  double *xS, double *yS, double *zS, double *qS, double *wS,
                                  double *xT, double *yT, double *zT, double *qT,
                                  double *xC, double *yC, double *zC, double *qC, double *wC,
                                  double *totalPotential, double *pointwisePotential, int interpolationOrder,
                                  int numSources, int numTargets, int numClusters,
                                  int batch_approx_offset, int batch_direct_offset,
                                  char *kernelName, double kappa)
{
        int rank, numProcs, ierr;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

        double *potentialDueToDirect, *potentialDueToApprox;
        make_vector(potentialDueToDirect, numTargets);
        make_vector(potentialDueToApprox, numTargets);

        for (int i = 0; i < numTargets; i++) {
            potentialDueToApprox[i] = 0.0;
            potentialDueToDirect[i] = 0.0;
        }

        int *ibegs = tree_array->ibeg;
        int *iends = tree_array->iend;
        int *clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], qT[0:numTargets], \
                            xC[0:numClusters], yC[0:numClusters], zC[0:numClusters], qC[0:numClusters], \
                            tree_inter_list[0:batch_approx_offset*batches->num], \
                            direct_inter_list[0:batch_direct_offset*batches->num], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) \
                            copy(potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
        {

        int numberOfInterpolationPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);

        for (int i = 0; i < batches->num; i++) {
            int batch_ibeg = batches->index[i][0];
            int batch_iend = batches->index[i][1];
            int numberOfClusterApproximations = batches->index[i][2];
            int numberOfDirectSums = batches->index[i][3];

            int numberOfTargets = batch_iend - batch_ibeg + 1;
            int batchStart =  batch_ibeg - 1;



/**********************************************************/
/************** POTENTIAL FROM APPROX *********************/
/**********************************************************/

            for (int j = 0; j < numberOfClusterApproximations; j++) {
                int node_index = tree_inter_list[i * batch_approx_offset + j];
                int clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                int streamID = j%3;



        /***********************************************/
        /***************** Coulomb *********************/
        /***********************************************/

                if (strcmp(kernelName, "coulomb") == 0) {
                
#ifdef OPENACC_ENABLED
                    #pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = 0; ii < numberOfTargets; ii++) {
                        double tempPotential = 0.0;
                        double xi = xT[ batchStart + ii];
                        double yi = yT[ batchStart + ii];
                        double zi = zT[ batchStart + ii];
                        double qi = qT[ batchStart + ii];

                        for (int jj = 0; jj < numberOfInterpolationPoints; jj++) {
                            tempPotential += coulombKernel(xi, yi, zi, qi,
                                                           xC[clusterStart + jj], yC[clusterStart + jj],
                                                           zC[clusterStart + jj], qC[clusterStart + jj],
                                                           wC[clusterStart + jj], kappa);

                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToApprox[batchStart + ii] += tempPotential;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /***************** Yukawa **********************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
                    #pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = 0; ii < numberOfTargets; ii++) {
                        double tempPotential = 0.0;
                        double xi = xT[ batchStart + ii];
                        double yi = yT[ batchStart + ii];
                        double zi = zT[ batchStart + ii];
                        double qi = qT[ batchStart + ii];

                        for (int jj = 0; jj < numberOfInterpolationPoints; jj++) {
                            tempPotential += yukawaKernel(xi, yi, zi, qi,
                                                          xC[clusterStart + jj], yC[clusterStart + jj],
                                                          zC[clusterStart + jj], qC[clusterStart + jj],
                                                          wC[clusterStart + jj], kappa);

                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToApprox[batchStart + ii] += tempPotential;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /************* Coulomb with SS *****************/
        /***********************************************/

                } else if (strcmp(kernelName, "coulomb_SS") == 0) {

#ifdef OPENACC_ENABLED
                    #pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = 0; ii < numberOfTargets; ii++) {
                        double tempPotential = 0.0;
                        double xi = xT[ batchStart + ii];
                        double yi = yT[ batchStart + ii];
                        double zi = zT[ batchStart + ii];
                        double qi = qT[ batchStart + ii];

                        for (int jj = 0; jj < numberOfInterpolationPoints; jj++) {
                            tempPotential += coulombKernel_SS_approx(xi, yi, zi, qi,
                                                          xC[clusterStart + jj], yC[clusterStart + jj],
                                                          zC[clusterStart + jj], qC[clusterStart + jj],
                                                          wC[clusterStart + jj], kappa);

                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToApprox[batchStart + ii] += tempPotential;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /************** Yukawa with SS *****************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa_SS") == 0) {

#ifdef OPENACC_ENABLED
                    #pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = 0; ii < numberOfTargets; ii++) {
                        double tempPotential = 0.0;
                        double xi = xT[ batchStart + ii];
                        double yi = yT[ batchStart + ii];
                        double zi = zT[ batchStart + ii];
                        double qi = qT[ batchStart + ii];

                        for (int jj = 0; jj < numberOfInterpolationPoints; jj++) {
                            tempPotential += yukawaKernel_SS_approx(xi, yi, zi, qi,
                                                          xC[clusterStart + jj], yC[clusterStart + jj],
                                                          zC[clusterStart + jj], qC[clusterStart + jj],
                                                          wC[clusterStart + jj], kappa);

                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToApprox[batchStart + ii] += tempPotential;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif
                } // end kernel selection

            } // end loop over cluster approximations



/**********************************************************/
/************** POTENTIAL FROM DIRECT *********************/
/**********************************************************/

            for (int j = 0; j < numberOfDirectSums; j++) {

                int node_index = direct_inter_list[i * batch_direct_offset + j];

                int source_start = ibegs[node_index]-1;
                int source_end = iends[node_index];

                int streamID = j%3;



        /***********************************************/
        /***************** Coulomb *********************/
        /***********************************************/

                if (strcmp(kernelName, "coulomb") == 0) {

#ifdef OPENACC_ENABLED
                    # pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                        double d_peng = 0.0;

                        for (int jj = source_start; jj < source_end; jj++) {
                            d_peng += coulombKernel(xT[ii], yT[ii], zT[ii], qT[ii],
                                                    xS[jj], yS[jj], zS[jj], qS[jj], wS[jj], kappa);
                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToDirect[ii] += d_peng;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /***************** Yukawa **********************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
                    # pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                        double d_peng = 0.0;

                        for (int jj = source_start; jj < source_end; jj++) {
                            d_peng += yukawaKernel(xT[ii], yT[ii], zT[ii], qT[ii],
                                                   xS[jj], yS[jj], zS[jj], qS[jj], wS[jj], kappa);
                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToDirect[ii] += d_peng;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /************* Coulomb with SS *****************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
                    # pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                        double d_peng = 0.0;

                        for (int jj = source_start; jj < source_end; jj++) {
                            d_peng += coulombKernel_SS_direct(xT[ii], yT[ii], zT[ii], qT[ii],
                                                              xS[jj], yS[jj], zS[jj], qS[jj], wS[jj], kappa);
                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToDirect[ii] += d_peng;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif



        /***********************************************/
        /************** Yukawa with SS *****************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
                    # pragma acc kernels async(streamID)
                    {
                    #pragma acc loop independent
#endif
                    for (int ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                        double d_peng = 0.0;

                        for (int jj = source_start; jj < source_end; jj++) {
                            d_peng += yukawaKernel_SS_direct(xT[ii], yT[ii], zT[ii], qT[ii],
                                                             xS[jj], yS[jj], zS[jj], qS[jj], wS[jj], kappa);
                        }
#ifdef OPENACC_ENABLED
                        #pragma acc atomic
#endif
                        potentialDueToDirect[ii] += d_peng;
                    }
#ifdef OPENACC_ENABLED
                    } // end kernel
#endif
                } // end kernel selection

            } // end loop over number of leaves

        } // end loop over target batches
#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
        } // end acc data region


        double totalDueToApprox = sum(potentialDueToApprox, numTargets);
        double totalDueToDirect = sum(potentialDueToDirect, numTargets);

        printf("Total due to direct = %f\n", totalDueToDirect);
        printf("Total due to approx = %f\n", totalDueToApprox);

        for (int k = 0; k < numTargets; k++) {
            pointwisePotential[k] += potentialDueToDirect[k];
            pointwisePotential[k] += potentialDueToApprox[k];
         }

        free_vector(potentialDueToDirect);
        free_vector(potentialDueToApprox);

        *totalPotential = sum(pointwisePotential, numTargets);

        return;

    } /* END of function pc_treecode */
