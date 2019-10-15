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


void fill_in_cluster_data_SS(struct particles *clusters, struct particles *sources, struct tnode *troot, int order, struct tnode_array * tree_array){

	int tree_numnodes = tree_array->numnodes;
    int pointsPerCluster = (order+1)*(order+1)*(order+1);
    int numInterpPoints = tree_numnodes * pointsPerCluster;
    make_vector(clusters->x, numInterpPoints);
    make_vector(clusters->y, numInterpPoints);
    make_vector(clusters->z, numInterpPoints);
    make_vector(clusters->q, numInterpPoints);
    make_vector(clusters->w, numInterpPoints);  // will be used in singularity subtraction
    clusters->num=numInterpPoints;

    for (int i = 0; i < numInterpPoints; i++) {
        clusters->x[i]=0.0;
        clusters->y[i]=0.0;
        clusters->z[i]=0.0;
        clusters->q[i]=0.0;
        clusters->w[i]=0.0;
    }

        double *tempQ, *tempX, *tempY, *tempZ, *tempW;
        make_vector(tempX,clusters->num);
        make_vector(tempY,clusters->num);
        make_vector(tempZ,clusters->num);
        make_vector(tempQ,clusters->num);
        make_vector(tempW,clusters->num);

        for (int i = 0; i < clusters->num; i++) {
            tempX[i] = 0.0;
            tempY[i] = 0.0;
            tempZ[i] = 0.0;
            tempQ[i] = 0.0;
            tempW[i] = 0.0;
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
        double *qW = clusters->w;

        int clusterNum = clusters->num;
        int sourceNum = sources->num;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(tt[0:torderlim], \
        xS[0:sourceNum], yS[0:sourceNum], zS[0:sourceNum], qS[0:sourceNum], wS[0:sourceNum]) \
        copy(tempX[0:clusterNum], tempY[0:clusterNum], tempZ[0:clusterNum], tempQ[0:clusterNum], tempW[0:clusterNum])
        {
#endif
            for (int i = 0; i < tree_numnodes; i++) {  // start from i=1, don't need to compute root moments
                pc_comp_ms_modifiedF_SS(tree_array, i, xS, yS, zS, qS, wS,
                                     tempX, tempY, tempZ, tempQ, tempW);
            }
#ifdef OPENACC_ENABLED
            #pragma acc wait
        } // end ACC DATA REGION
#endif

        int counter=0;
        for (int j = 0; j < clusters->num; j++)
        {
            if (tempQ[j]!=0.0){
                clusters->x[j] = tempX[j];
                clusters->y[j] = tempY[j];
                clusters->z[j] = tempZ[j];
                clusters->q[j] += tempQ[j];
                clusters->w[j] += tempW[j];
            }
        } // end j loop

        free_vector(tempX);
        free_vector(tempY);
        free_vector(tempZ);
        free_vector(tempQ);
        free_vector(tempW);


    return;
}



//void addNodeToArray_SS(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster)
//{
//	int torderlim = order+1;
//	int startingIndex = p->node_index * pointsPerCluster;
//	int i;
//
//
////	pc_comp_ms_modifiedF_SS(p, sources->x, sources->y, sources->z, sources->q, sources->w, \
////					clusters->x,clusters->y,clusters->z,clusters->q,clusters->w);
//
//
//	for (i = 0; i < p->num_children; i++) {
//		addNodeToArray_SS(p->child[i],sources,clusters,order,numInterpPoints,pointsPerCluster);
//	}
//
//	return;
//}




void pc_comp_ms_modifiedF_SS(struct tnode_array * tree_array, int idx,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{
    int i,j,k;
    int pointsPerCluster = torderlim*torderlim*torderlim;
    int pointsInNode = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClusters = idx * pointsPerCluster;
    int startingIndexInSources = tree_array->ibeg[idx]-1;

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
#ifdef OPENACC_ENABLED
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ, clusterW,tt) \
    create(modifiedF[0:pointsInNode],modifiedF2[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
            nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],weights[0:torderlim],dj[0:torderlim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < pointsInNode; j++) {
        modifiedF[j] = qS[startingIndexInSources+j] * wS[startingIndexInSources+j];
        modifiedF2[j] = wS[startingIndexInSources+j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (i = 0; i < torderlim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < torder+1; j++){
        dj[j] = 1.0;
        if (j==0) dj[j] = 0.5;
        if (j==torder) dj[j]=0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < torderlim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (i = 0; i < pointsInNode; i++) { // loop through the source points

        sumX=0.0;
        sumY=0.0;
        sumZ=0.0;

        sx = xS[startingIndexInSources+i];
        sy = yS[startingIndexInSources+i];
        sz = zS[startingIndexInSources+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (j = 0; j < torderlim; j++) {  // loop through the degree

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

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < pointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
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
        temp2=0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (i = 0; i < pointsInNode; i++) {  // loop over source points
            sx = xS[startingIndexInSources+i];
            sy = yS[startingIndexInSources+i];
            sz = zS[startingIndexInSources+i];

            numerator=1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (exactIndX[i]==-1) {
                numerator *=  w1 / (sx - cx);
            } else {
                if (exactIndX[i]!=k1) numerator *= 0;
            }

            if (exactIndY[i]==-1) {
                numerator *=  w2 / (sy - cy);
            } else {
                if (exactIndY[i]!=k2) numerator *= 0;
            }

            if (exactIndZ[i]==-1) {
                numerator *=  w3 / (sz - cz);
            } else {
                if (exactIndZ[i]!=k3) numerator *= 0;
            }

            temp += numerator * modifiedF[i];
            temp2 += numerator*modifiedF2[i];
        }
        clusterQ[startingIndexInClusters + j] += temp;
		clusterW[startingIndexInClusters + j] += temp2;
    }
#ifdef OPENACC_ENABLED
    }
#endif

    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}




void pc_interaction_list_treecode_Coulomb_SS(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double kappa, double *EnP)
{
        int i, j;
        int rank; int numProcs;	int ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

        double kappaSq=kappa*kappa;
        for (i = 0; i < targets->num; i++){
			EnP[i] = 2.0*M_PI*kappaSq*targets->q[i];
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
                            xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num], \
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
                    qi = qT[ batchStart + ii];

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
                        // Compute x, y, and z distances between target i and interpolation point j
                        dxt = xi - xC[clusterStart + jj];
                        dyt = yi - yC[clusterStart + jj];
                        dzt = zi - zC[clusterStart + jj];
                        r=sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
//                        tempPotential += qC[clusterStart + jj]*exp(-kappa*r) / r;
                        tempPotential += ( qC[clusterStart + jj]- qi*wC[clusterStart + jj]*exp(-r*r/kappaSq) ) / r;
//            			tempPotential += (clusterM[clusterStart + j]-qi*clusterM2[clusterStart + j]* exp(-r*r/kappaSq) )  / r;

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
                            d_peng += ( qS[jj] - qT[ii]*exp(-r*r/kappaSq) ) * wS[jj] / r;
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




