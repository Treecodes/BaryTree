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
    #pragma omp critical 
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


    #pragma omp critical
    {
        if (maxlevel < level) maxlevel = level;
    }

    (*p)->num_children = 0;
    for (i = 0; i < 8; i++)
        (*p)->child[i] = NULL;


    //printf("Node #%d. level = %d, ibeg = %d, iend = %d.\n", (*p)->node_index, level, ibeg, iend);

    
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
                idx = (*p)->num_children - 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                struct tnode **paddress = &((*p)->child[idx]);

                #pragma omp task firstprivate(paddress,ind,i,lxyzmm,loclev), \
                                 shared(sources,maxparnode)
                pc_create_tree_n0(paddress,
                                  sources, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev);

            }
        }

        #pragma omp taskwait
        
    } else {
   
        #pragma omp critical
        {     
            if (level < minlevel) minlevel = level;
            if (minpars > (*p)->numpar) minpars = (*p)->numpar;
            if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
            /* increment number of leaves */
            numleaves++;
        }
    }

    return;

} /* END of function create_tree_n0 */


int pc_set_tree_index(struct tnode *p, int index)
{
        int current_index = index;
        p->node_index = current_index;

        for (int i = 0; i < p->num_children; i++)
            current_index = pc_set_tree_index(p->child[i], current_index + 1);

        return current_index;
}




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



void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array)
{
    //    printf("Entering pc_create_tree_array.\n");
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
    tree_array->numpar[p->node_index] = p->numpar;
    tree_array->level[p->node_index] = p->level;
    tree_array->radius[p->node_index] = p->radius;
    tree_array->cluster_ind[p->node_index] = p->node_index;

    for (i = 0; i < p->num_children; i++) {
        pc_create_tree_array(p->child[i], tree_array);
    }
    
    return;
    
} /* END of function create_tree_n0 */



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order, int numDevices, int numThreads, struct tnode_array * tree_array){

    int pointsPerCluster = (order+1)*(order+1)*(order+1);
    int numInterpPoints = numnodes * pointsPerCluster;
    make_vector(clusters->x, numInterpPoints);
    make_vector(clusters->y, numInterpPoints);
    make_vector(clusters->z, numInterpPoints);
    make_vector(clusters->q, numInterpPoints);
    make_vector(clusters->w, numInterpPoints);  // will be used in singularity subtraction
    clusters->num=numInterpPoints;

    for (int i = 0; i < numInterpPoints; i++) {
        clusters->q[i]=0.0;
        clusters->w[i]=0.0;
    }

#pragma omp parallel num_threads(numThreads)
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

        for (int i = 0; i < clusters->num; i++) {
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
            #pragma omp for schedule(guided)
            for (int i = 0; i < numnodes; i++) {  // start from i=1, don't need to compute root moments
                pc_comp_ms_modifiedF(tree_array, i, xS, yS, zS, qS, wS,
                                     tempX, tempY, tempZ, tempQ);
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

        free_vector(tempX);
        free_vector(tempY);
        free_vector(tempZ);
        free_vector(tempQ);

    } // end OMP PARALLEL REGION

    return;
}




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




void pc_comp_ms_modifiedF(struct tnode_array * tree_array, int idx,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{
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
#pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ,tt) \
    create(modifiedF[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
            nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],weights[0:torderlim],dj[0:torderlim])
    {

    #pragma acc loop independent
    for (j = 0; j < pointsInNode; j++) {
        modifiedF[j] = qS[startingIndexInSources+j] * wS[startingIndexInSources+j];
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
    for (i = 0; i < pointsInNode; i++) { // loop through the source points

        sumX=0.0;
        sumY=0.0;
        sumZ=0.0;

        sx = xS[startingIndexInSources+i];
        sy = yS[startingIndexInSources+i];
        sz = zS[startingIndexInSources+i];

        #pragma acc loop independent
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

    }

    // Compute moments for each interpolation point
    double numerator, xn, yn, zn, temp;
    int k1, k2, k3, kk;
    double w1,w2,w3;

    #pragma acc loop independent
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
        #pragma acc loop independent
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




void pc_make_interaction_list(const struct tnode_array *tree_array, struct batch *batches,
                              int *tree_inter_list, int *direct_inter_list)
{
    /* local variables */
    int i, j;

    int **batches_ind;
    double **batches_center;
    double *batches_radius;
    
    int tree_numnodes;
    const int *tree_numpar, *tree_level;
    const double *tree_x_mid, *tree_y_mid, *tree_z_mid, *tree_radius;

    batches_ind = batches->index;
    batches_center = batches->center;
    batches_radius = batches->radius;

    tree_numnodes = tree_array->numnodes;
    tree_numpar = tree_array->numpar;
    tree_level = tree_array->level;
    tree_radius = tree_array->radius;
    tree_x_mid = tree_array->x_mid;
    tree_y_mid = tree_array->y_mid;
    tree_z_mid = tree_array->z_mid;

    for (i = 0; i < batches->num * numnodes; i++)
        tree_inter_list[i] = -1;

    for (i = 0; i < batches->num * numleaves; i++)
        direct_inter_list[i] = -1;
    
    for (i = 0; i < batches->num; i++)
        pc_compute_interaction_list(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,
                batches_ind[i], batches_center[i], batches_radius[i],
                &(tree_inter_list[i*numnodes]), &(direct_inter_list[i*numleaves]));

    return;

} /* END of function pc_treecode */



void pc_compute_interaction_list(int tree_numnodes, const int *tree_level, 
                const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list)
{
    /* local variables */
    double tx, ty, tz, dist;
    int j, current_level;
    int tree_index_counter, direct_index_counter;
    
    tree_index_counter = 0;
    direct_index_counter = 0;
    current_level = 0;
    
    for (j = 0; j < tree_numnodes; j++) {
        if (tree_level[j] <= current_level) {
            
            /* determine DIST for MAC test */
            tx = batch_mid[0] - tree_x_mid[j];
            ty = batch_mid[1] - tree_y_mid[j];
            tz = batch_mid[2] - tree_z_mid[j];
            dist = sqrt(tx*tx + ty*ty + tz*tz);

            if (((tree_radius[j] + batch_rad) < dist * sqrt(thetasq))
                && (tree_radius[j] > 0.00)
                && (torder*torder*torder < tree_numpar[j])) {
                current_level = tree_level[j];
            /*
             * If MAC is accepted and there is more than 1 particle
             * in the box, use the expansion for the approximation.
             */
        
                batch_tree_list[tree_index_counter] = j;
                tree_index_counter++;
        
            } else {
            /*
             * If MAC fails check to see if there are children. If not, perform direct
             * calculation. If there are children, call routine recursively for each.
             */

                if (tree_level[j+1] <= tree_level[j]) {
                    
                    batch_direct_list[direct_index_counter] = j;
                    direct_index_counter++;
            
                } else {
                    
                    current_level = tree_level[j+1];
                    
                }
            }
        }
    }

    // Setting tree and direct index counter for batch
    batch_ind[2] = tree_index_counter;
    batch_ind[3] = direct_index_counter;
    
    return;
}


void pc_compute_interaction_list_remote(int tree_numnodes, const int *tree_level,
                                        const int *tree_numpar, const double *tree_radius,
                                        const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                                        int *batch_ind, double *batch_mid, double batch_rad,
                                        int *batch_tree_list, int *batch_direct_list)
{
    /* local variables */
    double tx, ty, tz, dist;
    int j, current_level;
    
    current_level = 0;
    
    for (j = 0; j < tree_numnodes; j++) {
        if (tree_level[j] <= current_level) {
            
            /* determine DIST for MAC test */
            tx = batch_mid[0] - tree_x_mid[j];
            ty = batch_mid[1] - tree_y_mid[j];
            tz = batch_mid[2] - tree_z_mid[j];
            dist = sqrt(tx*tx + ty*ty + tz*tz);
            
            if (((tree_radius[j] + batch_rad) < dist * sqrt(thetasq))
                && (tree_radius[j] > 0.00)
                && (torder*torder*torder < tree_numpar[j])) {
                current_level = tree_level[j];
                /*
                 * If MAC is accepted and there is more than 1 particle
                 * in the box, use the expansion for the approximation.
                 */
                
                batch_tree_list[j] = j;
                
            } else {
                /*
                 * If MAC fails check to see if there are children. If not, perform direct
                 * calculation. If there are children, call routine recursively for each.
                 */
                
                if (tree_level[j+1] <= tree_level[j]) {
                    
                    batch_direct_list[j] = j;
                    
                } else {
                    
                    current_level = tree_level[j+1];
                    
                }
            }
        }
    }
    
    return;
}



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

//          printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//          printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

            int * ibegs = tree_array->ibeg;
            int * iends = tree_array->iend;

            int * clusterInd = tree_array->cluster_ind;

    #pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], \
                            qS[0:sources->num], wS[0:sources->num], \
                            xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
                            xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num], \
                            tree_inter_list[0:numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
                            ibegs[0:numnodes], iends[0:numnodes]) copy(EnP3[0:targets->num], EnP2[0:targets->num])
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

        #pragma omp for private(j,ii,jj,batch_ibeg,batch_iend,numberOfClusterApproximations,\
                                numberOfDirectSums,numberOfTargets,batchStart,node_index,clusterStart,streamID)
        for (i = 0; i < batches->num; i++) {
            batch_ibeg = batches->index[i][0];
            batch_iend = batches->index[i][1];
            numberOfClusterApproximations = batches->index[i][2];
            numberOfDirectSums = batches->index[i][3];

            numberOfTargets = batch_iend - batch_ibeg + 1;
            batchStart =  batch_ibeg - 1;

            for (j = 0; j < numberOfClusterApproximations; j++) {
                node_index = tree_inter_list[i * numnodes + j];
//                clusterStart = numberOfInterpolationPoints*node_index;
                clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                streamID = j%3;
                #pragma acc kernels async(streamID) //present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
                {
                #pragma acc loop independent
                for (ii = 0; ii < numberOfTargets; ii++) {
                    tempPotential = 0.0;
                    xi = xT[ batchStart + ii];
                    yi = yT[ batchStart + ii];
                    zi = zT[ batchStart + ii];

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
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
        } // end loop over target batches

        #pragma acc wait
        #pragma omp barrier
        } // end acc data region
        double totalDueToApprox=0.0; double totalDueToDirect=0.0;
        totalDueToApprox = sum(EnP3,targets->num);
        totalDueToDirect = sum(EnP2,targets->num);
        printf("Potential due to approximations: %f\n",totalDueToApprox);
        printf("Potential due to direct: %f\n",totalDueToDirect);
        for (int k = 0; k < targets->num; k++) {
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

        *tpeng = sum(EnP, targets->num);

        return;

    } /* END of function pc_treecode */
