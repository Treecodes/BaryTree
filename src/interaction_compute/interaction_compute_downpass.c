#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"
#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "interaction_compute.h"


static void cp_comp_pot(struct Tree *tree, int idx, double *potential, int interp_degree,
                        double *xT, double *yT, double *zT, double *clusterQ);

static void cp_comp_pot_parent_to_child(struct Tree *tree, int parent_index, int child_index, int interp_degree,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q);

static void cp_comp_pot_SS(struct Tree *tree, int idx, double *potential, int interp_degree,
                        double *xT, double *yT, double *zT, double *qT,
                        double *clusterQ, double *clusterW);

static void cp_comp_pot_SS_parent_to_child(struct Tree *tree, int parent_index, int child_index, int interp_degree,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q, double *cluster_w);

static void cp_comp_pot_hermite(struct Tree *tree, int idx, double *potential, int interp_degree,
                        double *xT, double *yT, double *zT, double *clusterQ, double *clusterW);

//static void cp_comp_pot_hermite_SS(struct Tree *tree, int idx, int interp_degree,
//                      int totalNumberInterpolationPoints,
//                      double *xT, double *yT, double *zT, double *qT,
//                      double *clusterQ, double *clusterW);


void InteractionCompute_Downpass(double *potential, struct Tree *tree,
                                 struct Particles *targets, struct Clusters *clusters,
                                 struct RunParams *run_params)
{
    int num_targets   = targets->num;
    double *target_x  = targets->x;
    double *target_y  = targets->y;
    double *target_z  = targets->z;
    double *target_q  = targets->q;

    double *cluster_x = clusters->x;
    double *cluster_y = clusters->y;
    double *cluster_z = clusters->z;
    double *cluster_q = clusters->q;
    double *cluster_w = clusters->w;

    int total_num_interp_pts = clusters->num;
    int total_num_interp_charges = clusters->num_charges;
    int total_num_interp_weights = clusters->num_weights;

    int tree_numnodes = tree->numnodes;
    int interp_degree = run_params->interp_degree;


    if ((run_params->approximation == LAGRANGE) && (run_params->singularity == SKIPPING)) {

        // interpolate up clusters, level by level
        for (int level = 0; level < tree->max_depth; ++level) {

            for (int cluster_index = 0; cluster_index < tree->levels_list_num[level]; ++cluster_index) {

                int parent_index = tree->levels_list[level][cluster_index];

                for (int child_counter=0; child_counter < tree->num_children[parent_index]; ++child_counter) {

                    int child_index = tree->children[8*parent_index + child_counter];

                    cp_comp_pot_parent_to_child(tree, parent_index, child_index, interp_degree,
                        cluster_x, cluster_y, cluster_z, cluster_q);
                }
            }
        }

        // interpolate from leaf cluster interpolation points to target particles

        for (int i = 0; i < tree->leaves_list_num; ++i) {
            int leaf_index = tree->leaves_list[i];
            cp_comp_pot(tree, leaf_index, potential, interp_degree,
                        target_x, target_y, target_z, cluster_q);
        }


    } else if ((run_params->approximation == LAGRANGE) && (run_params->singularity == SUBTRACTION)) {

        // interpolate up clusters, level by level
        for (int level = 0; level < tree->max_depth; ++level) {

            for (int cluster_index = 0; cluster_index < tree->levels_list_num[level]; ++cluster_index) {

                int parent_index = tree->levels_list[level][cluster_index];

                for (int child_counter=0; child_counter < tree->num_children[parent_index]; ++child_counter) {

                    int child_index = tree->children[8*parent_index + child_counter];

                    cp_comp_pot_SS_parent_to_child(tree, parent_index, child_index, interp_degree,
                        cluster_x, cluster_y, cluster_z, cluster_q, cluster_w);
                }
            }
        }

        // interpolate from leaf cluster interpolation points to target particles

        for (int i = 0; i < tree->leaves_list_num; ++i) {
            int leaf_index = tree->leaves_list[i];
            cp_comp_pot_SS(tree, leaf_index, potential, interp_degree,
                           target_x, target_y, target_z, target_q, cluster_q, cluster_w);
        }

//        for (int i = 0; i < tree_numnodes; i++){
//            cp_comp_pot_SS(tree, i, potential, interp_degree,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);
//        }

    } else if ((run_params->approximation == HERMITE) && (run_params->singularity == SKIPPING)) {
        for (int i = 0; i < tree_numnodes; i++)
            cp_comp_pot_hermite(tree, i, potential, interp_degree,
                        target_x, target_y, target_z, cluster_q, cluster_w);

    } else if ((run_params->approximation == HERMITE) && (run_params->singularity == SUBTRACTION)) {
        printf("Not set up to do Hermite SS downpass.\n");
        exit(-1);
//        for (int i = 0; i < tree_numnodes; i++)
//            cp_comp_pot_hermite_SS(tree, i, potential, interp_degree,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else {
        exit(1);
    }
        
#ifdef OPENACC_ENABLED
    #pragma acc wait
#endif
    
    return;
}




/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/

void cp_comp_pot_parent_to_child(struct Tree *tree, int parent_index, int child_index, int interp_degree,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q)
{
    int interp_degree_lim       = interp_degree + 1;
    int interp_pts_per_cluster = interp_degree_lim * interp_degree_lim * interp_degree_lim;

//    int num_targets_in_cluster = tree->iend[parent_index] - tree->ibeg[parent_index] + 1;
//    int target_start           = tree->ibeg[idx] - 1;

    int parent_cluster_start = parent_index * interp_pts_per_cluster;
    int child_cluster_start  = child_index  * interp_pts_per_cluster;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ;

    make_vector(weights, interp_degree_lim);
    make_vector(dj,      interp_degree_lim);
    make_vector(tt,      interp_degree_lim);
    make_vector(nodeX,   interp_degree_lim);
    make_vector(nodeY,   interp_degree_lim);
    make_vector(nodeZ,   interp_degree_lim);

    double x0 = tree->x_min[parent_index];
    double x1 = tree->x_max[parent_index];
    double y0 = tree->y_min[parent_index];
    double y1 = tree->y_max[parent_index];
    double z0 = tree->z_min[parent_index];
    double z1 = tree->z_max[parent_index];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(cluster_x, cluster_y, cluster_z, cluster_q) \
                create(nodeX[0:interp_degree_lim], nodeY[0:interp_degree_lim], nodeZ[0:interp_degree_lim], \
                       weights[0:interp_degree_lim], dj[0:interp_degree_lim], tt[0:interp_degree_lim])
    {
#endif


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_degree_lim; i++) {
        tt[i] = cos(i * M_PI / interp_degree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++){
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interp_degree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_pts_per_cluster; i++) { // loop through the child cluster points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = cluster_x[child_cluster_start+i];
        double ty = cluster_y[child_cluster_start+i];
        double tz = cluster_z[child_cluster_start+i];

        int eix = -1;
        int eiy = -1;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
        for (int j = 0; j < interp_degree_lim; j++) {  // loop through the degree

            double cx = tx - nodeX[j];
            double cy = ty - nodeY[j];
            double cz = tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) eix = j;
            if (fabs(cy)<DBL_MIN) eiy = j;
            if (fabs(cz)<DBL_MIN) eiz = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (eix==-1) denominator /= sumX;
        if (eiy==-1) denominator /= sumY;
        if (eiz==-1) denominator /= sumZ;

        double temp = 0.0;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over parent interpolation points

            int k1 = j%interp_degree_lim;
            int kk = (j-k1)/interp_degree_lim;
            int k2 = kk%interp_degree_lim;
            kk = kk - k2;
            int k3 = kk / interp_degree_lim;

            double w3 = weights[k3];
            double w2 = weights[k2];
            double w1 = weights[k1];

            double cx = nodeX[k1];
            double cy = nodeY[k2];
            double cz = nodeZ[k3];
            double cq = cluster_q[parent_cluster_start + j];

            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (eix == -1) {
                numerator *= w1 / (tx - cx);
            } else {
                if (eix != k1) numerator *= 0;
            }

            if (eiy == -1) {
                numerator *= w2 / (ty - cy);
            } else {
                if (eiy != k2) numerator *= 0;
            }

            if (eiz == -1) {
                numerator *= w3 / (tz - cz);
            } else {
                if (eiz != k3) numerator *= 0;
            }

            temp += numerator * denominator * cq;
        }

#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_q[i + child_cluster_start] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif

    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}


void cp_comp_pot(struct Tree *tree, int idx, double *potential, int interp_degree,
        double *target_x, double *target_y, double *target_z,
        double *cluster_q)
{
    int interp_degree_lim       = interp_degree + 1;
    int interp_pts_per_cluster = interp_degree_lim * interp_degree_lim * interp_degree_lim;

    int num_targets_in_cluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int target_start           = tree->ibeg[idx] - 1;
    int cluster_start          = idx * interp_pts_per_cluster;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ;

    make_vector(weights, interp_degree_lim);
    make_vector(dj,      interp_degree_lim);
    make_vector(tt,      interp_degree_lim);
    make_vector(nodeX,   interp_degree_lim);
    make_vector(nodeY,   interp_degree_lim);
    make_vector(nodeZ,   interp_degree_lim);

    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(target_x, target_y, target_z, cluster_q) \
                create(nodeX[0:interp_degree_lim], nodeY[0:interp_degree_lim], nodeZ[0:interp_degree_lim], \
                       weights[0:interp_degree_lim], dj[0:interp_degree_lim], tt[0:interp_degree_lim])
    {
#endif


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_degree_lim; i++) {
        tt[i] = cos(i * M_PI / interp_degree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++){
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interp_degree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < num_targets_in_cluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = target_x[target_start+i];
        double ty = target_y[target_start+i];
        double tz = target_z[target_start+i];

        int eix = -1;
        int eiy = -1;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
        for (int j = 0; j < interp_degree_lim; j++) {  // loop through the degree

            double cx = tx - nodeX[j];
            double cy = ty - nodeY[j];
            double cz = tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) eix = j;
            if (fabs(cy)<DBL_MIN) eiy = j;
            if (fabs(cz)<DBL_MIN) eiz = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (eix==-1) denominator /= sumX;
        if (eiy==-1) denominator /= sumY;
        if (eiz==-1) denominator /= sumZ;

        double temp = 0.0;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interp_degree_lim;
            int kk = (j-k1)/interp_degree_lim;
            int k2 = kk%interp_degree_lim;
            kk = kk - k2;
            int k3 = kk / interp_degree_lim;

            double w3 = weights[k3];
            double w2 = weights[k2];
            double w1 = weights[k1];

            double cx = nodeX[k1];
            double cy = nodeY[k2];
            double cz = nodeZ[k3];
            double cq = cluster_q[cluster_start + j];

            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (eix == -1) {
                numerator *= w1 / (tx - cx);
            } else {
                if (eix != k1) numerator *= 0;
            }

            if (eiy == -1) {
                numerator *= w2 / (ty - cy);
            } else {
                if (eiy != k2) numerator *= 0;
            }

            if (eiz == -1) {
                numerator *= w3 / (tz - cz);
            } else {
                if (eiz != k3) numerator *= 0;
            }

            temp += numerator * denominator * cq;
        }

#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[i + target_start] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif

    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}




void cp_comp_pot_SS_parent_to_child(struct Tree *tree, int parent_index, int child_index, int interp_degree,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q, double *cluster_w)
{
    int interp_degree_lim       = interp_degree + 1;
    int interp_pts_per_cluster = interp_degree_lim * interp_degree_lim * interp_degree_lim;

    int parent_cluster_start = parent_index * interp_pts_per_cluster;
    int child_cluster_start  = child_index  * interp_pts_per_cluster;
    
    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ;

    make_vector(weights, interp_degree_lim);
    make_vector(dj,      interp_degree_lim);
    make_vector(tt,      interp_degree_lim);
    make_vector(nodeX,   interp_degree_lim);
    make_vector(nodeY,   interp_degree_lim);
    make_vector(nodeZ,   interp_degree_lim);
    
    double x0 = tree->x_min[parent_index];
    double x1 = tree->x_max[parent_index];
    double y0 = tree->y_min[parent_index];
    double y1 = tree->y_max[parent_index];
    double z0 = tree->z_min[parent_index];
    double z1 = tree->z_max[parent_index];
    
#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(cluster_x, cluster_y, cluster_z, cluster_q, cluster_w) \
                create(nodeX[0:interp_degree_lim], nodeY[0:interp_degree_lim], nodeZ[0:interp_degree_lim], \
                       weights[0:interp_degree_lim], dj[0:interp_degree_lim], tt[0:interp_degree_lim])
    {
#endif
    

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_degree_lim; i++) {
        tt[i] = cos(i * M_PI / interp_degree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    
    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++){
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interp_degree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_pts_per_cluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = cluster_x[child_cluster_start+i];
        double ty = cluster_y[child_cluster_start+i];
        double tz = cluster_z[child_cluster_start+i];

        int eix = -1;
        int eiy = -1;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
        for (int j = 0; j < interp_degree_lim; j++) {  // loop through the degree

            double cx = tx - nodeX[j];
            double cy = ty - nodeY[j];
            double cz = tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) eix = j;
            if (fabs(cy)<DBL_MIN) eiy = j;
            if (fabs(cz)<DBL_MIN) eiz = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (eix==-1) denominator /= sumX;
        if (eiy==-1) denominator /= sumY;
        if (eiz==-1) denominator /= sumZ;
        
        double temp = 0.0;
        double temp2 = 0.0;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp) reduction(+:temp2)
#endif
        for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interp_degree_lim;
            int kk = (j-k1)/interp_degree_lim;
            int k2 = kk%interp_degree_lim;
            kk = kk - k2;
            int k3 = kk / interp_degree_lim;

            double w3 = weights[k3];
            double w2 = weights[k2];
            double w1 = weights[k1];
            
            double cx = nodeX[k1];
            double cy = nodeY[k2];
            double cz = nodeZ[k3];
            double cq = cluster_q[parent_cluster_start + j];
            double cw = cluster_w[parent_cluster_start + j];
        
            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (eix == -1) {
                numerator *= w1 / (tx - cx);
            } else {
                if (eix != k1) numerator *= 0;
            }

            if (eiy == -1) {
                numerator *= w2 / (ty - cy);
            } else {
                if (eiy != k2) numerator *= 0;
            }

            if (eiz == -1) {
                numerator *= w3 / (tz - cz);
            } else {
                if (eiz != k3) numerator *= 0;
            }

            temp  += numerator * denominator * cq;  
            temp2 += numerator * denominator * cw;
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_q[child_cluster_start + i] += temp;
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_w[child_cluster_start + i] += temp2;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}




void cp_comp_pot_SS(struct Tree *tree, int idx, double *potential, int interp_degree,
        double *target_x, double *target_y, double *target_z, double *target_q,
        double *cluster_q, double *cluster_w)
{
    int interp_degree_lim       = interp_degree + 1;
    int interp_pts_per_cluster = interp_degree_lim * interp_degree_lim * interp_degree_lim;

    int num_targets_in_cluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int target_start           = tree->ibeg[idx] - 1;
    int cluster_start          = idx * interp_pts_per_cluster;
    
    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ;

    make_vector(weights, interp_degree_lim);
    make_vector(dj,      interp_degree_lim);
    make_vector(tt,      interp_degree_lim);
    make_vector(nodeX,   interp_degree_lim);
    make_vector(nodeY,   interp_degree_lim);
    make_vector(nodeZ,   interp_degree_lim);
    
    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];
    
#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(target_x, target_y, target_z, target_q, cluster_q, cluster_w, potential) \
                create(nodeX[0:interp_degree_lim], nodeY[0:interp_degree_lim], nodeZ[0:interp_degree_lim], \
                       weights[0:interp_degree_lim], dj[0:interp_degree_lim], tt[0:interp_degree_lim])
    {
#endif
    

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_degree_lim; i++) {
        tt[i] = cos(i * M_PI / interp_degree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    
    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++){
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interp_degree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < num_targets_in_cluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = target_x[target_start+i];
        double ty = target_y[target_start+i];
        double tz = target_z[target_start+i];
        double tq = target_q[target_start+i];

        int eix = -1;
        int eiy = -1;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
        for (int j = 0; j < interp_degree_lim; j++) {  // loop through the degree

            double cx = tx - nodeX[j];
            double cy = ty - nodeY[j];
            double cz = tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) eix = j;
            if (fabs(cy)<DBL_MIN) eiy = j;
            if (fabs(cz)<DBL_MIN) eiz = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (eix==-1) denominator /= sumX;
        if (eiy==-1) denominator /= sumY;
        if (eiz==-1) denominator /= sumZ;
        
        double temp = 0.0;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interp_degree_lim;
            int kk = (j-k1)/interp_degree_lim;
            int k2 = kk%interp_degree_lim;
            kk = kk - k2;
            int k3 = kk / interp_degree_lim;

            double w3 = weights[k3];
            double w2 = weights[k2];
            double w1 = weights[k1];
            
            double cx = nodeX[k1];
            double cy = nodeY[k2];
            double cz = nodeZ[k3];
            double cq = cluster_q[cluster_start + j];
            double cw = cluster_w[cluster_start + j];
        
            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (eix == -1) {
                numerator *= w1 / (tx - cx);
            } else {
                if (eix != k1) numerator *= 0;
            }

            if (eiy == -1) {
                numerator *= w2 / (ty - cy);
            } else {
                if (eiy != k2) numerator *= 0;
            }

            if (eiz == -1) {
                numerator *= w3 / (tz - cz);
            } else {
                if (eiz != k3) numerator *= 0;
            }

            temp += numerator * denominator * (cq-tq*cw);  // subtract target_q*cluster_w for singularity subtraction
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[i + target_start] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}



void cp_comp_pot_hermite(struct Tree *tree, int idx, double *potential, int interp_degree,
        double *target_x, double *target_y, double *target_z, double *cluster_q, double *cluster_w)
{
    int interp_degree_lim       = interp_degree + 1;
    int interp_pts_per_cluster = interp_degree_lim * interp_degree_lim * interp_degree_lim;

    int num_targets_in_cluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int target_start           = tree->ibeg[idx] - 1;
    int cluster_start          = idx * interp_pts_per_cluster;
    
    double *dj, *tt, *ww, *wx, *wy, *wz, *nodeX, *nodeY, *nodeZ;

    make_vector(dj,      interp_degree_lim);
    make_vector(tt,      interp_degree_lim);
    make_vector(ww,      interp_degree_lim);
    make_vector(wx,      interp_degree_lim);
    make_vector(wy,      interp_degree_lim);
    make_vector(wz,      interp_degree_lim);
    make_vector(nodeX,   interp_degree_lim);
    make_vector(nodeY,   interp_degree_lim);
    make_vector(nodeZ,   interp_degree_lim);
    
    double *cluster_q_     = &cluster_q[8*cluster_start + 0*interp_pts_per_cluster];
    double *cluster_q_dx   = &cluster_q[8*cluster_start + 1*interp_pts_per_cluster];
    double *cluster_q_dy   = &cluster_q[8*cluster_start + 2*interp_pts_per_cluster];
    double *cluster_q_dz   = &cluster_q[8*cluster_start + 3*interp_pts_per_cluster];
    double *cluster_q_dxy  = &cluster_q[8*cluster_start + 4*interp_pts_per_cluster];
    double *cluster_q_dyz  = &cluster_q[8*cluster_start + 5*interp_pts_per_cluster];
    double *cluster_q_dxz  = &cluster_q[8*cluster_start + 6*interp_pts_per_cluster];
    double *cluster_q_dxyz = &cluster_q[8*cluster_start + 7*interp_pts_per_cluster];
    
    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];
    
#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(target_x, target_y, target_z, \
                                        cluster_q_, cluster_q_dx, cluster_q_dy, cluster_q_dz, \
                                        cluster_q_dxy, cluster_q_dyz, cluster_q_dxz, \
                                        cluster_q_dxyz) \
        create(nodeX[0:interp_degree_lim], nodeY[0:interp_degree_lim], nodeZ[0:interp_degree_lim], \
               dj[0:interp_degree_lim], tt[0:interp_degree_lim], ww[0:interp_degree_lim], \
               wx[0:interp_degree_lim], wy[0:interp_degree_lim], wz[0:interp_degree_lim])
    {
#endif
    
    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_degree_lim; i++) {
        double xx = i * M_PI / interp_degree;
        tt[i] =  cos(xx);
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    ww[0] = 0.25 * (interp_degree*interp_degree/3.0 + 1.0/6.0);
    ww[interp_degree] = -ww[0];
    
    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_degree_lim; j++){
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (x1 - x0);
        wy[j] = -4.0 * ww[j] / (y1 - y0);
        wz[j] = -4.0 * ww[j] / (z1 - z0);
    }
    dj[0] = 0.25;
    dj[interp_degree] = 0.25;

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < num_targets_in_cluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = target_x[target_start+i];
        double ty = target_y[target_start+i];
        double tz = target_z[target_start+i];

        int eix = -1;
        int eiy = -1;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
        for (int j = 0; j < interp_degree_lim; j++) {  // loop through the degree

            double cx =  tx - nodeX[j];
            double cy =  ty - nodeY[j];
            double cz =  tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) eix = j;
            if (fabs(cy)<DBL_MIN) eiy = j;
            if (fabs(cz)<DBL_MIN) eiz = j;

            // Increment the sums
            sumX += dj[j] / (cx*cx) + wx[j] / cx;
            sumY += dj[j] / (cy*cy) + wy[j] / cy;
            sumZ += dj[j] / (cz*cz) + wz[j] / cz;

        }

        double denominator = 1.0;
        if (eix==-1) denominator /= sumX;
        if (eiy==-1) denominator /= sumY;
        if (eiz==-1) denominator /= sumZ;
        
        double temp = 0.0;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interp_degree_lim;
            int kk = (j-k1)/interp_degree_lim;
            int k2 = kk%interp_degree_lim;
            kk = kk - k2;
            int k3 = kk / interp_degree_lim;
            
            double dx = tx - nodeX[k1];
            double dy = ty - nodeY[k2];
            double dz = tz - nodeZ[k3];
            
            double cq     = cluster_q_[j];
            double cqdx   = cluster_q_dx[j];
            double cqdy   = cluster_q_dy[j];
            double cqdz   = cluster_q_dz[j];
            double cqdxy  = cluster_q_dxy[j];
            double cqdyz  = cluster_q_dyz[j];
            double cqdxz  = cluster_q_dxz[j];
            double cqdxyz = cluster_q_dxyz[j];
        
            double numerator0 = 1.0, numerator1 = 1.0, numerator2 = 1.0, numerator3 = 1.0;
            double numerator4 = 1.0, numerator5 = 1.0, numerator6 = 1.0, numerator7 = 1.0;

            double Ax = dj[k1] / (dx*dx) + wx[k1] / dx;
            double Ay = dj[k2] / (dy*dy) + wy[k2] / dy;
            double Az = dj[k3] / (dz*dz) + wz[k3] / dz;
            double Bx = dj[k1] / dx;
            double By = dj[k2] / dy;
            double Bz = dj[k3] / dz;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (eix == -1) {
                numerator0 *=  Ax;                     // Aaa

                numerator1 *=  Bx;                     // Baa
                numerator2 *=  Ax;                     // Aba
                numerator3 *=  Ax;                     // Aab

                numerator4 *=  Bx;                     // Bba
                numerator5 *=  Ax;                     // Abb
                numerator6 *=  Bx;                     // Bab

                numerator7 *=  Bx;                     // Bbb
                
            } else {
                if (eix != k1) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                } else {
                    numerator1 *= 0; numerator4 *= 0; numerator6 *= 0; numerator7 *= 0;
                }
            }

            if (eiy == -1) {
                numerator0 *=  Ay;                     // aAa

                numerator1 *=  Ay;                     // bAa
                numerator2 *=  By;                     // aBa
                numerator3 *=  Ay;                     // aAb

                numerator4 *=  By;                     // bBa
                numerator5 *=  By;                     // aBb
                numerator6 *=  Ay;                     // bAb

                numerator7 *=  By;                     // bBb
                
            } else {
                if (eiy != k2) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                }  else {
                    numerator2 *= 0; numerator4 *= 0; numerator5 *= 0; numerator7 *= 0;
                }
            }

            if (eiz == -1) {
                numerator0 *=  Az;                    // aaA

                numerator1 *=  Az;                    // baA
                numerator2 *=  Az;                    // abA
                numerator3 *=  Bz;                    // aaB

                numerator4 *=  Az;                    // bbA
                numerator5 *=  Bz;                    // abB
                numerator6 *=  Bz;                    // baB
                
                numerator7 *=  Bz;                    // bbB
                
            } else {
                if (eiz != k3) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                } else {
                    numerator3 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                }
            }

            temp += denominator * (numerator0 * cq     +  numerator1 * cqdx   +  numerator2 * cqdy
                                +  numerator3 * cqdz   +  numerator4 * cqdxy  +  numerator5 * cqdyz
                                +  numerator6 * cqdxz  +  numerator7 * cqdxyz);
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[i + target_start] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(dj);
    free_vector(tt);
    free_vector(ww);
    free_vector(wx);
    free_vector(wy);
    free_vector(wz);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}
