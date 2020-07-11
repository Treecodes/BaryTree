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

static void cp_comp_pot(struct Tree *tree, int idx, double *potential, int interp_order,
        
        int target_x_low_ind,  int target_x_high_ind,
        int target_y_low_ind,  int target_y_high_ind,
        int target_z_low_ind,  int target_z_high_ind,
        double target_xmin,    double target_ymin,    double target_zmin,
        double target_xmax,    double target_ymax,    double target_zmax,
                    
        double target_xdd,     double target_ydd,     double target_zdd,
        int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,
        
        double *cluster_q);
        
static void cp_comp_pot_hermite(struct Tree *tree, int idx, double *potential, int interp_order,
                
        int target_x_low_ind,  int target_x_high_ind,
        int target_y_low_ind,  int target_y_high_ind,
        int target_z_low_ind,  int target_z_high_ind,
        double target_xmin,    double target_ymin,    double target_zmin,
        double target_xmax,    double target_ymax,    double target_zmax,
                    
        double target_xdd,     double target_ydd,     double target_zdd,
        int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,
        
        double *cluster_q);


void InteractionCompute_Downpass(double *potential, struct Tree *tree,
                                 struct Particles *targets, struct Clusters *clusters,
                                 struct RunParams *run_params)
{
    int total_num_interp_charges = clusters->num_charges;
    double *cluster_q = clusters->q;

    int tree_numnodes = tree->numnodes;
    int interp_order = run_params->interp_order;

    int num_targets = targets->num;
    
    int *target_tree_x_low_ind = tree->x_low_ind;
    int *target_tree_y_low_ind = tree->y_low_ind;
    int *target_tree_z_low_ind = tree->z_low_ind;
    
    int *target_tree_x_high_ind = tree->x_high_ind;
    int *target_tree_y_high_ind = tree->y_high_ind;
    int *target_tree_z_high_ind = tree->z_high_ind;
    
    double *target_tree_x_min = tree->x_min;
    double *target_tree_y_min = tree->y_min;
    double *target_tree_z_min = tree->z_min;
    
    double *target_tree_x_max = tree->x_max;
    double *target_tree_y_max = tree->y_max;
    double *target_tree_z_max = tree->z_max;
    
    
    int target_x_dim_glob = targets->xdim;
    int target_y_dim_glob = targets->ydim;
    int target_z_dim_glob = targets->zdim;
    
    double target_xdd = targets->xdd;
    double target_ydd = targets->ydd;
    double target_zdd = targets->zdd;


#ifdef OPENACC_ENABLED
    //#pragma acc data copyin(cluster_q[0:total_num_interp_charges]) \
    //                   copy(potential[0:num_targets])
    //#pragma acc data copy(potential[0:num_targets])
    {
#endif

    if ((run_params->approximation == LAGRANGE) && (run_params->singularity == SKIPPING)) {
        for (int i = 0; i < tree_numnodes; i++) {
            if (tree->used[i] == 1) {
                
                int target_x_low_ind = target_tree_x_low_ind[i];
                int target_y_low_ind = target_tree_y_low_ind[i];
                int target_z_low_ind = target_tree_z_low_ind[i];
    
                int target_x_high_ind = target_tree_x_high_ind[i];
                int target_y_high_ind = target_tree_y_high_ind[i];
                int target_z_high_ind = target_tree_z_high_ind[i];
    
                double target_x_min = target_tree_x_min[i];
                double target_y_min = target_tree_y_min[i];
                double target_z_min = target_tree_z_min[i];
                
                double target_x_max = target_tree_x_max[i];
                double target_y_max = target_tree_y_max[i];
                double target_z_max = target_tree_z_max[i];
        
                cp_comp_pot(tree, i, potential, interp_order,
                
                            target_x_low_ind, target_x_high_ind,
                            target_y_low_ind, target_y_high_ind,
                            target_z_low_ind, target_z_high_ind,
                            target_x_min,       target_y_min,       target_z_min,
                            target_x_max,       target_y_max,       target_z_max,
                                          
                            target_xdd,        target_ydd,        target_zdd,
                            target_x_dim_glob, target_y_dim_glob, target_z_dim_glob,
                            
                            cluster_q);
            }
        }

    } else if ((run_params->approximation == LAGRANGE) && (run_params->singularity == SUBTRACTION)) {
//        for (int i = 0; i < tree_numnodes; i++)
//            cp_comp_pot_SS(tree, i, potential interp_order,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else if ((run_params->approximation == HERMITE) && (run_params->singularity == SKIPPING)) {
        for (int i = 0; i < tree_numnodes; i++) {
            if (tree->used[i] == 1) {
            
                int target_x_low_ind = target_tree_x_low_ind[i];
                int target_y_low_ind = target_tree_y_low_ind[i];
                int target_z_low_ind = target_tree_z_low_ind[i];
    
                int target_x_high_ind = target_tree_x_high_ind[i];
                int target_y_high_ind = target_tree_y_high_ind[i];
                int target_z_high_ind = target_tree_z_high_ind[i];
    
                double target_x_min = target_tree_x_min[i];
                double target_y_min = target_tree_y_min[i];
                double target_z_min = target_tree_z_min[i];
                
                double target_x_max = target_tree_x_max[i];
                double target_y_max = target_tree_y_max[i];
                double target_z_max = target_tree_z_max[i];
                        
                cp_comp_pot_hermite(tree, i, potential, interp_order,
                                                        
                            target_x_low_ind, target_x_high_ind,
                            target_y_low_ind, target_y_high_ind,
                            target_z_low_ind, target_z_high_ind,
                            target_x_min,       target_y_min,       target_z_min,
                            target_x_max,       target_y_max,       target_z_max,
                                          
                            target_xdd,        target_ydd,        target_zdd,
                            target_x_dim_glob, target_y_dim_glob, target_z_dim_glob,
                            
                            cluster_q);
            }
        }

    } else if ((run_params->approximation == HERMITE) && (run_params->singularity == SUBTRACTION)) {
//        for (int i = 0; i < tree_numnodes; i++)
//            cp_comp_pot_hermite_SS(tree, i, potential, interp_order,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else {
        exit(1);
    }
        
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif
    
    return;
}




/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/

void cp_comp_pot(struct Tree *tree, int idx, double *potential, int interp_order,
        
        int target_x_low_ind,  int target_x_high_ind,
        int target_y_low_ind,  int target_y_high_ind,
        int target_z_low_ind,  int target_z_high_ind,
        double target_xmin,    double target_ymin,    double target_zmin,
        double target_xmax,    double target_ymax,    double target_zmax,
                    
        double target_xdd,     double target_ydd,     double target_zdd,
        int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,
        
        double *cluster_q)
{
    int interp_order_lim       = interp_order + 1;
    int interp_pts_per_cluster = interp_order_lim * interp_order_lim * interp_order_lim;
    
    int target_yz_dim = target_y_dim_glob * target_z_dim_glob;

    //int num_targets_in_cluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    //int target_start           = tree->ibeg[idx] - 1;
    int cluster_charge_start          = idx * interp_pts_per_cluster;

    int coeff_x_dim = (target_x_high_ind-target_x_low_ind + 1) * interp_order_lim;
    int coeff_y_dim = (target_y_high_ind-target_y_low_ind + 1) * interp_order_lim;
    int coeff_z_dim = (target_z_high_ind-target_z_low_ind + 1) * interp_order_lim;
    
    double *weights, *nodeX, *nodeY, *nodeZ;
    double *coeffX, *coeffY, *coeffZ;

    make_vector(weights, interp_order_lim);
    make_vector(nodeX,   interp_order_lim);
    make_vector(nodeY,   interp_order_lim);
    make_vector(nodeZ,   interp_order_lim);

    make_vector(coeffX,  coeff_x_dim);
    make_vector(coeffY,  coeff_y_dim);
    make_vector(coeffZ,  coeff_z_dim);
    
#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(potential, cluster_q) \
                create(nodeX[0:interp_order_lim], nodeY[0:interp_order_lim], nodeZ[0:interp_order_lim], \
                       weights[0:interp_order_lim], \
                       coeffX[0:coeff_x_dim], coeffY[0:coeff_y_dim], coeffZ[0:coeff_z_dim])
    {
#endif
    

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_order_lim; i++) {
        double tt = cos(i * M_PI / interp_order);
        nodeX[i] = target_xmin + (tt + 1.0)/2.0 * (target_xmax - target_xmin);
        nodeY[i] = target_ymin + (tt + 1.0)/2.0 * (target_ymax - target_ymin);
        nodeZ[i] = target_zmin + (tt + 1.0)/2.0 * (target_zmax - target_zmin);
        weights[i] = ((i % 2 == 0)? 1 : -1);
        if (i == 0 || i == interp_order) weights[i] = ((i % 2 == 0)? 1 : -1) * 0.5;
    }

    


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int ix = target_x_low_ind; ix <= target_x_high_ind; ix++) {
        double tx = target_xmin + (ix - target_x_low_ind) * target_xdd;
        double denominator = 0.0;
        int eix = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent reduction(+:denominator) reduction(max:eix)
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double cx = tx - nodeX[j];
            if (fabs(cx)<DBL_MIN) eix = j;
            denominator += weights[j] / cx;
        }

        if (eix!=-1) denominator = 1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double numerator = 1.0;
            if (eix == -1) {
                numerator *= weights[j] / (tx - nodeX[j]);
            } else {
                if (eix != j) numerator *= 0;
            }

            coeffX[(ix-target_x_low_ind) * interp_order_lim + j] = numerator / denominator;
        }
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int iy = target_y_low_ind; iy <= target_y_high_ind; iy++) {
        double ty = target_ymin + (iy - target_y_low_ind) * target_ydd;
        double denominator = 0.0;
        int eiy = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent reduction(+:denominator) reduction(max:eiy)
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double cy = ty - nodeY[j];
            if (fabs(cy)<DBL_MIN) eiy = j;
            denominator += weights[j] / cy;
        }

        if (eiy!=-1) denominator = 1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double numerator = 1.0;
            if (eiy == -1) {
                numerator *= weights[j] / (ty - nodeY[j]);
            } else {
                if (eiy != j) numerator *= 0;
            }

            coeffY[(iy-target_y_low_ind) * interp_order_lim + j] = numerator / denominator;
        }
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int iz = target_z_low_ind; iz <= target_z_high_ind; iz++) {
        double tz = target_zmin + (iz - target_z_low_ind) * target_zdd;
        double denominator = 0.0;
        int eiz = -1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent reduction(+:denominator) reduction(max:eiz)
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double cz = tz - nodeZ[j];
            if (fabs(cz)<DBL_MIN) eiz = j;
            denominator += weights[j] / cz;
        }

        if (eiz!=-1) denominator = 1;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) independent
#endif
        for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree
            double numerator = 1.0;
            if (eiz == -1) {
                numerator *= weights[j] / (tz - nodeZ[j]);
            } else {
                if (eiz != j) numerator *= 0;
            }

            coeffZ[(iz-target_z_low_ind) * interp_order_lim + j] = numerator / denominator;
        }
    }



#ifdef OPENACC_ENABLED
    #pragma acc loop gang collapse(3) independent
#endif
    for (int ix = target_x_low_ind; ix <= target_x_high_ind; ix++) {
        for (int iy = target_y_low_ind; iy <= target_y_high_ind; iy++) {
            for (int iz = target_z_low_ind; iz <= target_z_high_ind; iz++) {

                int ii = (ix * target_yz_dim) + (iy * target_z_dim_glob) + iz;
                int iix = (ix - target_x_low_ind) * interp_order_lim;
                int iiy = (iy - target_y_low_ind) * interp_order_lim;
                int iiz = (iz - target_z_low_ind) * interp_order_lim;
                
                double temp = 0.0;

        #ifdef OPENACC_ENABLED
                #pragma acc loop vector(32) independent reduction(+:temp)
        #endif
                for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
                    int k3 = j%interp_order_lim;
                    int kk = (j-k3)/interp_order_lim;
                    int k2 = kk%interp_order_lim;
                    kk = kk - k2;
                    int k1 = kk / interp_order_lim;

                    double cq = cluster_q[cluster_charge_start + j];
                    temp += coeffX[iix + k1] * coeffY[iiy + k2] * coeffZ[iiz + k3] * cq;
   
                }

#ifdef OPENACC_ENABLED
                #pragma acc atomic
#endif
                potential[ii] += temp;
            }
        }
    }


/*
#ifdef OPENACC_ENABLED
    #pragma acc loop collapse(3) independent
#endif
    for (int ix = target_x_low_ind; ix <= target_x_high_ind; ix++) {
        for (int iy = target_y_low_ind; iy <= target_y_high_ind; iy++) {
            for (int iz = target_z_low_ind; iz <= target_z_high_ind; iz++) {

                int ii = (ix * target_yz_dim) + (iy * target_z_dim_glob) + iz;

                double tx = target_xmin + (ix - target_x_low_ind) * target_xdd;
                double ty = target_ymin + (iy - target_y_low_ind) * target_ydd;
                double tz = target_zmin + (iz - target_z_low_ind) * target_zdd;

                double sumX = 0.0;
                double sumY = 0.0;
                double sumZ = 0.0;

                int eix = -1;
                int eiy = -1;
                int eiz = -1;

#ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
                for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree

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
        
//#ifdef OPENACC_ENABLED
//                #pragma acc loop collapse(3) independent reduction(+:temp)
//#endif
//                for (int k1 = 0; k1 < interp_order_lim; k1++) {
//                for (int k2 = 0; k2 < interp_order_lim; k2++) {
//                for (int k3 = 0; k3 < interp_order_lim; k3++) {
//                    double w3 = weights[k3];
//                    double w2 = weights[k2];
//                    double w1 = weights[k1];
//            
//                    double cx = nodeX[k1];
//                    double cy = nodeY[k2];
//                    double cz = nodeZ[k3];
//
//                    int jj = cluster_charge_start + k1 * interp_order_lim*interp_order_lim + k2 * interp_order_lim + k3;
//                    double cq = cluster_q[jj];
        #ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:temp)
        #endif
                for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
        
                    int k3 = j%interp_order_lim;
                    int kk = (j-k3)/interp_order_lim;
                    int k2 = kk%interp_order_lim;
                    kk = kk - k2;
                    int k1 = kk / interp_order_lim;

                    double w3 = weights[k3];
                    double w2 = weights[k2];
                    double w1 = weights[k1];
            
                    double cx = nodeX[k1];
                    double cy = nodeY[k2];
                    double cz = nodeZ[k3];

                    double cq = cluster_q[cluster_charge_start + j];
        
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
                //}
                //}
        
#ifdef OPENACC_ENABLED
                #pragma acc atomic
#endif
                potential[ii] += temp;
            }
        }
    }
*/

#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(weights);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);
    free_vector(coeffX);
    free_vector(coeffY);
    free_vector(coeffZ);

    return;
}



void cp_comp_pot_hermite(struct Tree *tree, int idx, double *potential, int interp_order,
                
        int target_x_low_ind,  int target_x_high_ind,
        int target_y_low_ind,  int target_y_high_ind,
        int target_z_low_ind,  int target_z_high_ind,
        double target_xmin,    double target_ymin,    double target_zmin,
        double target_xmax,    double target_ymax,    double target_zmax,
                    
        double target_xdd,     double target_ydd,     double target_zdd,
        int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,
        
        double *cluster_q)
{
    int interp_order_lim       = interp_order + 1;
    int interp_pts_per_cluster = interp_order_lim * interp_order_lim * interp_order_lim;
    int target_yz_dim = target_y_dim_glob * target_z_dim_glob;

    int num_targets_in_cluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int target_start           = tree->ibeg[idx] - 1;
    int cluster_start          = idx * interp_pts_per_cluster;
    
    double *dj, *tt, *ww, *wx, *wy, *wz, *nodeX, *nodeY, *nodeZ;

    make_vector(dj,      interp_order_lim);
    make_vector(tt,      interp_order_lim);
    make_vector(ww,      interp_order_lim);
    make_vector(wx,      interp_order_lim);
    make_vector(wy,      interp_order_lim);
    make_vector(wz,      interp_order_lim);
    make_vector(nodeX,   interp_order_lim);
    make_vector(nodeY,   interp_order_lim);
    make_vector(nodeZ,   interp_order_lim);
    
    double *cluster_q_     = &cluster_q[8*cluster_start + 0*interp_pts_per_cluster];
    double *cluster_q_dx   = &cluster_q[8*cluster_start + 1*interp_pts_per_cluster];
    double *cluster_q_dy   = &cluster_q[8*cluster_start + 2*interp_pts_per_cluster];
    double *cluster_q_dz   = &cluster_q[8*cluster_start + 3*interp_pts_per_cluster];
    double *cluster_q_dxy  = &cluster_q[8*cluster_start + 4*interp_pts_per_cluster];
    double *cluster_q_dyz  = &cluster_q[8*cluster_start + 5*interp_pts_per_cluster];
    double *cluster_q_dxz  = &cluster_q[8*cluster_start + 6*interp_pts_per_cluster];
    double *cluster_q_dxyz = &cluster_q[8*cluster_start + 7*interp_pts_per_cluster];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(potential, \
                                        cluster_q_, cluster_q_dx, cluster_q_dy, cluster_q_dz, \
                                        cluster_q_dxy, cluster_q_dyz, cluster_q_dxz, \
                                        cluster_q_dxyz) \
        create(nodeX[0:interp_order_lim], nodeY[0:interp_order_lim], nodeZ[0:interp_order_lim], \
               dj[0:interp_order_lim], tt[0:interp_order_lim], ww[0:interp_order_lim], \
               wx[0:interp_order_lim], wy[0:interp_order_lim], wz[0:interp_order_lim])
    {
#endif
    
    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interp_order_lim; i++) {
        double xx = i * M_PI / interp_order;
        tt[i] =  cos(xx);
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
        nodeX[i] = target_xmin + (tt[i] + 1.0)/2.0 * (target_xmax - target_xmin);
        nodeY[i] = target_ymin + (tt[i] + 1.0)/2.0 * (target_ymax - target_ymin);
        nodeZ[i] = target_zmin + (tt[i] + 1.0)/2.0 * (target_zmax - target_zmin);
    }
    ww[0] = 0.25 * (interp_order*interp_order/3.0 + 1.0/6.0);
    ww[interp_order] = -ww[0];
    
    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interp_order_lim; j++){
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (target_xmax - target_xmin);
        wy[j] = -4.0 * ww[j] / (target_ymax - target_ymin);
        wz[j] = -4.0 * ww[j] / (target_zmax - target_zmin);
    }
    dj[0] = 0.25;
    dj[interp_order] = 0.25;

#ifdef OPENACC_ENABLED
    #pragma acc loop collapse(3) independent
#endif
    for (int ix = target_x_low_ind; ix <= target_x_high_ind; ix++) {
        for (int iy = target_y_low_ind; iy <= target_y_high_ind; iy++) {
            for (int iz = target_z_low_ind; iz <= target_z_high_ind; iz++) {

                int ii = (ix * target_yz_dim) + (iy * target_z_dim_glob) + iz;

                double tx = target_xmin + (ix - target_x_low_ind) * target_xdd;
                double ty = target_ymin + (iy - target_y_low_ind) * target_ydd;
                double tz = target_zmin + (iz - target_z_low_ind) * target_zdd;

                double sumX = 0.0;
                double sumY = 0.0;
                double sumZ = 0.0;

                int eix = -1;
                int eiy = -1;
                int eiz = -1;

#ifdef OPENACC_ENABLED
                #pragma acc loop independent reduction(+:sumX,sumY,sumZ) reduction(max:eix,eiy,eiz)
#endif
                for (int j = 0; j < interp_order_lim; j++) {  // loop through the degree

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
                #pragma acc loop vector independent reduction(+:temp)
#endif
                for (int j = 0; j < interp_pts_per_cluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

                    int k1 = j%interp_order_lim;
                    int kk = (j-k1)/interp_order_lim;
                    int k2 = kk%interp_order_lim;
                    kk = kk - k2;
                    int k3 = kk / interp_order_lim;
            
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
                potential[ii] += temp;
            }
        }
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
