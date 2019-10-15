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
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



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

    for (i = 0; i < batches->num * tree_numnodes; i++)
        tree_inter_list[i] = -1;

    for (i = 0; i < batches->num * tree_numnodes; i++)
        direct_inter_list[i] = -1;
    
    for (i = 0; i < batches->num; i++)
        pc_compute_interaction_list(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,
                batches_ind[i], batches_center[i], batches_radius[i],
                &(tree_inter_list[i*tree_numnodes]), &(direct_inter_list[i*tree_numnodes]));

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
                && (tree_radius[j] > 0.00)) {
//                && (torder*torder*torder < tree_numpar[j])) {
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

                if ( (j==tree_numnodes-1) || (tree_level[j+1] <= tree_level[j]) ) {
                    
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
                && (tree_radius[j] > 0.00)) {
//                && (torder*torder*torder < tree_numpar[j])) {
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

                if ( (j==tree_numnodes-1) || (tree_level[j+1] <= tree_level[j]) ) {

                    batch_direct_list[j] = j;

                } else {

                    current_level = tree_level[j+1];

                }
            }
        }
    }

    return;
}


