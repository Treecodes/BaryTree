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

#include "interaction_lists.h"


void pc_compute_interaction_list(int tree_numnodes, const int *tree_level, 
                const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list);

void pc_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     
                int *batch_ind, double *batch_mid, double batch_rad,
                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter);



void Interaction_MakeList(const struct tnode_array *tree_array, struct batch *batches,
                          int *tree_inter_list, int *direct_inter_list, int approx_offset,
                          int direct_offset)
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

    int sizeloop = batches->num * approx_offset;
    for (i = 0; i < sizeloop; i++) tree_inter_list[i] = -1;

    sizeloop = batches->num * direct_offset;
    for (i = 0; i < sizeloop; i++) direct_inter_list[i] = -1;
    
    for (i = 0; i < batches->num; i++)
        pc_compute_interaction_list(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,
                batches_ind[i], batches_center[i], batches_radius[i],
                &(tree_inter_list[i*approx_offset]), &(direct_inter_list[i*direct_offset]));

    return;

} /* END of function pc_treecode */



void Interaction_MakeListRemote(const struct tnode_array *tree_array, struct batch *batches,
                                int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                                int *num_batch_approx, int *num_batch_direct)
{
    /* local variables */
    int rank, numProcs;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    int i, j;

    int **batches_ind;
    double **batches_center;
    double *batches_radius;

    int tree_numnodes;
    const int *tree_numpar, *tree_level;
    const double *tree_x_mid, *tree_y_mid, *tree_z_mid, *tree_radius;
    const int *tree_num_children, *tree_children;

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

    tree_num_children = tree_array->num_children;
    tree_children = tree_array->children;

    for (int i = 0; i < tree_numnodes; i++) approx_list_unpacked[i] = -1;
    for (int i = 0; i < tree_numnodes; i++) approx_list_packed[i] = -1;
    for (int i = 0; i < tree_numnodes; i++) direct_list[i] = -1;


    int **temp_tree_inter_list, **temp_direct_inter_list;
    int *num_tree_inter, *num_direct_inter;
    int *sizeof_tree_inter_list, *sizeof_direct_inter_list;

    make_matrix(temp_tree_inter_list, batches->num, 50);
    make_matrix(temp_direct_inter_list, batches->num, 50);

    make_vector(num_tree_inter, batches->num);
    make_vector(num_direct_inter, batches->num);

    make_vector(sizeof_tree_inter_list, batches->num);
    make_vector(sizeof_direct_inter_list, batches->num);
   
    int loopsize = batches->num;
    for (int i = 0; i < loopsize; i++) num_tree_inter[i] = 0;
    for (int i = 0; i < loopsize; i++) num_direct_inter[i] = 0;

    for (int i = 0; i < loopsize; i++) sizeof_tree_inter_list[i] = 50;
    for (int i = 0; i < loopsize; i++) sizeof_direct_inter_list[i] = 50;
    
    for (int i = 0; i < loopsize; i++)
        for(int j = 0; j < 50; j++)
            temp_tree_inter_list[i][j] = -1;

    for (int i = 0; i < loopsize; i++)
        for(int j = 0; j < 50; j++)
            temp_direct_inter_list[i][j] = -1;
 
    
    // Fill interaction lists
    for (int i = 0; i < batches->num; i++) {

        pc_compute_interaction_list_remote(0,
                tree_numpar, tree_radius,
                tree_x_mid, tree_y_mid, tree_z_mid,
                tree_num_children, tree_children,     
                batches_ind[i], batches_center[i], batches_radius[i],
                &(temp_tree_inter_list[i]), &(temp_direct_inter_list[i]),
                &(sizeof_tree_inter_list[i]), &(sizeof_direct_inter_list[i]),
                &(num_tree_inter[i]), &(num_direct_inter[i]));

        num_batch_approx[i] += num_tree_inter[i];
        num_batch_direct[i] += num_direct_inter[i];

    }


    for (int j = 0; j < batches->num; j++) {
        for (int i = 0; i < num_tree_inter[j]; i++) {

            int node_index = temp_tree_inter_list[j][i];
            approx_list_unpacked[node_index] = node_index;
        }

        for (int i = 0; i < num_direct_inter[j]; i++) {

            int node_index = temp_direct_inter_list[j][i];
            direct_list[node_index] = node_index;
        }
    }


    int approx_counter = 0;
    for (int i = 0; i < tree_numnodes; i++) {
        if (approx_list_unpacked[i] > -1) {
            approx_list_packed[approx_counter] = i;
            approx_counter++;
        }
    }


    free_matrix(temp_tree_inter_list);
    free_matrix(temp_direct_inter_list);

    free_vector(sizeof_tree_inter_list);
    free_vector(sizeof_direct_inter_list);
    
    free_vector(num_tree_inter);
    free_vector(num_direct_inter);

    return;

} /* END of function pc_treecode */



/**********************************************/
/************ LOCAL FUNCTIONS *****************/
/**********************************************/


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



void pc_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     
                int *batch_ind, double *batch_mid, double batch_rad,
                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter)
{

    /* determine DIST for MAC test */
    double tx = batch_mid[0] - tree_x_mid[tree_node];
    double ty = batch_mid[1] - tree_y_mid[tree_node];
    double tz = batch_mid[2] - tree_z_mid[tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((tree_radius[tree_node] + batch_rad) < dist * sqrt(thetasq))
      && (tree_radius[tree_node] != 0.00)) {
    /*
 *      * If MAC is accepted and there is more than 1 particle
 *           * in the box, use the expansion for the approximation.
 *                */

        if (*tree_index_counter >= *sizeof_tree_list) {
            (*sizeof_tree_list) *= 1.5;
            (*batch_tree_list) = realloc_vector(*batch_tree_list, *sizeof_tree_list);
        }

        (*batch_tree_list)[*tree_index_counter] = tree_node;
        (*tree_index_counter)++;

    } else {
    /*
 *      * If MAC fails check to see if there are children. If not, perform direct
 *           * calculation. If there are children, call routine recursively for each.
 *                */
        if (tree_num_children[tree_node] == 0) {

            if (*direct_index_counter >= *sizeof_direct_list) {
                (*sizeof_direct_list) *= 1.5;
                *batch_direct_list = realloc_vector(*batch_direct_list, *sizeof_direct_list);
            }

            (*batch_direct_list)[*direct_index_counter] = tree_node; 
            (*direct_index_counter)++;

        } else {
            for (int i = 0; i < tree_num_children[tree_node]; i++) {
                pc_compute_interaction_list_remote(tree_children[8*tree_node + i],
                           tree_numpar, tree_radius,
                           tree_x_mid, tree_y_mid, tree_z_mid,
                           tree_num_children, tree_children,     
                           batch_ind, batch_mid, batch_rad,
                           batch_tree_list, batch_direct_list,
                           sizeof_tree_list, sizeof_direct_list,
                           tree_index_counter, direct_index_counter);
            }
        }
    }

    return;

} 
