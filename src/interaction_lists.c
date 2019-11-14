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
#include "nodes_struct.h"
#include "particles_struct.h"
#include "tools.h"

#include "interaction_lists.h"


void pc_compute_interaction_list(int tree_numnodes, const int *tree_level, 
                const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int *batch_tree_list, int *batch_direct_list);


void pc_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter);



void Interaction_MakeList(const struct tnode_array *tree_array, struct tnode_array *batches,
                          int *tree_inter_list, int *direct_inter_list, int approx_offset,
                          int direct_offset)
{
    int batch_numnodes = batches->numnodes;
    const int *batch_numpar = batches->numpar;
    const double *batch_x_mid = batches->x_mid;
    const double *batch_y_mid = batches->y_mid;
    const double *batch_z_mid = batches->z_mid;
    const double *batch_radius = batches->radius;

    int *batch_num_direct = batches->numDirect;
    int *batch_num_approx = batches->numApprox;

    int tree_numnodes = tree_array->numnodes;
    const int *tree_numpar = tree_array->numpar;
    const int *tree_level = tree_array->level;
    const double *tree_radius = tree_array->radius;
    const double *tree_x_mid = tree_array->x_mid;
    const double *tree_y_mid = tree_array->y_mid;
    const double *tree_z_mid = tree_array->z_mid;

    int sizeloop = batches->numnodes * approx_offset;
    for (int i = 0; i < sizeloop; i++) tree_inter_list[i] = -1;

    sizeloop = batches->numnodes * direct_offset;
    for (int i = 0; i < sizeloop; i++) direct_inter_list[i] = -1;
    
    for (int i = 0; i < batches->numnodes; i++)
        pc_compute_interaction_list(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,

                &(batch_num_direct[i]), &(batch_num_approx[i]),
                batch_radius[i], batch_x_mid[i], batch_y_mid[i], batch_z_mid[i],

                &(tree_inter_list[i*approx_offset]), &(direct_inter_list[i*direct_offset]));

    return;

} /* END of function pc_treecode */



void Interaction_MakeListRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                                int *num_batch_approx, int *num_batch_direct)
{
    int batch_numnodes = batches->numnodes;
    const int *batch_numpar = batches->numpar;
    const double *batch_x_mid = batches->x_mid;
    const double *batch_y_mid = batches->y_mid;
    const double *batch_z_mid = batches->z_mid;
    const double *batch_radius = batches->radius;

    int *batch_num_direct = batches->numDirect;
    int *batch_num_approx = batches->numApprox;

    int tree_numnodes = tree_array->numnodes;
    const int *tree_numpar = tree_array->numpar;
    const int *tree_level = tree_array->level;
    const double *tree_radius = tree_array->radius;
    const double *tree_x_mid = tree_array->x_mid;
    const double *tree_y_mid = tree_array->y_mid;
    const double *tree_z_mid = tree_array->z_mid;

    const int *tree_num_children = tree_array->num_children;
    const int *tree_children = tree_array->children;


    for (int i = 0; i < tree_numnodes; i++) approx_list_unpacked[i] = -1;
    for (int i = 0; i < tree_numnodes; i++) approx_list_packed[i] = -1;
    for (int i = 0; i < tree_numnodes; i++) direct_list[i] = -1;


    int **temp_tree_inter_list, **temp_direct_inter_list;
    int *num_tree_inter, *num_direct_inter;
    int *sizeof_tree_inter_list, *sizeof_direct_inter_list;

    make_matrix(temp_tree_inter_list, batches->numnodes, 50);
    make_matrix(temp_direct_inter_list, batches->numnodes, 50);

    make_vector(num_tree_inter, batches->numnodes);
    make_vector(num_direct_inter, batches->numnodes);

    make_vector(sizeof_tree_inter_list, batches->numnodes);
    make_vector(sizeof_direct_inter_list, batches->numnodes);
   
    int loopsize = batches->numnodes;
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
    for (int i = 0; i < batches->numnodes; i++) {

        pc_compute_interaction_list_remote(0,
                tree_numpar, tree_radius,
                tree_x_mid, tree_y_mid, tree_z_mid,
                tree_num_children, tree_children,     

                &(batch_num_direct[i]), &(batch_num_approx[i]),
                batch_radius[i], batch_x_mid[i], batch_y_mid[i], batch_z_mid[i],

                &(temp_tree_inter_list[i]), &(temp_direct_inter_list[i]),
                &(sizeof_tree_inter_list[i]), &(sizeof_direct_inter_list[i]),
                &(num_tree_inter[i]), &(num_direct_inter[i]));

        num_batch_approx[i] += num_tree_inter[i];
        num_batch_direct[i] += num_direct_inter[i];

    }


    for (int j = 0; j < batches->numnodes; j++) {
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

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int *batch_tree_list, int *batch_direct_list)
{
    int tree_index_counter = 0;
    int direct_index_counter = 0;
    int current_level = 0;
    
    for (int j = 0; j < tree_numnodes; j++) {
        if (tree_level[j] <= current_level) {
            
            /* determine DIST for MAC test */
            double tx = batch_x_mid - tree_x_mid[j];
            double ty = batch_y_mid - tree_y_mid[j];
            double tz = batch_z_mid - tree_z_mid[j];
            double dist = sqrt(tx*tx + ty*ty + tz*tz);

            if (((tree_radius[j] + batch_radius) < dist * sqrt(thetasq))
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
    (*batch_num_direct) = direct_index_counter;
    (*batch_num_approx) = tree_index_counter;
    
    return;
}



void pc_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter)
{

    /* determine DIST for MAC test */
    double tx = batch_x_mid - tree_x_mid[tree_node];
    double ty = batch_y_mid - tree_y_mid[tree_node];
    double tz = batch_z_mid - tree_z_mid[tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((tree_radius[tree_node] + batch_radius) < dist * sqrt(thetasq))
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

                           batch_num_direct, batch_num_approx,
                           batch_radius, batch_x_mid, batch_y_mid, batch_z_mid,

                           batch_tree_list, batch_direct_list,
                           sizeof_tree_list, sizeof_direct_list,
                           tree_index_counter, direct_index_counter);
            }
        }
    }

    return;

} 
