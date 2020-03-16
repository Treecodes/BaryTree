/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/tools.h"
#include "../utilities/array.h"

#include "../tree/struct_nodes.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "interaction_lists.h"


void pc_compute_interaction_list(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     

                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter,
                struct RunParams *run_params);
                
                
void InteractionList_Make(const struct tnode_array *tree_array,
                          struct tnode_array *batches,
                          int ***approx_inter_list_addr, int ***direct_inter_list_addr,
                          struct RunParams *run_params)
{

    int batch_numnodes = batches->numnodes;
    const int *batch_numpar = batches->numpar;
    const double *batch_x_mid = batches->x_mid;
    const double *batch_y_mid = batches->y_mid;
    const double *batch_z_mid = batches->z_mid;
    const double *batch_radius = batches->radius;

    int *num_approx_inter = batches->numApprox;
    int *num_direct_inter = batches->numDirect;

    int tree_numnodes = tree_array->numnodes;
    const int *tree_numpar = tree_array->numpar;
    const int *tree_level = tree_array->level;
    const double *tree_radius = tree_array->radius;
    const double *tree_x_mid = tree_array->x_mid;
    const double *tree_y_mid = tree_array->y_mid;
    const double *tree_z_mid = tree_array->z_mid;
    
    const int *tree_num_children = tree_array->num_children;
    const int *tree_children = tree_array->children;
    
    
    make_matrix(*approx_inter_list_addr, batch_numnodes, 50);
    make_matrix(*direct_inter_list_addr, batch_numnodes, 50);
    int **approx_inter_list = *approx_inter_list_addr;
    int **direct_inter_list = *direct_inter_list_addr;


    int *sizeof_approx_inter_list, *sizeof_direct_inter_list;
    make_vector(sizeof_approx_inter_list, batch_numnodes);
    make_vector(sizeof_direct_inter_list, batch_numnodes);
    

    for (int i = 0; i < batch_numnodes; i++) sizeof_approx_inter_list[i] = 50;
    for (int i = 0; i < batch_numnodes; i++) sizeof_direct_inter_list[i] = 50;
    
    for (int i = 0; i < batch_numnodes; i++)
        for(int j = 0; j < 50; j++)
            approx_inter_list[i][j] = -1;

    for (int i = 0; i < batch_numnodes; i++)
        for(int j = 0; j < 50; j++)
            direct_inter_list[i][j] = -1;
            
    for (int i = 0; i < batch_numnodes; i++) num_approx_inter[i] = 0;
    for (int i = 0; i < batch_numnodes; i++) num_direct_inter[i] = 0;
    
    
    for (int i = 0; i < batch_numnodes; i++)
        pc_compute_interaction_list(
                    0, tree_numpar, tree_radius,
                    tree_x_mid, tree_y_mid, tree_z_mid,
                    tree_num_children, tree_children,

                    batch_radius[i], batch_x_mid[i], batch_y_mid[i], batch_z_mid[i],

                    &(approx_inter_list[i]), &(direct_inter_list[i]),
                    &(sizeof_approx_inter_list[i]), &(sizeof_direct_inter_list[i]),
                    &(num_approx_inter[i]), &(num_direct_inter[i]),
                    run_params);
                    
    free_vector(sizeof_approx_inter_list);
    free_vector(sizeof_direct_inter_list);
    

    return;

} /* END of function Interaction_MakeList */



void InteractionList_PC_MakeRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                                struct RunParams *run_params)
{
    int batch_numnodes = batches->numnodes;
    const int *batch_numpar = batches->numpar;
    const double *batch_x_mid = batches->x_mid;
    const double *batch_y_mid = batches->y_mid;
    const double *batch_z_mid = batches->z_mid;
    const double *batch_radius = batches->radius;


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


    int **temp_approx_inter_list, **temp_direct_inter_list;
    int *sizeof_approx_inter_list, *sizeof_direct_inter_list;
    int *num_approx_inter, *num_direct_inter;

    make_matrix(temp_approx_inter_list, batch_numnodes, 50);
    make_matrix(temp_direct_inter_list, batch_numnodes, 50);

    make_vector(num_approx_inter, batch_numnodes);
    make_vector(num_direct_inter, batch_numnodes);

    make_vector(sizeof_approx_inter_list, batch_numnodes);
    make_vector(sizeof_direct_inter_list, batch_numnodes);


    for (int i = 0; i < batch_numnodes; i++) sizeof_approx_inter_list[i] = 50;
    for (int i = 0; i < batch_numnodes; i++) sizeof_direct_inter_list[i] = 50;
    
    for (int i = 0; i < batch_numnodes; i++)
        for(int j = 0; j < 50; j++)
            temp_approx_inter_list[i][j] = -1;

    for (int i = 0; i < batch_numnodes; i++)
        for(int j = 0; j < 50; j++)
            temp_direct_inter_list[i][j] = -1;
            
    for (int i = 0; i < batch_numnodes; i++) num_approx_inter[i] = 0;
    for (int i = 0; i < batch_numnodes; i++) num_direct_inter[i] = 0;

    
    for (int i = 0; i < batches->numnodes; i++) {

        pc_compute_interaction_list(
                    0, tree_numpar, tree_radius,
                    tree_x_mid, tree_y_mid, tree_z_mid,
                    tree_num_children, tree_children,

                    batch_radius[i], batch_x_mid[i], batch_y_mid[i], batch_z_mid[i],

                    &(temp_approx_inter_list[i]), &(temp_direct_inter_list[i]),
                    &(sizeof_approx_inter_list[i]), &(sizeof_direct_inter_list[i]),
                    &(num_approx_inter[i]), &(num_direct_inter[i]),
                    run_params);

    }


    for (int j = 0; j < batch_numnodes; j++) {
        for (int i = 0; i < num_approx_inter[j]; i++) {

            int node_index = temp_approx_inter_list[j][i];
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


    free_matrix(temp_approx_inter_list);
    free_matrix(temp_direct_inter_list);

    free_vector(sizeof_approx_inter_list);
    free_vector(sizeof_direct_inter_list);
    
    free_vector(num_approx_inter);
    free_vector(num_direct_inter);

    return;

} /* END of function Interaction_MakeListRemote */



/**********************************************/
/************ LOCAL FUNCTIONS *****************/
/**********************************************/


void pc_compute_interaction_list(
                int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     

                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter,
                struct RunParams *run_params)
{


    /* determine DIST for MAC test */
    double tx = batch_x_mid - tree_x_mid[tree_node];
    double ty = batch_y_mid - tree_y_mid[tree_node];
    double tz = batch_z_mid - tree_z_mid[tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);


    if (((tree_radius[tree_node] + batch_radius) < dist * run_params->theta)
      && (tree_radius[tree_node] != 0.00)
      && (run_params->size_check_factor * run_params->interp_pts_per_cluster < tree_numpar[tree_node])) {
    /*
     * If MAC is accepted use the expansion for the approximation.
     */

        if (*tree_index_counter >= *sizeof_tree_list) {
            (*sizeof_tree_list) *= 1.5;
            (*batch_tree_list) = realloc_vector(*batch_tree_list, *sizeof_tree_list);
        }

        (*batch_tree_list)[*tree_index_counter] = tree_node;
        (*tree_index_counter)++;

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (tree_num_children[tree_node] == 0) {

            if (*direct_index_counter >= *sizeof_direct_list) {
                (*sizeof_direct_list) *= 1.5;
                *batch_direct_list = realloc_vector(*batch_direct_list, *sizeof_direct_list);
            }

            (*batch_direct_list)[*direct_index_counter] = tree_node; 
            (*direct_index_counter)++;

        } else {
            for (int i = 0; i < tree_num_children[tree_node]; i++) {
                pc_compute_interaction_list(tree_children[8*tree_node + i],
                           tree_numpar, tree_radius,
                           tree_x_mid, tree_y_mid, tree_z_mid,
                           tree_num_children, tree_children,     

                           batch_radius, batch_x_mid, batch_y_mid, batch_z_mid,

                           batch_tree_list, batch_direct_list,
                           sizeof_tree_list, sizeof_direct_list,
                           tree_index_counter, direct_index_counter,
                           run_params);
            }
        }
    }

    return;

} 
