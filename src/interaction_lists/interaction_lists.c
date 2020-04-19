#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../utilities/tools.h"
#include "../utilities/array.h"

#include "../tree/struct_tree.h"
#include "../run_params/struct_run_params.h"

#include "struct_interaction_lists.h"
#include "interaction_lists.h"
                

void pc_compute_interaction_list(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,     

                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int **batch_tree_list, int **batch_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter,
                const struct RunParams *run_params);
                

void cc_compute_interaction_list_1(
                int source_tree_node, const int *source_tree_numpar, const double *source_tree_radius,
                const double *source_tree_x_mid, const double *source_tree_y_mid, const double *source_tree_z_mid,
                const int *source_tree_num_children, const int *source_tree_children,

                int target_tree_node, const int *target_tree_numpar, const double *target_tree_radius,
                const double *target_tree_x_mid, const double *target_tree_y_mid, const double *target_tree_z_mid,
                const int *target_tree_num_children, const int *target_tree_children,

                int **target_approx_list, int **target_direct_list,
                int *sizeof_approx_list, int *sizeof_direct_list,
                int *approx_index_counter, int *direct_index_counter,
                const struct RunParams *run_params);


void cc_compute_interaction_list_2(
                int target_tree_node, const int *target_tree_numpar, const double *target_tree_radius,
                const double *target_tree_x_mid, const double *target_tree_y_mid, const double *target_tree_z_mid,
                const int *target_tree_num_children, const int *target_tree_children,

                int source_tree_node, const int *source_tree_numpar, const double *source_tree_radius,
                const double *source_tree_x_mid, const double *source_tree_y_mid, const double *source_tree_z_mid,
                const int *source_tree_num_children, const int *source_tree_children,

                int **target_approx_list, int **target_direct_list,
                int *sizeof_approx_list, int *sizeof_direct_list,
                int *approx_index_counter, int *direct_index_counter,
                const struct RunParams *run_params);
                
                
                
void InteractionLists_Make(struct InteractionLists **interaction_list_addr,
                          const struct Tree *source_tree,
                          const struct Tree *target_tree,
                          const struct RunParams *run_params)
{

    *interaction_list_addr = malloc(sizeof(struct InteractionLists));
    struct InteractionLists *interaction_list = *interaction_list_addr;

    int ***approx_inter_list_addr = &(interaction_list->approx_interactions);
    int ***direct_inter_list_addr = &(interaction_list->direct_interactions);
    
    int **num_approx_addr = &(interaction_list->num_approx);
    int **num_direct_addr = &(interaction_list->num_direct);
    
    
   int source_tree_numnodes = source_tree->numnodes;
   const int *source_tree_numpar = source_tree->numpar;
   const double *source_tree_radius = source_tree->radius;
   const double *source_tree_x_mid = source_tree->x_mid;
   const double *source_tree_y_mid = source_tree->y_mid;
   const double *source_tree_z_mid = source_tree->z_mid;

   const int *source_tree_num_children = source_tree->num_children;
   const int *source_tree_children = source_tree->children;
   
   
   int target_tree_numnodes = target_tree->numnodes;
   const int *target_tree_numpar = target_tree->numpar;
   const double *target_tree_radius = target_tree->radius;
   const double *target_tree_x_mid = target_tree->x_mid;
   const double *target_tree_y_mid = target_tree->y_mid;
   const double *target_tree_z_mid = target_tree->z_mid;
   
   const int *target_tree_num_children = target_tree->num_children;
   const int *target_tree_children = target_tree->children;
   

   make_matrix(*approx_inter_list_addr, target_tree_numnodes, 50);
   make_matrix(*direct_inter_list_addr, target_tree_numnodes, 50);
   int **approx_inter_list = *approx_inter_list_addr;
   int **direct_inter_list = *direct_inter_list_addr;
   
   make_vector(*num_approx_addr, target_tree_numnodes);
   make_vector(*num_direct_addr, target_tree_numnodes);
   int *num_approx_inter = *num_approx_addr;
   int *num_direct_inter = *num_direct_addr;
   

   int *sizeof_approx_inter_list, *sizeof_direct_inter_list;
   make_vector(sizeof_approx_inter_list, target_tree_numnodes);
   make_vector(sizeof_direct_inter_list, target_tree_numnodes);
   

   for (int i = 0; i < target_tree_numnodes; i++) sizeof_approx_inter_list[i] = 50;
   for (int i = 0; i < target_tree_numnodes; i++) sizeof_direct_inter_list[i] = 50;
   
   for (int i = 0; i < target_tree_numnodes; i++)
       for (int j = 0; j < 50; j++)
           approx_inter_list[i][j] = -1;

   for (int i = 0; i < target_tree_numnodes; i++)
       for (int j = 0; j < 50; j++)
           direct_inter_list[i][j] = -1;
           
   for (int i = 0; i < target_tree_numnodes; i++) num_approx_inter[i] = 0;
   for (int i = 0; i < target_tree_numnodes; i++) num_direct_inter[i] = 0;
    
    
    if (run_params->compute_type == PARTICLE_CLUSTER || run_params->compute_type == CLUSTER_PARTICLE) {
    
        for (int i = 0; i < target_tree_numnodes; i++) {
            pc_compute_interaction_list(
                    0, source_tree_numpar, source_tree_radius,
                    source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                    source_tree_num_children, source_tree_children,

                    target_tree_radius[i], target_tree_x_mid[i], target_tree_y_mid[i], target_tree_z_mid[i],

                    &(approx_inter_list[i]), &(direct_inter_list[i]),
                    &(sizeof_approx_inter_list[i]), &(sizeof_direct_inter_list[i]),
                    &(num_approx_inter[i]), &(num_direct_inter[i]),
                    run_params);
        }
    
    } else if (run_params->compute_type == CLUSTER_CLUSTER) {
    
        cc_compute_interaction_list_1(
                    0, source_tree_numpar, source_tree_radius,
                    source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                    source_tree_num_children, source_tree_children,

                    0, target_tree_numpar, target_tree_radius,
                    target_tree_x_mid, target_tree_y_mid, target_tree_z_mid,
                    target_tree_num_children, target_tree_children,

                    approx_inter_list, direct_inter_list,
                    sizeof_approx_inter_list, sizeof_direct_inter_list,
                    num_approx_inter, num_direct_inter,
                    run_params);
    }
                    
    free_vector(sizeof_approx_inter_list);
    free_vector(sizeof_direct_inter_list);
    

    return;

} /* END of function Interaction_MakeList */



void InteractionLists_Free(struct InteractionLists **interaction_list_addr)
{
    struct InteractionLists *interaction_list = *interaction_list_addr;

    free_matrix(interaction_list->approx_interactions);
    free_matrix(interaction_list->direct_interactions);
    free_vector(interaction_list->num_approx);
    free_vector(interaction_list->num_direct);
    free(interaction_list);
    
    interaction_list = NULL;

    return;
}



void InteractionLists_MakeRemote(const struct Tree *source_tree,
                                 const struct Tree *target_tree,
                                 int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                                 const struct RunParams *run_params)
{
    int source_tree_numnodes = source_tree->numnodes;
    const int *source_tree_numpar = source_tree->numpar;
    const double *source_tree_radius = source_tree->radius;
    const double *source_tree_x_mid = source_tree->x_mid;
    const double *source_tree_y_mid = source_tree->y_mid;
    const double *source_tree_z_mid = source_tree->z_mid;

    const int *source_tree_num_children = source_tree->num_children;
    const int *source_tree_children = source_tree->children;
    
    
    int target_tree_numnodes = target_tree->numnodes;
    const int *target_tree_numpar = target_tree->numpar;
    const double *target_tree_radius = target_tree->radius;
    const double *target_tree_x_mid = target_tree->x_mid;
    const double *target_tree_y_mid = target_tree->y_mid;
    const double *target_tree_z_mid = target_tree->z_mid;
    
    const int *target_tree_num_children = target_tree->num_children;
    const int *target_tree_children = target_tree->children;


    for (int i = 0; i < source_tree_numnodes; i++) approx_list_unpacked[i] = -1;
    for (int i = 0; i < source_tree_numnodes; i++) approx_list_packed[i] = -1;
    for (int i = 0; i < source_tree_numnodes; i++) direct_list[i] = -1;


    int **temp_approx_inter_list, **temp_direct_inter_list;
    int *sizeof_approx_inter_list, *sizeof_direct_inter_list;
    int *num_approx_inter, *num_direct_inter;

    make_matrix(temp_approx_inter_list, target_tree_numnodes, 50);
    make_matrix(temp_direct_inter_list, target_tree_numnodes, 50);

    make_vector(sizeof_approx_inter_list, target_tree_numnodes);
    make_vector(sizeof_direct_inter_list, target_tree_numnodes);
    
    make_vector(num_approx_inter, target_tree_numnodes);
    make_vector(num_direct_inter, target_tree_numnodes);


    for (int i = 0; i < target_tree_numnodes; i++) sizeof_approx_inter_list[i] = 50;
    for (int i = 0; i < target_tree_numnodes; i++) sizeof_direct_inter_list[i] = 50;
    
    for (int i = 0; i < target_tree_numnodes; i++)
        for(int j = 0; j < 50; j++)
            temp_approx_inter_list[i][j] = -1;

    for (int i = 0; i < target_tree_numnodes; i++)
        for(int j = 0; j < 50; j++)
            temp_direct_inter_list[i][j] = -1;
            
    for (int i = 0; i < target_tree_numnodes; i++) num_approx_inter[i] = 0;
    for (int i = 0; i < target_tree_numnodes; i++) num_direct_inter[i] = 0;


    if (run_params->compute_type == PARTICLE_CLUSTER) {
        
    for (int i = 0; i < target_tree_numnodes; i++) {
        pc_compute_interaction_list(
                    0, source_tree_numpar, source_tree_radius,
                    source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                    source_tree_num_children, source_tree_children,

                    target_tree_radius[i], target_tree_x_mid[i], target_tree_y_mid[i], target_tree_z_mid[i],

                    &(temp_approx_inter_list[i]), &(temp_direct_inter_list[i]),
                    &(sizeof_approx_inter_list[i]), &(sizeof_direct_inter_list[i]),
                    &(num_approx_inter[i]), &(num_direct_inter[i]),
                    run_params);
        }

    } else if (run_params->compute_type == CLUSTER_CLUSTER) {
    
        cc_compute_interaction_list_1(
                    0, source_tree_numpar, source_tree_radius,
                    source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                    source_tree_num_children, source_tree_children,

                    0, target_tree_numpar, target_tree_radius,
                    target_tree_x_mid, target_tree_y_mid, target_tree_z_mid,
                    target_tree_num_children, target_tree_children,

                    temp_approx_inter_list, temp_direct_inter_list,
                    sizeof_approx_inter_list, sizeof_direct_inter_list,
                    num_approx_inter, num_direct_inter,
                    run_params);
    }


    for (int j = 0; j < target_tree_numnodes; j++) {
        for (int i = 0; i < num_approx_inter[j]; i++) {

            int source_node_index = temp_approx_inter_list[j][i];
            approx_list_unpacked[source_node_index] = source_node_index;
        }

        for (int i = 0; i < num_direct_inter[j]; i++) {

            int source_node_index = temp_direct_inter_list[j][i];
            direct_list[source_node_index] = source_node_index;
        }
    }


    int approx_counter = 0;
    for (int i = 0; i < source_tree_numnodes; i++) {
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
                const struct RunParams *run_params)
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



void cc_compute_interaction_list_1(
                int source_tree_node, const int *source_tree_numpar, const double *source_tree_radius,
                const double *source_tree_x_mid, const double *source_tree_y_mid, const double *source_tree_z_mid,
                const int *source_tree_num_children, const int *source_tree_children,

                int target_tree_node, const int *target_tree_numpar, const double *target_tree_radius,
                const double *target_tree_x_mid, const double *target_tree_y_mid, const double *target_tree_z_mid,
                const int *target_tree_num_children, const int *target_tree_children,

                int **target_tree_list, int **target_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter,
                const struct RunParams *run_params)
{

    /* determine DIST for MAC test */
    double tx = target_tree_x_mid[target_tree_node] - source_tree_x_mid[source_tree_node];
    double ty = target_tree_y_mid[target_tree_node] - source_tree_y_mid[source_tree_node];
    double tz = target_tree_z_mid[target_tree_node] - source_tree_z_mid[source_tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((source_tree_radius[source_tree_node] + target_tree_radius[target_tree_node]) < dist * run_params->theta)
      && (source_tree_radius[source_tree_node] != 0.00)
      && (pow(run_params->size_check_factor * run_params->interp_pts_per_cluster, 2)
          < source_tree_numpar[source_tree_node] * target_tree_numpar[target_tree_node])) {
   /*
    * If MAC is accepted and there is more than 1 particle
    * in the box, use the expansion for the approximation.
    */

        if (tree_index_counter[target_tree_node] >= sizeof_tree_list[target_tree_node]) {
            sizeof_tree_list[target_tree_node] *= 1.5;
            target_tree_list[target_tree_node] = realloc_vector(target_tree_list[target_tree_node],
                                                                sizeof_tree_list[target_tree_node]);
        }

        target_tree_list[target_tree_node][tree_index_counter[target_tree_node]] = source_tree_node;
        tree_index_counter[target_tree_node]++;

    } else {
   /*
    * If MAC fails check to see if there are children. If not, perform direct
    * calculation. If there are children, call routine recursively for each.
    */
        if (target_tree_num_children[target_tree_node] == 0) {

            if (direct_index_counter[target_tree_node] >= sizeof_direct_list[target_tree_node]) {
                sizeof_direct_list[target_tree_node] *= 1.5;
                target_direct_list[target_tree_node] = realloc_vector(target_direct_list[target_tree_node],
                                                                      sizeof_direct_list[target_tree_node]);
            }

            target_direct_list[target_tree_node][direct_index_counter[target_tree_node]] = source_tree_node;
            direct_index_counter[target_tree_node]++;

        } else {
            for (int i = 0; i < target_tree_num_children[target_tree_node]; i++) {
                cc_compute_interaction_list_2(target_tree_children[8*target_tree_node + i],
                           target_tree_numpar, target_tree_radius,
                           target_tree_x_mid, target_tree_y_mid, target_tree_z_mid,
                           target_tree_num_children, target_tree_children,

                           source_tree_node, source_tree_numpar, source_tree_radius,
                           source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                           source_tree_num_children, source_tree_children,

                           target_tree_list, target_direct_list,
                           sizeof_tree_list, sizeof_direct_list,
                           tree_index_counter, direct_index_counter,
                           run_params);
            }
        }
    }

    return;

}



void cc_compute_interaction_list_2(
                int target_tree_node, const int *target_tree_numpar, const double *target_tree_radius,
                const double *target_tree_x_mid, const double *target_tree_y_mid, const double *target_tree_z_mid,
                const int *target_tree_num_children, const int *target_tree_children,

                int source_tree_node, const int *source_tree_numpar, const double *source_tree_radius,
                const double *source_tree_x_mid, const double *source_tree_y_mid, const double *source_tree_z_mid,
                const int *source_tree_num_children, const int *source_tree_children,

                int **target_tree_list, int **target_direct_list,
                int *sizeof_tree_list, int *sizeof_direct_list,
                int *tree_index_counter, int *direct_index_counter,
                const struct RunParams *run_params)
{

    /* determine DIST for MAC test */
    double tx = target_tree_x_mid[target_tree_node] - source_tree_x_mid[source_tree_node];
    double ty = target_tree_y_mid[target_tree_node] - source_tree_y_mid[source_tree_node];
    double tz = target_tree_z_mid[target_tree_node] - source_tree_z_mid[source_tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((source_tree_radius[source_tree_node] + target_tree_radius[target_tree_node]) < dist * run_params->theta)
      && (target_tree_radius[source_tree_node] != 0.00)
      && (pow(run_params->size_check_factor * run_params->interp_pts_per_cluster, 2)
          < source_tree_numpar[source_tree_node] * target_tree_numpar[target_tree_node])) {
   /*
    * If MAC is accepted and there is more than 1 particle
    * in the box, use the expansion for the approximation.
    */

        if (tree_index_counter[target_tree_node] >= sizeof_tree_list[target_tree_node]) {
            sizeof_tree_list[target_tree_node] *= 1.5;
            target_tree_list[target_tree_node] = realloc_vector(target_tree_list[target_tree_node],
                                                                sizeof_tree_list[target_tree_node]);
        }

        target_tree_list[target_tree_node][tree_index_counter[target_tree_node]] = source_tree_node;
        tree_index_counter[target_tree_node]++;

    } else {
   /*
    * If MAC fails check to see if there are children. If not, perform direct
    * calculation. If there are children, call routine recursively for each.
    */
        if (source_tree_num_children[source_tree_node] == 0) {

            if (direct_index_counter[target_tree_node] >= sizeof_direct_list[target_tree_node]) {
                sizeof_direct_list[target_tree_node] *= 1.5;
                target_direct_list[target_tree_node] = realloc_vector(target_direct_list[target_tree_node],
                                                                      sizeof_direct_list[target_tree_node]);
            }

            target_direct_list[target_tree_node][direct_index_counter[target_tree_node]] = source_tree_node;
            direct_index_counter[target_tree_node]++;

        } else {
            for (int i = 0; i < source_tree_num_children[source_tree_node]; i++) {
                cc_compute_interaction_list_1(source_tree_children[8*source_tree_node + i],
                           source_tree_numpar, source_tree_radius,
                           source_tree_x_mid, source_tree_y_mid, source_tree_z_mid,
                           source_tree_num_children, source_tree_children,
                           
                           target_tree_node, target_tree_numpar, target_tree_radius,
                           target_tree_x_mid, target_tree_y_mid, target_tree_z_mid,
                           target_tree_num_children, target_tree_children,

                           target_tree_list, target_direct_list,
                           sizeof_tree_list, sizeof_direct_list,
                           tree_index_counter, direct_index_counter,
                           run_params);
            }
        }
    }

    return;

}
