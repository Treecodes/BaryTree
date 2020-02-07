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
#include "struct_nodes.h"
#include "struct_particles.h"
#include "tools.h"

#include "interaction_lists.h"


void cp_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int interpolationOrder, double sizeCheckFactor);



void Interaction_CP_MakeListRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                   int *direct_list, int interpolationOrder, double sizeCheckFactor)
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


    for (int i = 0; i < batch_numnodes; i++) direct_list[i] = -1;
    
    for (int i = 0; i < batch_numnodes; i++) {
        cp_compute_interaction_list_remote(0,
                tree_numpar, tree_radius,
                tree_x_mid, tree_y_mid, tree_z_mid,
                tree_num_children, tree_children,

                &(batch_num_direct[i]), &(batch_num_approx[i]),
                batch_radius[i], batch_x_mid[i], batch_y_mid[i], batch_z_mid[i],

                interpolationOrder, sizeCheckFactor);
                
        if ((batch_num_direct[i] > 0) || (batch_num_approx[i] > 0)) direct_list[i] = 1;

    }

    return;

} /* END of function pc_treecode */



/**********************************************/
/************ LOCAL FUNCTIONS *****************/
/**********************************************/


void cp_compute_interaction_list_remote(int tree_node, const int *tree_numpar, const double *tree_radius,
                const double *tree_x_mid, const double *tree_y_mid, const double *tree_z_mid,
                const int *tree_num_children, const int *tree_children,

                int *batch_num_direct, int *batch_num_approx,
                double batch_radius, double batch_x_mid, double batch_y_mid, double batch_z_mid,

                int interpolationOrder, double sizeCheckFactor)
{

    /* determine DIST for MAC test */
    double tx = batch_x_mid - tree_x_mid[tree_node];
    double ty = batch_y_mid - tree_y_mid[tree_node];
    double tz = batch_z_mid - tree_z_mid[tree_node];
    double dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((tree_radius[tree_node] + batch_radius) < dist * sqrt(thetasq))
      && (tree_radius[tree_node] != 0.00) //) {
      && (sizeCheckFactor*(interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1) < tree_numpar[tree_node])) {
    /*
 *      * If MAC is accepted and there is more than 1 particle
 *           * in the box, use the expansion for the approximation.
 *                */

        (*tree_index_counter)++;

    } else {
    /*
 *      * If MAC fails check to see if there are children. If not, perform direct
 *           * calculation. If there are children, call routine recursively for each.
 *                */
        if (tree_num_children[tree_node] == 0) {

            (*direct_index_counter)++;

        } else {
            for (int i = 0; i < tree_num_children[tree_node]; i++) {
                cp_compute_interaction_list_remote(tree_children[8*tree_node + i],
                           tree_numpar, tree_radius,
                           tree_x_mid, tree_y_mid, tree_z_mid,
                           tree_num_children, tree_children,

                           batch_num_direct, batch_num_approx,
                           batch_radius, batch_x_mid, batch_y_mid, batch_z_mid,

                           interpolationOrder, sizeCheckFactor);
            }
        }
    }

    return;

}
