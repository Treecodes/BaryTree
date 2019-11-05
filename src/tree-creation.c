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

                pc_create_tree_n0(paddress,
                                  sources, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev);

            }
        }
        
    } else {
   
        {     
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

    tree_array->num_children[p->node_index] = p->num_children;

    for (i = 0; i < p->num_children; i++) {
        tree_array->children[8*p->node_index+i] = (p->child[i])->node_index;
        pc_create_tree_array(p->child[i], tree_array);
    }
    
    return;
    
} /* END of function create_tree_n0 */



