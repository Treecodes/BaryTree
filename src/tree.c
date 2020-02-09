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
#include "partition.h"

#include "tree.h"


double thetasq;
double *tt, *ww;

static void remove_node(struct tnode *p);


void Tree_Setup(struct particles *particles1, struct particles *particles2, 
                int order, double theta, double *xyzminmax)
{
    /* changing values of our extern variables */
    thetasq = theta * theta;
    make_vector(tt, order+1);
    make_vector(ww, order+1);

    /* initializing array for Chev points */
    for (int i = 0; i < order+1; i++)
        tt[i] = cos(i * M_PI / order);

    ww[0] = 0.25 * (order*order/3.0 + 1.0/6.0);
    ww[order] = -ww[0];

    for (int i = 1; i < order; i++) {
        double xx = i * M_PI / order;
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
    }

    int particles1_num = particles1->num;
    int particles2_num = particles2->num;

    /* find bounds of Cartesian box enclosing the particles */
    xyzminmax[0] = minval(particles1->x, particles1_num);
    xyzminmax[1] = maxval(particles1->x, particles1_num);
    xyzminmax[2] = minval(particles1->y, particles1_num);
    xyzminmax[3] = maxval(particles1->y, particles1_num);
    xyzminmax[4] = minval(particles1->z, particles1_num);
    xyzminmax[5] = maxval(particles1->z, particles1_num);

    /* setting up ordering vectors */
    make_vector(particles1->order, particles1_num);
    make_vector(particles2->order, particles2_num);
    for (int i = 0; i < particles1_num; i++) particles1->order[i] = i+1;
    for (int i = 0; i < particles2_num; i++) particles2->order[i] = i+1;

    return;
    
} /* END of function setup */




void Tree_CP_Create(struct tnode **p, struct particles *targets,
                    int ibeg, int iend, int maxparnode,
                    double *xyzmm, int level, int *numnodes, int *numleaves)
{
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }

    (*p) = malloc(sizeof(struct tnode));

    (*numnodes)++;

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;

    (*p)->x_min = minval(targets->x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(targets->x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(targets->y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(targets->y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(targets->z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(targets->z + ibeg - 1, (*p)->numpar);
    

    /*compute aspect ratio*/
    double xl = (*p)->x_max - (*p)->x_min;
    double yl = (*p)->y_max - (*p)->y_min;
    double zl = (*p)->z_max - (*p)->z_min;
        
    double tmax = max3(xl, yl, zl);
    double tmin = min3(xl, yl, zl);


    if (tmin != 0.0)
        (*p)->aspect = tmax/tmin;
    else
        (*p)->aspect = 0.0;

    
    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    double t1 = (*p)->x_max - (*p)->x_mid;
    double t2 = (*p)->y_max - (*p)->y_mid;
    double t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    
    /*set particle limits, tree level of node, and nullify child pointers*/
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    (*p)->num_children = 0;
    for (int i = 0; i < 8; i++)
        (*p)->child[i] = NULL;
    

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

        double x_mid = (*p)->x_mid;
        double y_mid = (*p)->y_mid;
        double z_mid = (*p)->z_mid;
        int numposchild;

        cp_partition_8(targets->x, targets->y, targets->z, targets->q, targets->order,
                       xyzmms, xl, yl, zl, tmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        int loclev = level + 1;

        for (int i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                (*p)->num_children = (*p)->num_children + 1;
                int idx = (*p)->num_children - 1;

                for (int j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                struct tnode **paddress = &((*p)->child[idx]);

                Tree_CP_Create(paddress,
                               targets, ind[i][0], ind[i][1],
                               maxparnode, lxyzmm, loclev, numnodes, numleaves);
            }
        }

    } else {

        (*numleaves)++;
    }

    return;

} /* end of function create_tree_n0 */




void Tree_PC_Create(struct tnode **p, struct particles *sources,
                    int ibeg, int iend, int maxparnode, double *xyzmm,
                    int level, int *numnodes, int *numleaves)
{
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];
    
    
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }
                        

    (*p) = malloc(sizeof(struct tnode));


    /* increment number of nodes */
    (*numnodes)++;

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
    double xl = (*p)->x_max - (*p)->x_min;
    double yl = (*p)->y_max - (*p)->y_min;
    double zl = (*p)->z_max - (*p)->z_min;
        
    double tmax = max3(xl, yl, zl);
    double tmin = min3(xl, yl, zl);


    if (tmin != 0.0)
        (*p)->aspect = tmax/tmin;
    else
        (*p)->aspect = 0.0;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    double t1 = (*p)->x_max - (*p)->x_mid;
    double t2 = (*p)->y_max - (*p)->y_mid;
    double t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    /*set particle limits, tree level of node, and nullify child pointers*/
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    (*p)->num_children = 0;
    for (int i = 0; i < 8; i++)
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

        double x_mid = (*p)->x_mid;
        double y_mid = (*p)->y_mid;
        double z_mid = (*p)->z_mid;
        int numposchild;

        pc_partition_8(sources->x, sources->y, sources->z, sources->q, sources->w, sources->order,
                       xyzmms, xl, yl, zl, tmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        int loclev = level + 1;

        for (int i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                (*p)->num_children = (*p)->num_children + 1;
                int idx = (*p)->num_children - 1;

                for (int j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                struct tnode **paddress = &((*p)->child[idx]);

                Tree_PC_Create(paddress,
                               sources, ind[i][0], ind[i][1],
                               maxparnode, lxyzmm, loclev,
                               numnodes, numleaves);

            }
        }
        
    } else {
   
        {     
            /* increment number of leaves */
            (*numleaves)++;
        }
    }

    return;

} /* END of function create_tree_n0 */




int Tree_SetIndex(struct tnode *p, int index)
{
        int current_index = index;
        p->node_index = current_index;

        for (int i = 0; i < p->num_children; i++)
            current_index = Tree_SetIndex(p->child[i], current_index + 1);

        return current_index;
}




void Tree_AllocArray(struct tnode_array **new_tree_array, int length)  {

    (*new_tree_array) = malloc(sizeof(struct tnode_array));
    struct tnode_array *tree_array = *new_tree_array;

    tree_array->numnodes = length;
    make_vector(tree_array->ibeg, length);
    make_vector(tree_array->iend, length);
    make_vector(tree_array->numpar, length);
    make_vector(tree_array->x_mid, length);
    make_vector(tree_array->y_mid, length);
    make_vector(tree_array->z_mid, length);
    make_vector(tree_array->x_min, length);
    make_vector(tree_array->y_min, length);
    make_vector(tree_array->z_min, length);
    make_vector(tree_array->x_max, length);
    make_vector(tree_array->y_max, length);
    make_vector(tree_array->z_max, length);
    make_vector(tree_array->level, length);
    make_vector(tree_array->cluster_ind, length);
    make_vector(tree_array->radius, length);
    make_vector(tree_array->num_children, length);
    make_vector(tree_array->children, 8*length);
    
    return;
}   /* END of function allocate_tree_array */




void Tree_FreeArray(struct tnode_array *tree_array)  {

    free_vector(tree_array->ibeg);
    free_vector(tree_array->iend);
    free_vector(tree_array->numpar);
    free_vector(tree_array->x_mid);
    free_vector(tree_array->y_mid);
    free_vector(tree_array->z_mid);
    free_vector(tree_array->x_min);
    free_vector(tree_array->y_min);
    free_vector(tree_array->z_min);
    free_vector(tree_array->x_max);
    free_vector(tree_array->y_max);
    free_vector(tree_array->z_max);
    free_vector(tree_array->level);
    free_vector(tree_array->cluster_ind);
    free_vector(tree_array->radius);
    free_vector(tree_array->num_children);
    free_vector(tree_array->children);

    free(tree_array);
    tree_array = NULL;

    return;
}   /* END of function allocate_tree_array */




void Tree_ReallocArray(struct tnode_array *tree_array, int newlength)  {

    tree_array->numnodes = newlength;
    realloc_vector(tree_array->ibeg, newlength);
    realloc_vector(tree_array->iend, newlength);
    realloc_vector(tree_array->numpar, newlength);
    realloc_vector(tree_array->x_mid, newlength);
    realloc_vector(tree_array->y_mid, newlength);
    realloc_vector(tree_array->z_mid, newlength);
    realloc_vector(tree_array->x_min, newlength);
    realloc_vector(tree_array->y_min, newlength);
    realloc_vector(tree_array->z_min, newlength);
    realloc_vector(tree_array->x_max, newlength);
    realloc_vector(tree_array->y_max, newlength);
    realloc_vector(tree_array->z_max, newlength);
    realloc_vector(tree_array->level, newlength);
    realloc_vector(tree_array->cluster_ind, newlength);
    realloc_vector(tree_array->radius, newlength);

    realloc_vector(tree_array->num_children, newlength);
    realloc_vector(tree_array->children, 8*newlength);

    return;

}   /* END of function allocate_tree_array */




void Tree_CreateArray(struct tnode *p, struct tnode_array *tree_array)
{
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

    for (int i = 0; i < p->num_children; i++) {
        tree_array->children[8*p->node_index+i] = (p->child[i])->node_index;
        Tree_CreateArray(p->child[i], tree_array);
    }
    
    return;
    
} /* END of function create_tree_n0 */




void Tree_Free(struct tnode *p)
{
    remove_node(p);
    free(p);

    return;

} /* END function cleanup */




/*********************************/
/******* LOCAL FUNCTIONS *********/
/*********************************/

static void remove_node(struct tnode *p)
{

    if (p->num_children > 0) {
        for (int i = 0; i < p->num_children; i++) {
            remove_node(p->child[i]);
            free(p->child[i]);
        }
    }

    return;

} /* END function remove_node */
