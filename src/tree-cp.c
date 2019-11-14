/*
 *Procedures for Cluster-Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "nodes_struct.h"
#include "particles_struct.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"


double thetasq;
double *tt, *ww;


void setup(struct particles *particles1, struct particles *particles2, int order, double theta,
           double *xyzminmax)
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

    /* find bounds of Cartesian box enclosing the particles */
    xyzminmax[0] = minval(particles1->x, particles1->num);
    xyzminmax[1] = maxval(particles1->x, particles1->num);
    xyzminmax[2] = minval(particles1->y, particles1->num);
    xyzminmax[3] = maxval(particles1->y, particles1->num);
    xyzminmax[4] = minval(particles1->z, particles1->num);
    xyzminmax[5] = maxval(particles1->z, particles1->num);

    /* setting up ordering vectors */
    make_vector(particles1->order, particles1->num);
    make_vector(particles2->order, particles2->num);
    for (int i = 0; i < particles1->num; i++) particles1->order[i] = i+1;
    for (int i = 0; i < particles2->num; i++) particles2->order[i] = i+1;

    return;
    
} /* END of function setup */




void cp_create_tree_n0(struct tnode **p, struct particles *targets,
                       int ibeg, int iend, int maxparnode,
                       double *xyzmm, int level, int *numnodes, int *numleaves)
{
	printf("Entering cp_create_tree_n0.\n");
    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild;
    
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

        cp_partition_8(targets->x, targets->y, targets->z, targets->q, targets->order,
                       xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                cp_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                                  targets, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev, numnodes, numleaves);
            }
        }

    } else {

        (*numleaves)++;
    }
    printf("Exiting cp_create_tree_n0.\n");
    return;

} /* end of function create_tree_n0 */




/*
 * cleanup deallocates allocated global variables and then calls
 * recursive function remove_node to delete the tree
 */
void cleanup(struct tnode *p)
{
    remove_node(p);
    free(p);

    return;

} /* END function cleanup */




void remove_node(struct tnode *p)
{
    /* local variables */
    int i;

//    if (p->exist_ms == 1)
//        free(p->ms);
//    	free(p->ms2);

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++) {
            remove_node(p->child[i]);
            free(p->child[i]);
        }
    }

    return;

} /* END function remove_node */
