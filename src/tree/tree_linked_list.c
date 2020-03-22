#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../utilities/tools.h"
#include "../particles/struct_particles.h"

#include "struct_tree_linked_list_node.h"
#include "tree_linked_list.h"
#include "partition.h"


static void remove_node(struct TreeLinkedListNode *p);


void TreeLinkedList_Targets_Construct(struct TreeLinkedListNode **p, struct Particles *targets,
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

    (*p) = malloc(sizeof(struct TreeLinkedListNode));

    (*numnodes)++;

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

    
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    double t1 = (*p)->x_max - (*p)->x_mid;
    double t2 = (*p)->y_max - (*p)->y_mid;
    double t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    (*p)->num_children = 0;
    for (int i = 0; i < 8; i++)
        (*p)->child[i] = NULL;
    

    if ((*p)->numpar > maxparnode) {
    /*
     * IND array holds indices of the eight new subregions.
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

                struct TreeLinkedListNode **paddress = &((*p)->child[idx]);

                TreeLinkedList_Targets_Construct(paddress,
                               targets, ind[i][0], ind[i][1],
                               maxparnode, lxyzmm, loclev, numnodes, numleaves);
            }
        }

    } else {

        (*numleaves)++;
    }

    return;

} /* end of function Tree_CP_Create */




void TreeLinkedList_Sources_Construct(struct TreeLinkedListNode **p, struct Particles *sources,
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
                        

    (*p) = malloc(sizeof(struct TreeLinkedListNode));


    (*numnodes)++;

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

    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    double t1 = (*p)->x_max - (*p)->x_mid;
    double t2 = (*p)->y_max - (*p)->y_mid;
    double t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    (*p)->num_children = 0;
    for (int i = 0; i < 8; i++)
        (*p)->child[i] = NULL;

    
    if ((*p)->numpar > maxparnode) {
    /*
     * IND array holds indices of the eight new subregions.
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

                struct TreeLinkedListNode **paddress = &((*p)->child[idx]);

                TreeLinkedList_Sources_Construct(paddress,
                               sources, ind[i][0], ind[i][1],
                               maxparnode, lxyzmm, loclev,
                               numnodes, numleaves);
            }
        }
        
    } else {
        (*numleaves)++;
    }

    return;
    
} /* END of function PC_Create_Tree */



int TreeLinkedList_SetIndex(struct TreeLinkedListNode *p, int index)
{
        int current_index = index;
        p->node_index = current_index;

        for (int i = 0; i < p->num_children; i++)
            current_index = TreeLinkedList_SetIndex(p->child[i], current_index + 1);

        return current_index;
}



void TreeLinkedList_Free(struct TreeLinkedListNode **p_addr)
{
    struct TreeLinkedListNode *p = *p_addr;
    
    if (p != NULL) {
        remove_node(p);
        free(p);
    }
    
    p = NULL;

    return;
    
} /* END function Tree_Free */



/*********************************/
/******* LOCAL FUNCTIONS *********/
/*********************************/

static void remove_node(struct TreeLinkedListNode *p)
{

    if (p->num_children > 0) {
        for (int i = 0; i < p->num_children; i++) {
            remove_node(p->child[i]);
            free(p->child[i]);
        }
    }

    return;
    
} /* END function remove_node */
