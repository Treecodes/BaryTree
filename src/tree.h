#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "struct_nodes.h"
#include "struct_particles.h"


/* declaration of treecode support functions */

/* used by cluster-particle and particle-cluster */
void Tree_Free(struct tnode *p);

void Tree_Setup(struct particles *particles1, struct particles *particles2,
                int order, double theta, double *xyzminmax);

int Tree_SetIndex(struct tnode *p, int index);


/* used by particle-cluster */
void Tree_PC_Create(struct tnode **p, struct particles *sources,
                    int ibeg, int iend, int maxparnode, double *xyzmm,
                    int level, int *numnodes, int *numleaves);


/* used by cluster-particle */
void Tree_CP_Create(struct tnode **p, struct particles *targets,
                    int ibeg, int iend, int maxparnode, double *xyzmm,
                    int level, int *numnodes, int * numleaves);


/* used for tree arrays */
void Tree_CreateArray(struct tnode *p, struct tnode_array *tree_array);

void Tree_AllocArray(struct tnode_array **new_tree_array, int length);

void Tree_ReallocArray(struct tnode_array *tree_array, int newlength);

void Tree_FreeArray(struct tnode_array *tree_array);

#endif /* H_TREEFUNCTIONS_H */
