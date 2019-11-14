#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "struct_nodes.h"
#include "struct_particles.h"


/* declaration of treecode support functions */

/* used by cluster-particle and particle-cluster */
void Tree_Free(struct tnode *p);

void Tree_Setup(struct particles *particles1, struct particles *particles2,
                int order, double theta, double *xyzminmax);


/* used by particle-cluster */
void Tree_PC_Create(struct tnode **p, struct particles *sources,
                    int ibeg, int iend, int maxparnode, double *xyzmm,
                    int level, int *numnodes, int *numleaves);

int Tree_PC_SetIndex(struct tnode *p, int index);

void Tree_PC_CreateArray(struct tnode *p, struct tnode_array *tree_array);


/* used by cluster-particle */
void Tree_CP_Create(struct tnode **p, struct particles *targets,
                    int ibeg, int iend, int maxparnode, double *xyzmm,
                    int level, int *numnodes, int * numleaves);

#endif /* H_TREEFUNCTIONS_H */
