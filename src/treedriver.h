#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

#include "particles_struct.h"

/* declaration of primary treecode driver */

void treedriver(struct particles *sources, struct particles *targets,
                int order, double theta, int maxparnode, int batch_size,
                char *kernel, double kappa, char *singularityHandling, char *approximationName,
                int tree_type, double *tEn, double *tpeng, double *timetree);

#endif /* H_TREEDRIVER_H */
