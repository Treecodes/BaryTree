#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

#include "struct_particles.h"
#include "struct_kernel.h"
#include "struct_output.h"

/* declaration of primary treecode driver */

void treedriver(struct particles *sources, struct particles *targets,
                int order, double theta, int maxparnode, int batch_size,
                struct kernel *kernel, char *singularityHandling, char *approximationName,
                int tree_type, struct output *output, double *timetree, double sizeCheckFactor,
                int verbosity);

#endif /* H_TREEDRIVER_H */
