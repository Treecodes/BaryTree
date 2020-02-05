#ifndef H_DIRECTDRIVER_H
#define H_DIRECTDRIVER_H

#include "struct_particles.h"
#include "struct_kernel.h"

void directdriver(struct particles *sources, struct particles *targets,
                  struct kernel *kernelName, char *singularityHandling,
                  char *approximationName, double *pointwisePotential,
                  double *time_direct);

#endif
