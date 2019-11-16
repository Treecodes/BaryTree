#ifndef H_DIRECTDRIVER_H
#define H_DIRECTDRIVER_H

#include "struct_particles.h"

void directdriver(struct particles *sources, struct particles *targets,
                  char *kernelName, double kappa, char *singularityHandling,
                  char *approximationName, double *pointwisePotential,
                  double *time_direct);

#endif
