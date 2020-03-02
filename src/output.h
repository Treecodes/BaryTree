#ifndef H_OUTPUTFUNCTIONS_H
#define H_OUTPUTFUNCTIONS_H

#include "struct_output.h"


void Output_Alloc(struct output *output, int numberOfParticles, int forces);

void Output_Free(struct output *output);


#endif /* H_OUTPUTFUNCTIONS_H */
