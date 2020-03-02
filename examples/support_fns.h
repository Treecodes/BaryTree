#ifndef H_SUPPORT_FUNCTIONS_H
#define H_SUPPORT_FUNCTIONS_H

#include <stdlib.h>

#include "../src/struct_run_params.h"

void Parse_Params(FILE *fp, struct RunParams **run_params, int *N, int *M, int *run_direct, int *slice);

#endif /* H_SUPPORT_FUNCTIONS_H */
