/*
 * C code for a safe memory allocation routine
 *
 * This C code written by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on the work of Rouben Rostamian, presented in
 * "Programming Projects in C for Students of Engineering,
 *  Science, and Mathematics"
 *
 * Last modified by Leighton Wilson, 06/23/2016
 */

#include <stdio.h>
#include "xmalloc.h"

void *malloc_or_exit(size_t nbytes, const char *file, int line)
{
    void *x;
    if ((x = malloc(nbytes)) == NULL || nbytes == 0) 
    {
        fprintf(stderr, "%s: line %d: malloc() of %zu bytes failed\n",
                file, line, nbytes);
        exit(EXIT_FAILURE);
    }
    else
        return x;
}
