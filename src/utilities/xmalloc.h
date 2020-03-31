/*
 * C header file for xmalloc.c memory allocation routine
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

#ifndef H_XMALLOC_H
#define H_XMALLOC_H
#include <stdlib.h>
void *malloc_or_exit(size_t nbytes, const char *file, int line);
#define xmalloc(nbytes)  malloc_or_exit((nbytes), __FILE__, __LINE__)
#endif /*H_XMALLOC_H*/
