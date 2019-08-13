/*
 * C header file containing macros for vector and array creation
 * using the xmalloc.c memory allocation routine
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


#ifndef H_ARRAY_H
#define H_ARRAY_H
#include "xmalloc.h"

#define make_vector(v,n) ((v) = xmalloc((n) * sizeof *(v)))

#define free_vector(v)  do { free(v); v = NULL; } while (0)

#define realloc_vector(v,n) ((v) = realloc(v, (n) * sizeof *(v)))


#define make_matrix(a, m, n) do {                              \
    size_t make_matrix_loop_counter;                           \
    make_vector(a, (m) + 1);                                   \
    for (make_matrix_loop_counter = 0;                         \
            make_matrix_loop_counter < (size_t)(m);            \
            make_matrix_loop_counter++)                        \
        make_vector((a)[make_matrix_loop_counter], (n));       \
    (a)[m] = NULL;                                             \
} while (0)

#define free_matrix(a) do {                                    \
    if (a != NULL) {                                           \
        size_t make_matrix_loop_counter;                       \
        for (make_matrix_loop_counter = 0;                     \
                (a)[make_matrix_loop_counter] != NULL;         \
                make_matrix_loop_counter++)                    \
            free_vector((a)[make_matrix_loop_counter]);        \
        free_vector(a);                                        \
        a = NULL;                                              \
    }                                                          \
} while (0)

#define make_3array(a, l, m, n) do {                           \
    size_t make_3array_loop_counter;                           \
    make_vector(a, (l) + 1);                                   \
    for (make_3array_loop_counter = 0;                         \
            make_3array_loop_counter < (size_t)(l);            \
            make_3array_loop_counter++)                        \
        make_matrix((a)[make_3array_loop_counter], (m), (n));  \
    (a)[l] = NULL;                                             \
} while (0)

#define free_3array(a) do {                                    \
    if (a != NULL) {                                           \
        size_t make_3array_loop_counter;                       \
        for (make_3array_loop_counter = 0;                     \
                (a)[make_3array_loop_counter] != NULL;         \
                make_3array_loop_counter++)                    \
            free_matrix((a)[make_3array_loop_counter]);        \
        free_vector(a);                                        \
        a = NULL;                                              \
    }                                                          \
} while (0)

#define print_vector(fmt, v, n) do {                           \
    size_t print_vector_loop_counter;                          \
    for (print_vector_loop_counter = 0;                        \
            print_vector_loop_counter < (size_t)(n);           \
            print_vector_loop_counter++)                       \
        printf(fmt, (v)[print_vector_loop_counter]);           \
    putchar('\n');                                             \
} while (0)

#define print_matrix(fmt, a, m, n) do {                        \
    size_t print_matrix_loop_counter;                          \
    for (print_matrix_loop_counter = 0;                        \
            print_matrix_loop_counter < (size_t)(m);           \
            print_matrix_loop_counter++)                       \
    print_vector(fmt, (a)[print_matrix_loop_counter], (n));    \
} while (0)



#define reallocate_sources(struct particles *sources, newlength) do {

	sources->num = newlength;
	realloc_vector(sources->x, newlength);
	realloc_vector(sources->y, newlength);
	realloc_vector(sources->z, newlength);
	realloc_vector(sources->q, newlength);
	realloc_vector(sources->w, newlength);

} while(0) /* END of function reallocate_sources */

#define allocate_sources(struct particles *sources, length) do {

	sources->num = length;
	make_vector(sources->x, length);
	make_vector(sources->y, length);
	make_vector(sources->z, length);
	make_vector(sources->q, length);
	make_vector(sources->w, length);

    return;
} while(0) /* END of function allocate_sources */


#define allocate_cluster(struct particles *clusters, length) do {

    make_vector(clusters->x, length);
    make_vector(clusters->y, length);
    make_vector(clusters->z, length);
    make_vector(clusters->q, length);
    make_vector(clusters->w, length);  // will be used in singularity subtraction
    clusters->num=length;

    return;
} while(0) /* END of function allocate_cluster */


#define reallocate_cluster(struct particles *clusters, newlength) do {

	realloc_vector(clusters->x, newlength);
	realloc_vector(clusters->y, newlength);
	realloc_vector(clusters->z, newlength);
	realloc_vector(clusters->q, newlength);
	realloc_vector(clusters->w, newlength);  // will be used in singularity subtraction
    clusters->num=newlength;

    return;
} while(0) /* END of function reallocate_cluster */

#define allocate_tree_array(tnode_array *let_tree_array, length) do {

	let_tree_array->numnodes = length;
	make_vector(let_tree_array->ibeg, length);
	make_vector(let_tree_array->iend, length);
	make_vector(let_tree_array->numpar, length);
	make_vector(let_tree_array->x_mid, length);
	make_vector(let_tree_array->y_mid, length);
	make_vector(let_tree_array->z_mid, length);
	make_vector(let_tree_array->x_min, length);
	make_vector(let_tree_array->y_min, length);
	make_vector(let_tree_array->z_min, length);
	make_vector(let_tree_array->x_max, length);
	make_vector(let_tree_array->y_max, length);
	make_vector(let_tree_array->z_max, length);
	make_vector(let_tree_array->level, length);
	make_vector(let_tree_array->cluster_ind, length);
	make_vector(let_tree_array->radius, length);

    return;
} while(0) /* END of function allocate_tree_array */

#define free_tree_array(tnode_array *let_tree_array, length) do {

	free_vector(let_tree_array->ibeg);
	free_vector(let_tree_array->iend, length);
	free_vector(let_tree_array->numpar, length);
	free_vector(let_tree_array->x_mid, length);
	free_vector(let_tree_array->y_mid, length);
	free_vector(let_tree_array->z_mid, length);
	free_vector(let_tree_array->x_min, length);
	free_vector(let_tree_array->y_min, length);
	free_vector(let_tree_array->z_min, length);
	free_vector(let_tree_array->x_max, length);
	free_vector(let_tree_array->y_max, length);
	free_vector(let_tree_array->z_max, length);
	free_vector(let_tree_array->level, length);
	free_vector(let_tree_array->cluster_ind, length);
	free_vector(let_tree_array->radius, length);
	let_tree_array=NULL

    return;
} while(0) /* END of function allocate_tree_array */


#define reallocate_tree_array(tnode_array *let_tree_array, newlength) do {

	let_tree_array->numnodes = newlength;
	realloc_vector(let_tree_array->ibeg, newlength);
	realloc_vector(let_tree_array->iend, newlength);
	realloc_vector(let_tree_array->numpar, newlength);
	realloc_vector(let_tree_array->x_mid, newlength);
	realloc_vector(let_tree_array->y_mid, newlength);
	realloc_vector(let_tree_array->z_mid, newlength);
	realloc_vector(let_tree_array->x_min, newlength);
	realloc_vector(let_tree_array->y_min, newlength);
	realloc_vector(let_tree_array->z_min, newlength);
	realloc_vector(let_tree_array->x_max, newlength);
	realloc_vector(let_tree_array->y_max, newlength);
	realloc_vector(let_tree_array->z_max, newlength);
	realloc_vector(let_tree_array->level, newlength);
	realloc_vector(let_tree_array->cluster_ind, newlength);
	realloc_vector(let_tree_array->radius, newlength);

} while(0) /* END of function allocate_tree_array */



#endif /*H_ARRAY_H*/
