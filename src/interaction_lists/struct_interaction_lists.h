#ifndef H_STRUCT_INTERACTION_LISTS_H
#define H_STRUCT_INTERACTION_LISTS_H


struct InteractionLists
{
    int **approx_interactions;
    int **direct_interactions;

    int *num_approx;
    int *num_direct;
};


#endif /* H_STRUCT_INTERACTION_LISTS_H */
