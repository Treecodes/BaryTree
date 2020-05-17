#ifndef H_STRUCT_INTERACTION_LISTS_H
#define H_STRUCT_INTERACTION_LISTS_H


struct InteractionLists
{
    int **approx_interactions;
    int **direct_interactions;

    int *num_approx;
    int *num_direct;
    
    int **cc_source_approx_interactions;
    int **cc_target_approx_interactions;
    
    int *num_cc_source_approx;
    int *num_cc_target_approx;
};


#endif /* H_STRUCT_INTERACTION_LISTS_H */
