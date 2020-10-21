#ifndef H_STRUCT_INTERACTION_LISTS_H
#define H_STRUCT_INTERACTION_LISTS_H


struct InteractionLists
{
    int *num_pp;
    int *num_cc;
    int *num_pc;
    int *num_cp;
    
    int **pp_interactions;
    int **cc_interactions;
    int **pc_interactions;
    int **cp_interactions;
};


#endif /* H_STRUCT_INTERACTION_LISTS_H */
