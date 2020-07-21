#ifndef H_STRUCT_TREE_LINKED_LIST_NODE_H
#define H_STRUCT_TREE_LINKED_LIST_NODE_H

struct TreeLinkedListNode
{
    int numpar, ibeg, iend;
    
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    
    double radius, aspect;
    
    int num_children;
    struct TreeLinkedListNode *child[8];
    struct TreeLinkedListNode *parent;

    int node_index;

    int level;
};

#endif /* H_STRUCT_TREE_LINKED_LIST_NODE_H */
