#ifndef H_STRUCT_TREE_LINKED_LIST_NODE_H
#define H_STRUCT_TREE_LINKED_LIST_NODE_H

struct TreeLinkedListNode
{
    int numpar, ibeg, iend;
    double x_min, y_min, z_min;
    double x_max, y_max, z_max;
    double x_mid, y_mid, z_mid;
    double radius, sqradius, aspect;
    int level, num_children, exist_ms;
    double *ms, *ms2;
    struct TreeLinkedListNode *child[8];

    int xdim, ydim, zdim;
    int xlind, ylind, zlind;
    int xhind, yhind, zhind;

    int node_index;

    double *tx, *ty, *tz;
};

#endif /* H_STRUCT_TREE_LINKED_LIST_NODE_H */
