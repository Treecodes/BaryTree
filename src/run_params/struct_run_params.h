#ifndef H_RUN_PARAMS_H
#define H_RUN_PARAMS_H

#include "../utilities/enums.h"


struct RunParams
{
    KERNEL kernel;
    int num_kernel_params;
    double *kernel_params;

    APPROXIMATION approximation;
    SINGULARITY singularity;
    COMPUTE_TYPE compute_type;

    BOUNDARY_CONDITION boundary_type_x;
    BOUNDARY_CONDITION boundary_type_y;
    BOUNDARY_CONDITION boundary_type_z;

    double boundary_length_x;
    double boundary_length_y;
    double boundary_length_z;

    double theta;
    double size_check_factor;

    int interp_order;
    int interp_pts_per_cluster;
    int interp_charges_per_cluster;
    int interp_weights_per_cluster;

    int max_per_source_leaf;
    int max_per_target_leaf;

    int verbosity;
};


#endif /* H_RUN_PARAMS_H */
