#ifndef H_STRUCT_COMM_TYPES_H
#define H_STRUCT_COMM_TYPES_H

#include <mpi.h>


struct CommTypes
{
        int *num_remote_approx_array;
        int *previous_let_clusters_length_array;

        int let_clusters_length;
        int let_clusters_num;

        MPI_Datatype *MPI_approx_type;
        MPI_Datatype *MPI_approx_charges_type; 
        MPI_Datatype *MPI_approx_weights_type;

        int *new_sources_length_array;
        int *previous_let_sources_length_array;

        int let_sources_length;

        MPI_Datatype *MPI_direct_type;
};


#endif /* H_STRUCT_COMM_TYPES_H */
