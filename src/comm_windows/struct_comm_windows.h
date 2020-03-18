#ifndef H_STRUCT_COMM_WINDOWS_H
#define H_STRUCT_COMM_WINDOWS_H

#include <mpi.h>


struct CommWindows
{
    MPI_Win win_clusters_x, win_clusters_y, win_clusters_z, win_clusters_q, win_clusters_w;
    MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
};


#endif /* H_STRUCT_COMM_WINDOWS_H */
