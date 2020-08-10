#ifndef H_COMM_WINDOWS_FUNCTIONS_H
#define H_COMM_WINDOWS_FUNCTIONS_H

#include "../clusters/struct_clusters.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"
#include "../comm_types/struct_comm_types.h"

#include "struct_comm_windows.h"


void CommWindows_Create(struct CommWindows **comm_windows_addr,
                        struct Clusters *clusters, struct Particles *sources, struct RunParams *run_params);

void CommWindows_Free(struct CommWindows **comm_windows_addr, struct RunParams *run_params);

void CommWindows_Lock(struct CommWindows *comm_windows, int get_from, struct RunParams *run_params);

void CommWindows_Unlock(struct CommWindows *comm_windows, int get_from, struct RunParams *run_params);

void CommWindows_GetData(struct Clusters *let_clusters, struct Particles *let_sources,
                         struct CommTypes *comm_types, struct CommWindows *comm_windows,
                         int get_from, struct RunParams *run_params);


#endif /* H_COMM_WINDOWS_FUNCTIONS_H */
