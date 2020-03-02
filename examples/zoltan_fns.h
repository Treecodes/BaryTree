#ifndef H_ZOLTAN_SUPPORT_FUNCTIONS_H
#define H_ZOLTAN_SUPPORT_FUNCTIONS_H

#include <zoltan.h>

typedef struct{
    int numGlobalPoints;
    int numMyPoints;
    ZOLTAN_ID_PTR myGlobalIDs;
    double *x;
    double *y;
    double *z;
    double *q;
    double *w;
    double *b;
} MESH_DATA;

typedef struct{
    ZOLTAN_ID_TYPE myGlobalID;
    double x;
    double y;
    double z;
    double q;
    double w;
    double b;
} SINGLE_MESH_DATA;


int ztn_get_number_of_objects(void *data, int *ierr);

void ztn_get_object_list(void *data, int sizeGID, int sizeLID,
     ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
     int wgt_dim, float *obj_wgts, int *ierr);

int ztn_get_num_geometry(void *data, int *ierr);

void ztn_get_geometry_list(void *data, int sizeGID, int sizeLID,
     int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
     int num_dim, double *geom_vec, int *ierr);

void ztn_pack(void *data, int num_gid_entries, int num_lid_entries,
     ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
     int dest, int size, char *buf, int *ierr);

void ztn_unpack(void *data, int num_gid_entries,
     ZOLTAN_ID_PTR global_id,
     int size, char *buf, int *ierr);

int ztn_obj_size(void *data, int num_gid_entries, int num_lid_entries,
     ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr);


#endif /* H_ZOLTAN_SUPPORT_FUNCTIONS_H */
