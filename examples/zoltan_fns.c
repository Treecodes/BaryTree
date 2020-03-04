#include <stdio.h>
#include <stdlib.h>
#include <zoltan.h>

#include "zoltan_fns.h"


int ztn_get_number_of_objects(void *data, int *ierr)
{
    MESH_DATA *mesh = (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return mesh->numMyPoints;
}


void ztn_get_object_list(void *data, int sizeGID, int sizeLID,
                     ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                     int wgt_dim, float *obj_wgts, int *ierr)
{
    int i;
    MESH_DATA *mesh = (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our objects, but no weights.
   * Zoltan will assume equally weighted objects.
   */

    for (i = 0; i < mesh->numMyPoints; i++) {
        globalID[i] = mesh->myGlobalIDs[i];
        localID[i] = i;
        obj_wgts[i] = mesh->b[i];
    }
}


int ztn_get_num_geometry(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    return 3;
}


void ztn_get_geometry_list(void *data, int sizeGID, int sizeLID, int num_obj,
                       ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                       int num_dim, double *geom_vec, int *ierr)
{
    int i;

    MESH_DATA *mesh = (MESH_DATA *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)) {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    for (i = 0;  i < num_obj ; i++){
        geom_vec[3*i] = (double)mesh->x[i];
        geom_vec[3*i + 1] = (double)mesh->y[i];
        geom_vec[3*i + 2] = (double)mesh->z[i];
    } 

    return;
}


void ztn_pack(void *data, int num_gid_entries, int num_lid_entries,
         ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
         int dest, int size, char *buf, int *ierr) {

    SINGLE_MESH_DATA *mesh_single = (SINGLE_MESH_DATA *)buf;
    MESH_DATA *mesh = (MESH_DATA *)data;

    mesh_single->x = mesh->x[(*local_id)];
    mesh_single->y = mesh->y[(*local_id)];
    mesh_single->z = mesh->z[(*local_id)];
    mesh_single->q = mesh->q[(*local_id)];
    mesh_single->w = mesh->w[(*local_id)];
    mesh_single->b = mesh->b[(*local_id)];
    mesh_single->myGlobalID = mesh->myGlobalIDs[(*local_id)];

    mesh->myGlobalIDs[(*local_id)] = (ZOLTAN_ID_TYPE)(-1); // Mark local particle as exported

    return;
}


void ztn_unpack(void *data, int num_gid_entries,
         ZOLTAN_ID_PTR global_id,
         int size, char *buf, int *ierr) {

    SINGLE_MESH_DATA *mesh_single = (SINGLE_MESH_DATA *)buf;
    MESH_DATA *mesh = (MESH_DATA *)data;

    mesh->numMyPoints += 1;

    mesh->myGlobalIDs = (ZOLTAN_ID_TYPE *)realloc(mesh->myGlobalIDs,
                        sizeof(ZOLTAN_ID_TYPE) * mesh->numMyPoints);
    mesh->x = (double *)realloc(mesh->x, sizeof(double) * mesh->numMyPoints);
    mesh->y = (double *)realloc(mesh->y, sizeof(double) * mesh->numMyPoints);
    mesh->z = (double *)realloc(mesh->z, sizeof(double) * mesh->numMyPoints);
    mesh->q = (double *)realloc(mesh->q, sizeof(double) * mesh->numMyPoints);
    mesh->w = (double *)realloc(mesh->w, sizeof(double) * mesh->numMyPoints);
    mesh->b = (double *)realloc(mesh->b, sizeof(double) * mesh->numMyPoints);

    mesh->x[mesh->numMyPoints-1] = mesh_single->x;
    mesh->y[mesh->numMyPoints-1] = mesh_single->y;
    mesh->z[mesh->numMyPoints-1] = mesh_single->z;
    mesh->q[mesh->numMyPoints-1] = mesh_single->q;
    mesh->w[mesh->numMyPoints-1] = mesh_single->w;
    mesh->b[mesh->numMyPoints-1] = mesh_single->b;
    mesh->myGlobalIDs[mesh->numMyPoints-1] = mesh_single->myGlobalID;

    return;
}


int ztn_obj_size(void *data, int num_gid_entries, int num_lid_entries, 
         ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
    return sizeof(SINGLE_MESH_DATA);
}
