#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include "array.h"
#include "tools.h"
#include "globvars.h"

#include "struct_nodes.h"
#include "struct_particles.h"
#include "struct_clusters.h"

#include "interaction_lists.h"
#include "interaction_compute.h"
#include "tree.h"
#include "batches.h"
#include "clusters.h"
#include "particles.h"

#include "treedriver.h"


/* definition of primary treecode driver */

void treedriver(struct particles *sources, struct particles *targets,
                int interpolationOrder, double theta, int maxparnode, int batch_size,
                char *kernelName, double kappa, char *singularityHandling,
                char *approximationName,
                int tree_type, double *tEn, double *time_tree)
{

    int verbosity = 1;

    int rank, numProcs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    if (verbosity > 0) printf("Set rank %i and numProcs %i.\n", rank, numProcs);

    /* date and time */
    double time1;
    
    int *tree_inter_list, *local_tree_inter_list;
    int *direct_inter_list, *local_direct_inter_list;

    struct tnode *troot = NULL;
    struct tnode_array *tree_array = NULL;
    double xyzminmax[6];
    int numnodes = 0, numleaves = 0;

    struct tnode_array *batches = NULL;
    double batch_lim[6];

    struct clusters *clusters = NULL;

    int max_batch_approx = 0;
    int max_batch_direct = 0;

    int totalNumberDirect=0;
    int totalNumberApprox=0;
    int totalNumberInteractions=0;
    int cumulativeNumberInteractions=0;
    int maxNumberInteractions=0;
    int minNumberInteractions=0;

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    if (tree_type == 0) {

        fprintf(stderr, "ERROR: Cluster-particle treecode currently disabled.\n");
        exit(1);

        //setup(targets, order, theta, xyzminmax);
        //cp_create_tree_n0(&troot, targets, 1, targets->num,
        //                  maxparnode, xyzminmax, level);
        //setup_batch(&batches, batch_lim, targets, batch_size);
        //create_source_batch(batches, sources, 1, sources->num,
        //                    batch_size, batch_lim);
        
    } else if (tree_type == 1) {

        time1 = MPI_Wtime();

        Tree_Setup(sources, targets, interpolationOrder, theta, xyzminmax);
        Tree_PC_Create(&troot, sources, 1, sources->num,
                       maxparnode, xyzminmax, 0, &numnodes, &numleaves);
        Tree_SetIndex(troot, 0);
        Tree_AllocArray(&tree_array, numnodes);
        Tree_CreateArray(troot, tree_array);

        time_tree[0] = MPI_Wtime() - time1; //time_maketreearray
        

        time1 = MPI_Wtime();

        Batches_Alloc(&batches, batch_lim, targets, batch_size);
        Batches_CreateTargetBatches(batches, targets, 1, targets->num, batch_size, batch_lim);

        time_tree[1] = MPI_Wtime() - time1; //time_createbatch
        

        time1 = MPI_Wtime();

        Clusters_PC_Setup(&clusters, sources, interpolationOrder, tree_array,
                          approximationName, singularityHandling);

        time_tree[2] = MPI_Wtime() - time1; //time_fillclusters
    }
    

    if (verbosity > 0) {
        printf("Tree information: \n\n");

        printf("                      numpar: %d\n", troot->numpar);
        printf("                       x_mid: %e\n", troot->x_mid);
        printf("                       y_mid: %e\n", troot->y_mid);
        printf("                       z_mid: %e\n\n", troot->z_mid);
        printf("                      radius: %f\n\n", troot->radius);
        printf("                       x_len: %e\n", troot->x_max - troot->x_min);
        printf("                       y_len: %e\n", troot->y_max - troot->y_min);
        printf("                       z_len: %e\n\n", troot->z_max - troot->z_min);
        printf("                      torder: %d\n", interpolationOrder);
        printf("                       theta: %f\n", theta);
        printf("                  maxparnode: %d\n", maxparnode);
        printf("            number of leaves: %d\n", numleaves);
        printf("             number of nodes: %d\n", numnodes);
        printf("           target batch size: %d\n", batch_size);
        printf("           number of batches: %d\n\n", batches->numnodes);
    }



    if (tree_type == 0) {
        fprintf(stderr, "ERROR: Cluster-particle treecode is currently not disabled.\n");
        exit(1);
//        if (pot_type == 0) {
//            cp_treecode(troot, batches, sources, targets,
//                        tpeng, tEn, &timetree[1]);
//        }
//        } else if (pot_type == 1) {
//            cp_treecode_yuk(troot, batches, sources, targets, kappa,
//                            tpeng, tEn, &timetree[1]);
//        }
    } else if (tree_type == 1) {

        MPI_Barrier(MPI_COMM_WORLD);
  
        time1 = MPI_Wtime();

        int numNodesOnProc[numProcs];
        MPI_Allgather(&numnodes, 1, MPI_INT, numNodesOnProc, 1, MPI_INT, MPI_COMM_WORLD);

        int pointsPerCluster = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
        int chargesPerCluster = pointsPerCluster;
        int weightsPerCluster = pointsPerCluster;

        if (strcmp(approximationName, "hermite") == 0)
            chargesPerCluster *= 8;

        if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0))
            weightsPerCluster *= 8;

        struct tnode_array *let_tree_array = NULL;
        int let_tree_array_length = 0;

        struct clusters *let_clusters = NULL;
        let_clusters = malloc(sizeof(struct clusters));
        int let_clusters_length = 0; // previously let_clusters included the local.  Now it should not

        struct particles *let_sources = NULL;
        let_sources = malloc(sizeof(struct particles));  // let_sources will hold all source nodes needed for direct interactions
        int let_sources_length = 0;  // previously let_sources included local.  Now it should not


        MPI_Win win_x_mid, win_y_mid, win_z_mid, win_radius, win_numpar, win_ibeg, win_iend, win_level;
        MPI_Win win_clusters_x, win_clusters_y, win_clusters_z, win_clusters_q, win_clusters_w;
        MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
        MPI_Win win_children, win_num_children;
        
        MPI_Win_create(tree_array->x_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_x_mid);
        MPI_Win_create(tree_array->y_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_y_mid);
        MPI_Win_create(tree_array->z_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_z_mid);
        MPI_Win_create(tree_array->radius, numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_radius);
        MPI_Win_create(tree_array->numpar, numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_numpar);
        MPI_Win_create(tree_array->ibeg,   numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_ibeg);
        MPI_Win_create(tree_array->iend,   numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_iend);
        MPI_Win_create(tree_array->level,  numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_level);
        MPI_Win_create(tree_array->num_children,  numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_num_children);
        MPI_Win_create(tree_array->children,    8*numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_children);

        MPI_Win_create(clusters->x, clusters->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_x);
        MPI_Win_create(clusters->y, clusters->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_y);
        MPI_Win_create(clusters->z, clusters->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_z);
        MPI_Win_create(clusters->q, clusters->num_charges*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_q);
        MPI_Win_create(clusters->w, clusters->num_weights*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_w);

        MPI_Win_create(sources->x, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
        MPI_Win_create(sources->y, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
        MPI_Win_create(sources->z, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
        MPI_Win_create(sources->q, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
        MPI_Win_create(sources->w, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);

        // Perform MPI round robin, filling LET with remote data
        int num_remote_approx_array[numProcs], new_sources_length_array[numProcs];
        int previous_let_clusters_length_array[numProcs], previous_let_sources_length_array[numProcs];
        MPI_Datatype approx_type[numProcs], approx_charges_type[numProcs], approx_weights_type[numProcs];
        MPI_Datatype direct_type[numProcs];
        int let_clusters_num = 0;

        int *sizeof_batch_approx, *sizeof_batch_direct;
        int *offset_batch_approx, *offset_batch_direct;

        make_vector(sizeof_batch_approx, batches->numnodes);
        make_vector(sizeof_batch_direct, batches->numnodes);
        make_vector(offset_batch_approx, batches->numnodes);
        make_vector(offset_batch_direct, batches->numnodes);

        int loopsize = batches->numnodes;
        for (int i = 0; i < loopsize; i++) sizeof_batch_approx[i] = 0;
        for (int i = 0; i < loopsize; i++) sizeof_batch_direct[i] = 0;
        for (int i = 0; i < loopsize; i++) offset_batch_approx[i] = 0;
        for (int i = 0; i < loopsize; i++) offset_batch_approx[i] = 0;

        int total_batch_direct = 0;
        int total_batch_approx = 0;
         

        for (int procID = 1; procID < numProcs; ++procID) {

            int getFrom = (numProcs+rank-procID) % numProcs;

            // Allocate remote_tree_array
            struct tnode_array *remote_tree_array = NULL;
            Tree_AllocArray(&remote_tree_array, numNodesOnProc[getFrom]);

            // Get remote_tree_array
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_x_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_y_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_z_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_radius);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_numpar);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_ibeg);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_iend);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_level);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_children);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_num_children);

            
            MPI_Get(remote_tree_array->x_mid, numNodesOnProc[getFrom], MPI_DOUBLE,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_DOUBLE, win_x_mid);
            MPI_Get(remote_tree_array->y_mid, numNodesOnProc[getFrom], MPI_DOUBLE,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_DOUBLE, win_y_mid);
            MPI_Get(remote_tree_array->z_mid, numNodesOnProc[getFrom], MPI_DOUBLE,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_DOUBLE, win_z_mid);
            MPI_Get(remote_tree_array->radius, numNodesOnProc[getFrom], MPI_DOUBLE,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_DOUBLE, win_radius);
            MPI_Get(remote_tree_array->numpar, numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_INT, win_numpar);
            MPI_Get(remote_tree_array->ibeg, numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_INT, win_ibeg);
            MPI_Get(remote_tree_array->iend, numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_INT, win_iend);
            MPI_Get(remote_tree_array->level, numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_INT, win_level);

            MPI_Get(remote_tree_array->children, 8*numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, 8*numNodesOnProc[getFrom], MPI_INT, win_children);
            MPI_Get(remote_tree_array->num_children, numNodesOnProc[getFrom], MPI_INT,
                    getFrom, 0, numNodesOnProc[getFrom], MPI_INT, win_num_children);
            
            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Win_unlock(getFrom, win_x_mid);
            MPI_Win_unlock(getFrom, win_y_mid);
            MPI_Win_unlock(getFrom, win_z_mid);
            MPI_Win_unlock(getFrom, win_radius);
            MPI_Win_unlock(getFrom, win_numpar);
            MPI_Win_unlock(getFrom, win_ibeg);
            MPI_Win_unlock(getFrom, win_iend);
            MPI_Win_unlock(getFrom, win_level);
            MPI_Win_unlock(getFrom, win_children);
            MPI_Win_unlock(getFrom, win_num_children);
            
            MPI_Barrier(MPI_COMM_WORLD);



            // Construct masks
            int *approx_list_packed, *approx_list_unpacked, *direct_list, *direct_ibeg_list, *direct_length_list;
            make_vector(approx_list_packed, numNodesOnProc[getFrom]);
            make_vector(approx_list_unpacked, numNodesOnProc[getFrom]);
            make_vector(direct_list, numNodesOnProc[getFrom]);
            make_vector(direct_ibeg_list, numNodesOnProc[getFrom]);
            make_vector(direct_length_list, numNodesOnProc[getFrom]);

            Interaction_MakeListRemote(remote_tree_array, batches, approx_list_unpacked, approx_list_packed,
                                       direct_list, sizeof_batch_approx, sizeof_batch_direct);

            //MPI_Barrier(MPI_COMM_WORLD);

            // Count number of unique clusters adding to LET
            int numberOfUniqueClusters =  0;
            int previousTreeArrayLength = let_tree_array_length;
            for (int i = 0; i < numNodesOnProc[getFrom]; ++i) {
                numberOfUniqueClusters++;
                let_tree_array_length++;
            }
            //MPI_Barrier(MPI_COMM_WORLD);
            

            if (procID == 1) {
                Tree_AllocArray(&let_tree_array, let_tree_array_length);
            } else {
                Tree_ReallocArray(let_tree_array, let_tree_array_length);
            }
            //MPI_Barrier(MPI_COMM_WORLD);

            int numberOfRemoteApprox = 0;
            int previous_let_clusters_length = let_clusters_length;

            int numberOfRemoteDirect = 0;
            int previous_let_sources_length = let_sources_length;



            // Fill in LET tree array from Remote tree array.
            int appendCounter = 0;
            for (int i = 0; i < numNodesOnProc[getFrom]; ++i) {

                let_tree_array->x_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->x_mid[i];
                let_tree_array->y_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->y_mid[i];
                let_tree_array->z_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->z_mid[i];
                let_tree_array->radius[previousTreeArrayLength + appendCounter] = remote_tree_array->radius[i];
                let_tree_array->numpar[previousTreeArrayLength + appendCounter] = remote_tree_array->numpar[i];
                let_tree_array->level[previousTreeArrayLength + appendCounter] = remote_tree_array->level[i];
                    
                if (approx_list_unpacked[i] != -1) {
                    let_tree_array->cluster_ind[previousTreeArrayLength + appendCounter] = let_clusters_num;
                    let_clusters_length += pointsPerCluster;
                    let_clusters_num++;
                    numberOfRemoteApprox++;
                }
                    
                if (direct_list[i] != -1) {
                        
                    // Set the beginning and ending particle indices for the associated nodes in the local sources list
                    let_tree_array->ibeg[previousTreeArrayLength + appendCounter] = let_sources_length + 1;  // These are one-index based!!!
                    let_tree_array->iend[previousTreeArrayLength + appendCounter] = let_sources_length + remote_tree_array->numpar[i];
                    let_sources_length += remote_tree_array->numpar[i];
                        
                    // Determine displacements and lengths for getting prticles from remote sources list
                    direct_ibeg_list[numberOfRemoteDirect] = remote_tree_array->ibeg[i] - 1; // These are zero-index based!!!
                    direct_length_list[numberOfRemoteDirect] = remote_tree_array->numpar[i];
                    numberOfRemoteDirect++;
                }
                
                appendCounter++;
            }

            
            num_remote_approx_array[getFrom] = numberOfRemoteApprox;
            new_sources_length_array[getFrom] = let_sources_length - previous_let_sources_length;
            previous_let_clusters_length_array[getFrom] = previous_let_clusters_length; 
            previous_let_sources_length_array[getFrom] = previous_let_sources_length; 
            
            //MPI_Barrier(MPI_COMM_WORLD);
            
            int *approx_list_displacements, *approx_charges_list_displacements, *approx_weights_list_displacements;
            make_vector(approx_list_displacements, numNodesOnProc[getFrom]);
            make_vector(approx_charges_list_displacements, numNodesOnProc[getFrom]);
            make_vector(approx_weights_list_displacements, numNodesOnProc[getFrom]);

            // Use masks to get remote data
            for (int ii = 0; ii < numberOfRemoteApprox; ++ii) {
                approx_list_displacements[ii] = approx_list_packed[ii] * pointsPerCluster;
                approx_charges_list_displacements[ii] = approx_list_packed[ii] * chargesPerCluster;
                approx_weights_list_displacements[ii] = approx_list_packed[ii] * weightsPerCluster;
            }
            
            MPI_Type_create_indexed_block(numberOfRemoteApprox, pointsPerCluster, approx_list_displacements,
                                          MPI_DOUBLE, &approx_type[getFrom]);
            MPI_Type_commit(&approx_type[getFrom]);

            MPI_Type_create_indexed_block(numberOfRemoteApprox, chargesPerCluster, approx_charges_list_displacements,
                                          MPI_DOUBLE, &approx_charges_type[getFrom]);
            MPI_Type_commit(&approx_charges_type[getFrom]);

            MPI_Type_create_indexed_block(numberOfRemoteApprox, weightsPerCluster, approx_weights_list_displacements,
                                          MPI_DOUBLE, &approx_weights_type[getFrom]);
            MPI_Type_commit(&approx_weights_type[getFrom]);

            MPI_Type_indexed(numberOfRemoteDirect, direct_length_list, direct_ibeg_list, 
                                          MPI_DOUBLE, &direct_type[getFrom]);
            MPI_Type_commit(&direct_type[getFrom]);

            free_vector(approx_list_packed);
            free_vector(approx_list_unpacked);
            free_vector(approx_list_displacements);
            free_vector(approx_charges_list_displacements);
            free_vector(approx_weights_list_displacements);
            free_vector(direct_list);
            free_vector(direct_ibeg_list);
            free_vector(direct_length_list);
            Tree_FreeArray(remote_tree_array);
        } //end loop over numProcs

        if (let_sources_length > 0) Particles_AllocSources(let_sources, let_sources_length);
        if (let_clusters_length > 0) Clusters_Alloc(let_clusters, let_clusters_length,
                                                    approximationName, singularityHandling);
    
        for (int procID = 1; procID < numProcs; ++procID) {

            int getFrom = (numProcs+rank-procID) % numProcs;

            MPI_Barrier(MPI_COMM_WORLD);

            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_clusters_x);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_clusters_y);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_clusters_z);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_clusters_w);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_clusters_q);

            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_x);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_y);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_z);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_q);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_w);


            MPI_Get(&(let_clusters->x[previous_let_clusters_length_array[getFrom]]),
                    num_remote_approx_array[getFrom] * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type[getFrom], win_clusters_x);
            MPI_Get(&(let_clusters->y[previous_let_clusters_length_array[getFrom]]),
                    num_remote_approx_array[getFrom] * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type[getFrom], win_clusters_y);
            MPI_Get(&(let_clusters->z[previous_let_clusters_length_array[getFrom]]),
                    num_remote_approx_array[getFrom] * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type[getFrom], win_clusters_z);

            MPI_Get(&(let_clusters->q[(chargesPerCluster/pointsPerCluster) * previous_let_clusters_length_array[getFrom]]),
                    num_remote_approx_array[getFrom] * chargesPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_charges_type[getFrom], win_clusters_q);
            MPI_Get(&(let_clusters->w[(weightsPerCluster/pointsPerCluster) * previous_let_clusters_length_array[getFrom]]),
                    num_remote_approx_array[getFrom] * weightsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_weights_type[getFrom], win_clusters_w);

            MPI_Get(&(let_sources->x[previous_let_sources_length_array[getFrom]]),
                    new_sources_length_array[getFrom], MPI_DOUBLE,
                    getFrom, 0, 1, direct_type[getFrom], win_sources_x);
            MPI_Get(&(let_sources->y[previous_let_sources_length_array[getFrom]]),
                    new_sources_length_array[getFrom], MPI_DOUBLE,
                    getFrom, 0, 1, direct_type[getFrom], win_sources_y);
            MPI_Get(&(let_sources->z[previous_let_sources_length_array[getFrom]]),
                    new_sources_length_array[getFrom], MPI_DOUBLE,
                    getFrom, 0, 1, direct_type[getFrom], win_sources_z);
            MPI_Get(&(let_sources->q[previous_let_sources_length_array[getFrom]]), 
                    new_sources_length_array[getFrom], MPI_DOUBLE,
                    getFrom, 0, 1, direct_type[getFrom], win_sources_q);
            MPI_Get(&(let_sources->w[previous_let_sources_length_array[getFrom]]),
                    new_sources_length_array[getFrom], MPI_DOUBLE,
                    getFrom, 0, 1, direct_type[getFrom], win_sources_w);

            MPI_Barrier(MPI_COMM_WORLD);
            
            MPI_Win_unlock(getFrom, win_clusters_x);
            MPI_Win_unlock(getFrom, win_clusters_y);
            MPI_Win_unlock(getFrom, win_clusters_z);
            MPI_Win_unlock(getFrom, win_clusters_q);
            MPI_Win_unlock(getFrom, win_clusters_w);
            
            MPI_Win_unlock(getFrom, win_sources_x);
            MPI_Win_unlock(getFrom, win_sources_y);
            MPI_Win_unlock(getFrom, win_sources_z);
            MPI_Win_unlock(getFrom, win_sources_q);
            MPI_Win_unlock(getFrom, win_sources_w);

            MPI_Barrier(MPI_COMM_WORLD);

        } // end loop over numProcs

        time_tree[3] = MPI_Wtime() - time1;


        time1 = MPI_Wtime();

        // Local particles
//MAKE THIS BETTER!
        make_vector(local_tree_inter_list, batches->numnodes * tree_array->numnodes);
        make_vector(local_direct_inter_list, batches->numnodes * tree_array->numnodes);
        Interaction_MakeList(tree_array, batches, local_tree_inter_list, local_direct_inter_list,
                             tree_array->numnodes, tree_array->numnodes);

        time_tree[4] = MPI_Wtime() - time1; //time_constructlet


        if (verbosity>0){
            for (int j = 0; j < batches->numnodes; j++){
                totalNumberApprox += batches->numApprox[j];
                totalNumberDirect += batches->numDirect[j];
            }
        }
        time1 = MPI_Wtime();

        Interaction_PC_Compute(tree_array, batches,
                        local_tree_inter_list, local_direct_inter_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                        tEn, interpolationOrder,
                        sources->num, targets->num, clusters->num,
                        tree_array->numnodes, tree_array->numnodes,
                        kernelName, kappa, singularityHandling,
                        approximationName);

        time_tree[5] = MPI_Wtime() - time1; //time_constructlet


        // Compute interaction lists based on LET
        if (numProcs > 1) {
            time1 = MPI_Wtime();

            max_batch_approx = maxval_int(sizeof_batch_approx, batches->numnodes);
            max_batch_direct = maxval_int(sizeof_batch_direct, batches->numnodes);



            if (max_batch_approx > 0) make_vector(tree_inter_list, batches->numnodes * max_batch_approx);
            if (max_batch_direct > 0) make_vector(direct_inter_list, batches->numnodes * max_batch_direct);
            Interaction_MakeList(let_tree_array, batches, tree_inter_list, direct_inter_list,
                                 max_batch_approx, max_batch_direct);

            // Count number of interactions

            if (verbosity>0){
                for (int j = 0; j < batches->numnodes; j++){
                    totalNumberApprox += batches->numApprox[j];
                    totalNumberDirect += batches->numDirect[j];
                }
            }
            time_tree[6] = MPI_Wtime() - time1; //time_makeglobintlist

            // After filling LET, call interaction_list_treecode
            time1 = MPI_Wtime(); // start timer for tree evaluation

            Interaction_PC_Compute(let_tree_array, batches,
                                   tree_inter_list, direct_inter_list,
                                   let_sources->x, let_sources->y, let_sources->z, let_sources->q, let_sources->w,
                                   targets->x, targets->y, targets->z, targets->q,
                                   let_clusters->x, let_clusters->y, let_clusters->z, let_clusters->q, let_clusters->w,
                                   tEn, interpolationOrder,
                                   let_sources->num, targets->num, let_clusters->num,
                                   max_batch_approx, max_batch_direct,
                                   kernelName, kappa, singularityHandling, approximationName);

            time_tree[7] = MPI_Wtime() - time1;
        }

        // Free let clusters and sources, if they were used.
        if (numProcs > 1) {
            Clusters_Free(let_clusters);
            Particles_FreeSources(let_sources);
            Tree_FreeArray(let_tree_array);
        }


        time1 = MPI_Wtime();
        
        Interaction_SubtractionPotentialCorrection(tEn, targets->q, targets->num,
                                  kernelName, kappa, singularityHandling);

        Particles_ReorderTargetsAndPotential(targets, tEn);

        time_tree[8] = MPI_Wtime() - time1;
    }

    time1 = MPI_Wtime();
    if (verbosity > 0) {
        totalNumberInteractions=totalNumberDirect+totalNumberApprox;
        printf("Interaction information: \n");
        printf("rank %d: number of direct batch-cluster interactions: %d\n", rank, totalNumberApprox);
        printf("rank %d: number of approx batch-cluster interactions: %d\n", rank, totalNumberDirect);
        printf("rank %d:  total number of batch-cluster interactions: %d\n\n", rank, totalNumberInteractions);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    if (verbosity > 0) {
        MPI_Reduce(&totalNumberInteractions,&cumulativeNumberInteractions, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&totalNumberInteractions,&maxNumberInteractions, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&totalNumberInteractions,&minNumberInteractions, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        if (rank==0){
            printf("Cumulative number of interactions across all ranks: %d\n", cumulativeNumberInteractions);
            printf("   Maximum number of interactions across all ranks: %d\n", maxNumberInteractions);
            printf("   Minimum number of interactions across all ranks: %d\n", minNumberInteractions);
            printf("                                             Ratio: %f\n\n", (double)maxNumberInteractions/(double)minNumberInteractions );
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    /***********************************************/
    /***************** Cleanup *********************/
    /***********************************************/

    free_vector(targets->order); // free particle order arrays
    free_vector(sources->order); // free particle order arrays
    Tree_Free(troot); // free tree
    free_vector(local_tree_inter_list); // free interaction lists
    free_vector(local_direct_inter_list); // free interaction lists
    if (numProcs > 1) { // free remote interaction lists, if they were built
        if (max_batch_approx > 0) free_vector(tree_inter_list);
        if (max_batch_direct > 0) free_vector(direct_inter_list);
    }
    Clusters_Free(clusters); // free local clusters
    Tree_FreeArray(tree_array); // free tree array
    Batches_Free(batches); // free target batches

    time_tree[9] = MPI_Wtime() - time1; //time_cleanup
    
//    MPI_Barrier(MPI_COMM_WORLD);

    return;
} /* END function treecode */

