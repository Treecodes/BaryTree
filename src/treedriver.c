#include <stdio.h>
#include <mpi.h>
#include <limits.h>
#include <omp.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"
#include "tree.h"
#include "structAllocations.h"
#include "interaction-masks.h"

#include "treedriver.h"


/* definition of primary treecode driver */

void treedriver(struct particles *sources, struct particles *targets,
                int order, double theta, int maxparnode, int batch_size,
                int pot_type, double kappa, int tree_type,
                double *tEn, double *tpeng, double *timetree, int numDevices, int numThreads)
{

	int rank; int numProcs;	int ierr;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    int verbosity=1;
    /* local variables */
    struct tnode *troot = NULL;
    int level;
    double xyzminmax[6];
    
//    int i, j;
    
    /* batch variables */
    struct batch *batches = NULL;
    double batch_lim[6];
    
    /* date and time */
    double time1, time2, time3, timeFillClusters1, timeFillClusters2;

    
    time1 = MPI_Wtime();
    
    level = 0;
    numleaves = 0;
    minlevel = INT_MAX;
    minpars = INT_MAX;
    maxlevel = 0;
    maxpars = 0;
    
    int *tree_inter_list, *tree_inter_list2;
    int *direct_inter_list, *direct_inter_list2;
    struct tnode_array *tree_array = NULL;
    numnodes = 0;
    struct particles *clusters = NULL;
    clusters = malloc(sizeof(struct particles));


    time2 = MPI_Wtime();
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
        if (verbosity>0) printf("Treetype %i: entering setup.\n", tree_type);

        time1 = MPI_Wtime();
        setup(sources, order, theta, xyzminmax);
        time2 = MPI_Wtime();

        if (verbosity>0) printf("Time to setup: %f\n", time2-time1);

        time1 = MPI_Wtime();
        
        #pragma omp parallel
        {
            #pragma omp single
            { 
                pc_create_tree_n0(&troot, sources, 1, sources->num,
                                  maxparnode, xyzminmax, level);
            }
        }

        int final_index = pc_set_tree_index(troot, 0);

        time2 = MPI_Wtime();
        printf("Time to pc_create_tree_n0: %f\n", time2-time1);
        printf("Rank %i got past pc_create_tree_n0.\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);
        time1 = MPI_Wtime();
        tree_array = malloc(sizeof(struct tnode_array));
        tree_array->numnodes = numnodes;
        make_vector(tree_array->ibeg, numnodes);
        make_vector(tree_array->iend, numnodes);
        make_vector(tree_array->numpar, numnodes);
        make_vector(tree_array->x_mid, numnodes);
        make_vector(tree_array->y_mid, numnodes);
        make_vector(tree_array->z_mid, numnodes);
        make_vector(tree_array->x_min, numnodes);
        make_vector(tree_array->y_min, numnodes);
        make_vector(tree_array->z_min, numnodes);
        make_vector(tree_array->x_max, numnodes);
        make_vector(tree_array->y_max, numnodes);
        make_vector(tree_array->z_max, numnodes);
        make_vector(tree_array->level, numnodes);
        make_vector(tree_array->cluster_ind, numnodes);
        make_vector(tree_array->radius, numnodes);
        time2 = MPI_Wtime();
//		printf("Time to make tree_array: %f\n", time2-time1);

        time1 = MPI_Wtime();
        pc_create_tree_array(troot, tree_array);
        time2 = MPI_Wtime();
        if (verbosity>0) printf("Rank %i exiting pc_create_tree_array.\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);
//        printf("Time to pc_create_tree_array: %f\n", time2-time1);

        time1 = MPI_Wtime();
        setup_batch(&batches, batch_lim, targets, batch_size);
        time2 = MPI_Wtime();
        if (verbosity>0) printf("Rank %i exiting setup_batch.\n", rank);
		MPI_Barrier(MPI_COMM_WORLD);
//        printf("Time to setup_batch: %f\n", time2-time1);

        time1 = MPI_Wtime();
        printf("Rank %i:  targets->num %i, batch_size %i, batch_lim %i .\n", rank, targets->num, batch_size, batch_lim);

//        if (rank==1){
//        	for (int i=0;i<targets->num;i++){
//        		printf("%f\n",targets->x);
//        	}
//        }
        create_target_batch(batches, targets, 1, targets->num, batch_size, batch_lim);
        time2 = MPI_Wtime();
//        printf("Time to create_target_batch: %f\n", time2-time1);
        if (verbosity>0) printf("Rank %i exiting create_target_batch.\n", rank);
        MPI_Barrier(MPI_COMM_WORLD);

        timeFillClusters1 = MPI_Wtime();
        if         ((pot_type == 0) || (pot_type==1))  {
            fill_in_cluster_data(   clusters, sources, troot, order,
                                  numDevices, numThreads, tree_array);

        } else if  ((pot_type == 2) || (pot_type==3)) {
            fill_in_cluster_data_SS(clusters, sources, troot, order);

        } else if  ((pot_type == 4) || (pot_type==5)) {
            fill_in_cluster_data_hermite(clusters, sources, troot, order);

        } else if  ((pot_type == 6) || (pot_type==7)) {
            fill_in_cluster_data_hermite(clusters, sources, troot, order);
        }
        timeFillClusters2 = MPI_Wtime();
        timeFillClusters1 = timeFillClusters2-timeFillClusters1;
        printf("Time to compute modified weights(s):  %f\n", timeFillClusters1);

    }

//	MPI_Barrier(MPI_COMM_WORLD);
    printf("Rank %i got past fill_in_cluster_data.\n", rank);
	MPI_Barrier(MPI_COMM_WORLD);

    time2 = MPI_Wtime();
    timetree[0] = time2-time1;

    if (verbosity>1) {
        printf("Tree creation (s):  %f\n\n", time2-time1);
        printf("Tree information: \n\n");

        printf("                      numpar: %d\n", troot->numpar);
        printf("                       x_mid: %e\n", troot->x_mid);
        printf("                       y_mid: %e\n", troot->y_mid);
        printf("                       z_mid: %e\n\n", troot->z_mid);
        printf("                      radius: %f\n\n", troot->radius);
        printf("                       x_len: %e\n", troot->x_max - troot->x_min);
        printf("                       y_len: %e\n", troot->y_max - troot->y_min);
        printf("                       z_len: %e\n\n", troot->z_max - troot->z_min);
        printf("                      torder: %d\n", torder);
        printf("                       theta: %f\n", theta);
        printf("                  maxparnode: %d\n", maxparnode);
        printf("               tree maxlevel: %d\n", maxlevel);
        printf("               tree minlevel: %d\n", minlevel);
        printf("                tree maxpars: %d\n", maxpars);
        printf("                tree minpars: %d\n", minpars);
        printf("            number of leaves: %d\n", numleaves);
        printf("             number of nodes: %d\n", numnodes);
        printf("           number of devices: %d\n", numDevices);
        printf("           target batch size: %d\n", batch_size);
        printf("           number of batches: %d\n\n", batches->num);
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




		printf("Got  MPI_Comm_rank.\n");
		MPI_Barrier(MPI_COMM_WORLD);
    	// Get number of nodes for each processor's LET
		int numNodesOnProc[numProcs];
    	MPI_Allgather(&numnodes, 1, MPI_INT,
    			numNodesOnProc, 1, MPI_INT,
				MPI_COMM_WORLD);

    	int pointsPerCluster = (order+1)*(order+1)*(order+1);
		int let_sources_length;
    	int let_clusters_length;
    	// Allocate structures before MPI round robin

    	struct tnode_array *let_tree_array = NULL;
    	int let_tree_array_length=numnodes;
    	let_tree_array = malloc(sizeof(struct tnode_array));
    	allocate_tree_array(let_tree_array,let_tree_array_length); // start by allocating let_tree_array with size numnodes

		struct particles *let_clusters = NULL;
		let_clusters_length=numnodes*pointsPerCluster;
		let_clusters = malloc(sizeof(struct particles)); // let_clusters will hold all cluster data for LET
		allocate_cluster(let_clusters,let_clusters_length);

		struct particles *let_sources = NULL;
		let_sources = malloc(sizeof(struct particles));  // let_sources will hold all source nodes needed for direct interactions
		let_sources_length = troot->numpar;
		allocate_sources(let_sources,let_sources_length);



		// Fill LET with own data



		for (int i=0; i<numnodes; i++){
			let_tree_array->ibeg[i] = tree_array->ibeg[i];
			let_tree_array->iend[i] = tree_array->iend[i];
			let_tree_array->numpar[i] = tree_array->numpar[i];
			let_tree_array->x_mid[i] = tree_array->x_mid[i];
			let_tree_array->y_mid[i] = tree_array->y_mid[i];
			let_tree_array->z_mid[i] = tree_array->z_mid[i];
			let_tree_array->level[i] = tree_array->level[i];
			let_tree_array->cluster_ind[i] = tree_array->cluster_ind[i];
			let_tree_array->radius[i] = tree_array->radius[i];
		}

		for (int i=0; i<numnodes*pointsPerCluster; i++){
			// Fill in clusters
			let_clusters->x[i] = clusters->x[i];
			let_clusters->y[i] = clusters->y[i];
			let_clusters->z[i] = clusters->z[i];
			let_clusters->q[i] = clusters->q[i];
			let_clusters->w[i] = clusters->w[i];
		}
		for (int i=0; i<troot->numpar; i++){
			// Fill in own sources
			let_sources->x[i] = sources->x[i];
			let_sources->y[i] = sources->y[i];
			let_sources->z[i] = sources->z[i];
			let_sources->q[i] = sources->q[i];
			let_sources->w[i] = sources->w[i];
		}
        
        MPI_Win win_x_mid, win_y_mid, win_z_mid, win_radius, win_numpar, win_ibeg, win_iend, win_level;
        MPI_Win win_clusters_x, win_clusters_y, win_clusters_z, win_clusters_q, win_clusters_w;
        MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
        
        MPI_Win_create(tree_array->x_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_x_mid);
//        printf("rank %i created first tree_array windows.\n", rank);
        MPI_Win_create(tree_array->y_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_y_mid);
        MPI_Win_create(tree_array->z_mid,  numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_z_mid);
        MPI_Win_create(tree_array->radius, numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_radius);
//        printf("rank %i created first tree_array windows.\n", rank);
        MPI_Win_create(tree_array->numpar, numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_numpar);
//        printf("rank %i created first tree_array windows.\n", rank);
        MPI_Win_create(tree_array->ibeg,   numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_iend);
        MPI_Win_create(tree_array->iend,   numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_ibeg);
        MPI_Win_create(tree_array->level,  numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_level);
        printf("rank %i created tree_array windows.\n", rank);

        MPI_Win_create(clusters->x, numnodes*pointsPerCluster*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_x);
        MPI_Win_create(clusters->y, numnodes*pointsPerCluster*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_y);
        MPI_Win_create(clusters->z, numnodes*pointsPerCluster*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_z);
        MPI_Win_create(clusters->q, numnodes*pointsPerCluster*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_q);
        MPI_Win_create(clusters->w, numnodes*pointsPerCluster*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_clusters_w);
        printf("rank %i created clusters windows.\n", rank);

        MPI_Win_create(sources->x, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
        MPI_Win_create(sources->y, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
        MPI_Win_create(sources->z, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
        MPI_Win_create(sources->q, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
        MPI_Win_create(sources->w, troot->numpar*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);
        printf("rank %i created source windows.\n", rank);


    	// Perform MPI round robin, filling LET with remote data
		for (int procID = 1; procID < numProcs; ++procID) {

			int getFrom = (numProcs+rank-procID) % numProcs;

			// Allocate remote_tree_array
			struct tnode_array *remote_tree_array = NULL;
			remote_tree_array = malloc(sizeof(struct tnode_array));
			allocate_tree_array(remote_tree_array, numNodesOnProc[getFrom]); // start by allocating let_tree_array with size numnodes

			// Get remote_tree_array
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_x_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_y_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_z_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_radius);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_numpar);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_ibeg);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_iend);
            MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_level);
            
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
            
            MPI_Win_unlock(getFrom, win_x_mid);
            MPI_Win_unlock(getFrom, win_y_mid);
            MPI_Win_unlock(getFrom, win_z_mid);
            MPI_Win_unlock(getFrom, win_radius);
            MPI_Win_unlock(getFrom, win_numpar);
            MPI_Win_unlock(getFrom, win_ibeg);
            MPI_Win_unlock(getFrom, win_iend);
            MPI_Win_unlock(getFrom, win_level);
            printf("rank %i gotfrom %i.\n", rank,getFrom);
            MPI_Barrier(MPI_COMM_WORLD);
            printf("Rank %i remote_tree_array->x_mid[0] = %1.2e\n", rank, remote_tree_array->x_mid[0]);
            printf("Rank %i remote_tree_array->numpar[0] = %i\n", rank, remote_tree_array->numpar[0]);
            printf("Rank %i remote_tree_array->numpar[1] = %i\n", rank, remote_tree_array->numpar[1]);
            printf("Rank %i remote_tree_array->numpar[2] = %i\n", rank, remote_tree_array->numpar[2]);

			// Construct masks
			int *approx_list; int *direct_list; int *direct_ibeg_list; int *direct_length_list;
			make_vector(approx_list, numNodesOnProc[getFrom]);
			make_vector(direct_list, numNodesOnProc[getFrom]);
            make_vector(direct_ibeg_list, numNodesOnProc[getFrom]);
            make_vector(direct_length_list, numNodesOnProc[getFrom]);
            
			remote_interaction_lists(remote_tree_array, batches, approx_list, direct_list, numNodesOnProc[getFrom]);

			MPI_Barrier(MPI_COMM_WORLD);
			printf("Rank %i, approx_list[0]=%i\n", rank, approx_list[0]);
			printf("Rank %i, approx_list[1]=%i\n", rank, approx_list[1]);
			printf("Rank %i, direct_list[0]=%i\n", rank, direct_list[0]);
			printf("Rank %i, direct_list[1]=%i\n", rank, direct_list[1]);
			MPI_Barrier(MPI_COMM_WORLD);
		    // Reallocate structs to hold new data

			// Count number of unique clusters adding to LET
			int numberOfUniqueClusters =  0;
			int previousTreeArrayLength = let_tree_array_length;
			for (int i = 0; i < numNodesOnProc[getFrom]; ++i) {
				if ((approx_list[i] != -1) || (direct_list[i] != -1)) {
					numberOfUniqueClusters+=1;
					let_tree_array_length+=1;
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
			printf("Rank %i got here.  New let_tree_array_length = %i\n", rank, let_tree_array_length);
            reallocate_tree_array(let_tree_array, let_tree_array_length);
            MPI_Barrier(MPI_COMM_WORLD);
			printf("a Rank %i got here.\n", rank);
            int numberOfRemoteApprox = 0;
            int previous_let_clusters_length;
            
            int numberOfRemoteDirect = 0;
            int previous_let_sources_length;

            // Fill in LET tree array from Remote tree array.
            int appendCounter = 0;
            for (int i = 0; i < numNodesOnProc[getFrom]; ++i) {

            	previous_let_clusters_length = let_clusters_length;
            	previous_let_sources_length = let_sources_length;


				if ((approx_list[i] != -1) || (direct_list[i] != -1)) {
					let_tree_array->x_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->x_mid[i];
					let_tree_array->y_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->y_mid[i];
					let_tree_array->z_mid[previousTreeArrayLength + appendCounter] = remote_tree_array->z_mid[i];
					let_tree_array->radius[previousTreeArrayLength + appendCounter] = remote_tree_array->radius[i];
					let_tree_array->numpar[previousTreeArrayLength + appendCounter] = remote_tree_array->numpar[i];
                    let_tree_array->level[previousTreeArrayLength + appendCounter] = remote_tree_array->level[i];
                    
                    if (approx_list[i] != -1) {
                        let_tree_array->cluster_ind[previousTreeArrayLength + appendCounter] = previousTreeArrayLength + numberOfRemoteApprox;
                        let_clusters_length += pointsPerCluster;
                        numberOfRemoteApprox++;
                    }
                    
                    if (direct_list[i] != -1) {
                        
                        // Set the beginning and ending particle indices for the associated nodes in the local sources list
                        let_tree_array->ibeg[previousTreeArrayLength + appendCounter] = let_sources_length;
                        let_tree_array->iend[previousTreeArrayLength + appendCounter] = let_sources_length + remote_tree_array->numpar[i] - 1; // this minus 1 s
                        let_sources_length += remote_tree_array->numpar[i];
                        
                        // Determine displacements and lengths for getting prticles from remote sources list
                        direct_ibeg_list[numberOfRemoteDirect] = remote_tree_array->ibeg[i];
                        direct_length_list[numberOfRemoteDirect] = remote_tree_array->numpar[i];
                        numberOfRemoteDirect++;
                    }
                    
					appendCounter+=1;
				}
			}
            MPI_Barrier(MPI_COMM_WORLD);
            printf("b Rank %i got here.\n", rank);
            reallocate_cluster(let_clusters, let_clusters_length);
            reallocate_sources(let_sources, let_sources_length);
            
			// Use masks to get remote data
            
            // Leighton, approx_list is of form [a0, a1, a2, -1, -1, -1, ...] where we need clusters a0, a1, a2 for the let_clusters array.
            //           The first numberOfRemoteApprox are needed.  They should be stored at let_clusters[previous_let_clusters_length].
            //              direct_list is of form [d0, d1, d2, -1, -1, -1, ...] where we need sources d0, d1, d2 for the let_sources array.
            //           The first numberOfRemoteDirect are needed.  They should be stored at let_sources[previous_let_sources_length].

            
            MPI_Datatype approx_type, direct_type;
            
            MPI_Type_create_indexed_block(numberOfRemoteApprox, pointsPerCluster, approx_list, MPI_DOUBLE, &approx_type);
            MPI_Type_commit(&approx_type);
            int new_sources_length = let_sources_length - previous_let_sources_length;
            MPI_Type_indexed(numberOfRemoteDirect, direct_length_list, direct_ibeg_list, MPI_DOUBLE, &direct_type);
            MPI_Type_commit(&direct_type);

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
            

            MPI_Barrier(MPI_COMM_WORLD);
			printf("f Rank %i got here.\n", rank);
			fflush(stdout);
			MPI_Barrier(MPI_COMM_WORLD);


            MPI_Get(&(let_clusters->x[previous_let_clusters_length]), numberOfRemoteApprox * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type, win_clusters_x);
            MPI_Get(&(let_clusters->y[previous_let_clusters_length]), numberOfRemoteApprox * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type, win_clusters_y);
            MPI_Get(&(let_clusters->z[previous_let_clusters_length]), numberOfRemoteApprox * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type, win_clusters_z);
            MPI_Get(&(let_clusters->q[previous_let_clusters_length]), numberOfRemoteApprox * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type, win_clusters_q);
            MPI_Get(&(let_clusters->w[previous_let_clusters_length]), numberOfRemoteApprox * pointsPerCluster, MPI_DOUBLE,
                    getFrom, 0, 1, approx_type, win_clusters_w);
            
            MPI_Get(&(let_sources->x[previous_let_sources_length]), new_sources_length, MPI_DOUBLE,
                    getFrom, 0, 1, direct_type, win_sources_x);
            MPI_Get(&(let_sources->y[previous_let_sources_length]), new_sources_length, MPI_DOUBLE,
                    getFrom, 0, 1, direct_type, win_sources_y);
            MPI_Get(&(let_sources->z[previous_let_sources_length]), new_sources_length, MPI_DOUBLE,
                    getFrom, 0, 1, direct_type, win_sources_z);
            MPI_Get(&(let_sources->q[previous_let_sources_length]), new_sources_length, MPI_DOUBLE,
                    getFrom, 0, 1, direct_type, win_sources_q);
            MPI_Get(&(let_sources->w[previous_let_sources_length]), new_sources_length, MPI_DOUBLE,
                    getFrom, 0, 1, direct_type, win_sources_w);
            
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
            printf("Rank %i, clusters->x[previous_let_clusters_length+1] = %1.2e\n", rank, clusters->x[previous_let_clusters_length+1]);
//            printf("Rank %i, clusters->y[previous_let_clusters_length+2] = %1.2e\n", rank, clusters->y[previous_let_clusters_length+2]);
            printf("Rank %i, sources->x[previous_let_sources_length+1] = %1.2e\n", rank, sources->x[previous_let_sources_length+1]);
            MPI_Barrier(MPI_COMM_WORLD);

            
			free_vector(approx_list);
			free_vector(direct_list);
            free_vector(direct_ibeg_list);
            free_vector(direct_length_list);

			free_tree_array(remote_tree_array);
		} // end loop over numProcs


    	// Compute interaction lists based on LET
        make_vector(tree_inter_list, batches->num * let_tree_array->numnodes);
        make_vector(direct_inter_list, batches->num * let_tree_array->numnodes);
    	pc_make_interaction_list(let_tree_array, batches, tree_inter_list,  direct_inter_list);

    	printf("Rank %i just made interaction list based off of its LET.\n", rank);



        // After filling LET, call interaction_list_treecode
    	time1 = MPI_Wtime(); // start timer for tree evaluation
        if (pot_type == 0) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=0 (Coulomb).\n");
            pc_interaction_list_treecode(let_tree_array, let_clusters, batches,
                                         tree_inter_list, direct_inter_list, let_sources, targets,
                                         tpeng, tEn, numDevices, numThreads);
            if (verbosity>0) printf("Exiting particle-cluster, pot_type=0 (Coulomb).\n");

        } else if (pot_type == 1) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=1 (Yukawa).  Need to modify to call on let_tree_array, let_clusters, etc. ????\n");
            pc_interaction_list_treecode_yuk(let_tree_array, clusters, batches,
                                             tree_inter_list, direct_inter_list, sources, targets,
                                             tpeng, kappa, tEn, numDevices, numThreads);

        } else if (pot_type == 2) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=2 (Coulomb w/ singularity subtraction).\n");
            pc_treecode_coulomb_SS(troot, batches, sources, targets,clusters,
                                   kappa, tpeng, tEn, numDevices, numThreads);
//          pc_interaction_list_treecode_Coulomb_SS(tree_array, clusters, batches,
//                                                  tree_inter_list, direct_inter_list, sources, targets,
//                                                  tpeng, kappa, tEn, numDevices, numThreads);

        } else if (pot_type == 3) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=3 (Yukawa w/ singularity subtraction).\n");
            pc_treecode_yuk_SS(troot, batches, sources, targets,clusters,
                               kappa, tpeng, tEn, numDevices, numThreads);

        } else if (pot_type == 4) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=4 (Coulomb Hermite).\n");
            pc_interaction_list_treecode_hermite_coulomb(let_tree_array, clusters, batches,
                                                    tree_inter_list, direct_inter_list, sources, targets,
                                                    tpeng, tEn, numDevices, numThreads);

        }else if (pot_type == 5) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=4 (Yukawa Hermite).\n");
            pc_interaction_list_treecode_hermite_yukawa(let_tree_array, clusters, batches,
                                                    tree_inter_list, direct_inter_list, sources, targets,
                                                    tpeng, kappa, tEn, numDevices, numThreads);

        }else if (pot_type == 6) {
            if (verbosity>0) printf("Entering particle-cluster, pot_type=6 (Coulomb Hermite w/ singularity subtraction).\n");
            pc_treecode_hermite_coulomb_SS(troot, batches, sources, targets,clusters,
                                           kappa, tpeng, tEn, numDevices, numThreads);
        }
        
        reorder_energies(batches->reorder, targets->num, tEn);
    }


    time2 = MPI_Wtime();  // end time for tree evaluation
    timetree[3] = time2-time1; //+ timetree[0];

    if (verbosity>0) printf("Time to compute: %f\n", time2-time1);


    cleanup(troot);

    // free interaction lists
    free_vector(tree_inter_list);
    free_vector(direct_inter_list);

    // free clusters
    free_vector(clusters->x);
    free_vector(clusters->y);
    free_vector(clusters->z);
    free_vector(clusters->q);
    free_vector(clusters->w);
    free(clusters);

    // free tree_array
    free_vector(tree_array->ibeg);
    free_vector(tree_array->iend);
    free_vector(tree_array->x_mid);
    free_vector(tree_array->y_mid);
    free_vector(tree_array->z_mid);
    free_vector(tree_array->x_min);
    free_vector(tree_array->y_min);
    free_vector(tree_array->z_min);
    free_vector(tree_array->x_max);
    free_vector(tree_array->y_max);
    free_vector(tree_array->z_max);
    free(tree_array);

    // free target batches
    free_vector(batches->reorder);
    free_matrix(batches->index);
    free_matrix(batches->center);
    free_vector(batches->radius);
    free(batches);

    MPI_Barrier(MPI_COMM_WORLD);
    printf("Exiting treedriver.\n");
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);

    return;

} /* END function treecode */

