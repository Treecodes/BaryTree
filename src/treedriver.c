#include <stdio.h>
#include <mpi.h>
#include <limits.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"


/* definition of primary treecode driver */

void treedriver(struct particles *sources, struct particles *targets,
                int order, double theta, int maxparnode, int batch_size,
                int pot_type, double kappa, int tree_type,
                double *tEn, double *tpeng, double *timetree, int numDevices)
{

    /* local variables */
    struct tnode *troot = NULL;
    int level;
    double xyzminmax[6];
    
//    int i, j;
    
    /* batch variables */
    struct batch *batches = NULL;
    double batch_lim[6];
    
    /* date and time */
    double time1, time2, timeFillClusters1, timeFillClusters2;

    
    time1 = MPI_Wtime();
    
    level = 0;
    numleaves = 0;
    minlevel = INT_MAX;
    minpars = INT_MAX;
    maxlevel = 0;
    maxpars = 0;
    
    int **tree_inter_list;
    int **direct_inter_list;
    struct tnode_array *tree_array = NULL;
    numnodes = 0;
    struct particles *clusters = NULL;
	clusters = malloc(sizeof(struct particles));

#pragma omp parallel num_threads(numDevices)
	{
        acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
        acc_init(acc_get_device_type());
	}
    
    printf("Creating tree... \n\n");

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    if (tree_type == 0) {
        if (pot_type == 0) {
            setup(targets, order, theta, xyzminmax);
        }
//        } else if (pot_type == 1) {
//            setup_yuk(targets, order, theta, xyzminmax);
//        }
        
        cp_create_tree_n0(&troot, targets, 1, targets->num,
                          maxparnode, xyzminmax, level);
        
        setup_batch(&batches, batch_lim, targets, batch_size);
        create_source_batch(batches, sources, 1, sources->num,
                            batch_size, batch_lim);
        
    } else if (tree_type == 1) {
    	printf("Treetype %i\n", tree_type);
//        if (pot_type == 0) {
//            setup(sources, order, theta, xyzminmax);
//        } else if (pot_type == 1) {
////            setup_yuk(sources, order, theta, xyzminmax);
//            setup(sources, order, theta, xyzminmax);  // call the non-Yukawa setup.  This has the Chebyshev parts.
//        }

    	setup(sources, order, theta, xyzminmax);
        
        pc_create_tree_n0(&troot, sources, 1, sources->num,
                          maxparnode, xyzminmax, level);
        
        tree_array = malloc(sizeof(struct tnode_array));
        tree_array->numnodes = numnodes;
        make_vector(tree_array->ibeg, numnodes);
        make_vector(tree_array->iend, numnodes);
        make_vector(tree_array->x_mid, numnodes);
        make_vector(tree_array->y_mid, numnodes);
        make_vector(tree_array->z_mid, numnodes);


        pc_create_tree_array(troot, tree_array);

//        printf("Entering setup_batch.\n");
        setup_batch(&batches, batch_lim, targets, batch_size);
//        printf("Exiting setup_batch.\n");
//        printf("Entering create_target_batch.\n");
        create_target_batch(batches, targets, 1, targets->num,batch_size, batch_lim);
//        printf("Exiting create_target_batch.\n");


//#pragma acc data region copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
//		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num])
//        {

        timeFillClusters1 = MPI_Wtime();
        if (        (pot_type == 0) || (pot_type==1)) {
        	fill_in_cluster_data(clusters, sources, troot, order, numDevices);
        }else if  ( (pot_type == 2) || (pot_type==3)){
        	printf("Calling fill_in_cluster_data_SS().\n");
			fill_in_cluster_data_SS(clusters, sources, troot, order);
        }else if  ( (pot_type == 4) || (pot_type==5)){
        	printf("Calling fill_in_cluster_data_hermite().\n");
			fill_in_cluster_data_hermite(clusters, sources, troot, order);
		}else if  ( (pot_type == 6) || (pot_type==7)){
			printf("Calling fill_in_cluster_data_hermite_SS().\n");
			fill_in_cluster_data_hermite_SS(clusters, sources, troot, order);
		}
        timeFillClusters2 = MPI_Wtime();
        timeFillClusters1 = timeFillClusters2-timeFillClusters1;
        printf("Time to compute modified weights(s):  %f\n", timeFillClusters1);

    }

    time2 = MPI_Wtime();
    timetree[0] = time2-time1;

    printf("Tree created.\n\n");
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


    time1 = MPI_Wtime();

    /* Copy source and target arrays to GPU */
//#pragma acc data copyin(xS[numparsS], yS[numparsS], zS[numparsS], qS[numparsS], wS[numparsS], \
//		xT[numparsT], yT[numparsT], zT[numparsT], qT[numparsT])


    if (tree_type == 0) {
    	printf("\nTook hooks out for cp_treecode...\n");
//        if (pot_type == 0) {
//            cp_treecode(troot, batches, sources, targets,
//                        tpeng, tEn, &timetree[1]);
//        }
//        } else if (pot_type == 1) {
//            cp_treecode_yuk(troot, batches, sources, targets, kappa,
//                            tpeng, tEn, &timetree[1]);
//        }
    } else if (tree_type == 1) {
    	make_matrix(tree_inter_list, batches->num, numnodes);
		make_matrix(direct_inter_list, batches->num, numleaves);

		pc_make_interaction_list(troot, batches, tree_inter_list, direct_inter_list);
        if (pot_type == 0) {
        	printf("Entering tree_type=1 (particle-cluster), pot_type=0 (Coulomb).\n");
            pc_treecode(troot, batches, sources, targets, clusters, tpeng, tEn, numDevices);
        } else if (pot_type == 1) {
        	printf("Entering tree_type=1 (particle-cluster), pot_type=1 (Yukawa).\n");
            pc_treecode_yuk(troot, batches, sources, targets, clusters,
                            kappa, tpeng, tEn, numDevices);
        }else if (pot_type == 2) {
        	printf("Entering tree_type=1 (particle-cluster), pot_type=2 (Coulomb w/ singularity subtraction).\n");
        	pc_treecode_coulomb_SS(troot, batches, sources, targets,clusters,
        	                            kappa, tpeng, tEn, numDevices);
        }else if (pot_type == 3) {
        	printf("Entering tree_type=1 (particle-cluster), pot_type=3 (Yukawa w/ singularity subtraction).\n");
        	pc_treecode_yuk_SS(troot, batches, sources, targets,clusters,
        	                            kappa, tpeng, tEn, numDevices);
        }else if (pot_type == 4) {
        	printf("Entering tree_type=1 (particle-cluster), pot_type=4 (Coulomb Hermite).\n");
        	pc_treecode_hermite(troot, batches, sources, targets,clusters, tpeng, tEn, numDevices);
        }
        
        reorder_energies(batches->reorder, targets->num, tEn);
    }


    time2 = MPI_Wtime();
    timetree[3] = time2-time1 + timetree[0];

//    printf("Deallocating tree structure... \n\n");
//    printf("Time1: %12.5f\n",time1);
//    printf("Time2: %12.5f\n",time2);
//    printf("Timetree[0]: %12.5f\n",timetree[0]);
//    printf("Time: %12.5f\n\n",timetree[3]);

    cleanup(troot);
//    printf("Finished cleanup of troot.\n");
    return;

} /* END function treecode */

