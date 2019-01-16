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
                double *tEn, double *tpeng, double *timetree)
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
    double time1, time2;

    
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
    
    printf("Creating tree... \n\n");

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    if (tree_type == 0) {
        if (pot_type == 0) {
            setup(targets, order, theta, xyzminmax);
        } else if (pot_type == 1) {
            setup_yuk(targets, order, theta, xyzminmax);
        }
        
        cp_create_tree_n0(&troot, targets, 1, targets->num,
                          maxparnode, xyzminmax, level);
        
        setup_batch(&batches, batch_lim, targets, batch_size);
        create_source_batch(batches, sources, 1, sources->num,
                            batch_size, batch_lim);
        
    } else if (tree_type == 1) {
    
        if (pot_type == 0) {
            setup(sources, order, theta, xyzminmax);
        } else if (pot_type == 1) {
            setup_yuk(sources, order, theta, xyzminmax);
        }
        
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

        setup_batch(&batches, batch_lim, targets, batch_size);
        create_target_batch(batches, targets, 1, targets->num,
                            batch_size, batch_lim);
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

    time1 = MPI_Wtime();

    /* Copy source arrays to GPU */
	//#pragma acc data copyin(xS[numparsS], yS[numparsS], zS[numparsS], qS[numparsS])


    if (tree_type == 0) {
        if (pot_type == 0) {
            cp_treecode(troot, batches, sources, targets,
                        tpeng, tEn, &timetree[1]);
        } else if (pot_type == 1) {
            cp_treecode_yuk(troot, batches, sources, targets, kappa,
                            tpeng, tEn, &timetree[1]);
        }
    } else if (tree_type == 1) {
        if (pot_type == 0) {
        
            make_matrix(tree_inter_list, batches->num, numnodes);
            make_matrix(direct_inter_list, batches->num, numleaves);
            
            pc_make_interaction_list(troot, batches, tree_inter_list, direct_inter_list);
            
            pc_treecode(troot, batches, sources, targets, tpeng, tEn);
            
        } else if (pot_type == 1) {
            pc_treecode_yuk(troot, batches, sources, targets,
                            kappa, tpeng, tEn);
        }
        
        reorder_energies(batches->reorder, targets->num, tEn);
    }


    time2 = MPI_Wtime();
    timetree[3] = time2-time1 + timetree[0];

    printf("Deallocating tree structure... \n\n");

    cleanup(troot);

    return;

} /* END function treecode */

