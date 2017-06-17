#include <stdio.h>
#include <mpi.h>
#include <limits.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"


/* definition of primary treecode driver */

void treecode_grid(double *xS, double *yS, double *zS, double *qS, 
                   double *xyzminmax, int *xyzdim, int numparsS,
                   double *tEn, double *tpeng, int order, double theta, 
                   int maxparnode, double *timetree, int pot_type, double kappa)
{

    /* local variables */
    struct tnode *troot = NULL;
    int level;
    int numparsT;
    int xyzind[6];

    /* date and time */
    double time1, time2;

    
    /* set global variables to track tree levels during construction */
    level = 0;
    numleaves = 0;
    minlevel = INT_MAX;
    minpars = INT_MAX;
    xdiv = 0; ydiv = 0; zdiv = 0;

    maxlevel = 0;
    maxpars = 0;

    numparsT = xyzdim[0] * xyzdim[1] * xyzdim[2];

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    time1 = MPI_Wtime();
    setup_grid(xyzminmax, xyzdim, xyzind, order, theta);
    
    printf("Creating tree... \n\n");
    cp_create_tree_n0_grid(&troot, maxparnode, xyzminmax, xyzdim, xyzind, level);

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
    printf("        maxparnode (iflag 0): %d\n", maxparnode);
    printf("               tree maxlevel: %d\n", maxlevel);
    printf("               tree minlevel: %d\n", minlevel);
    printf("                tree maxpars: %d\n", maxpars);
    printf("                tree minpars: %d\n", minpars);
    printf("            number of leaves: %d\n", numleaves);
    printf("        number of x,y,z divs: %d, %d, %d\n", xdiv, ydiv, zdiv);


    time1 = MPI_Wtime();

    printf("Calling tree computation... \n\n");
    cp_treecode_grid(troot, xS, yS, zS, qS, tpeng, tEn, 
                     numparsS, numparsT, &timetree[1]);

    time2 = MPI_Wtime();
    timetree[3] = time2-time1 + timetree[0];

    //printf("       Tree building time (s): %f\n", *timetree - totaltime);
    //printf("    Tree computation time (s): %f\n\n", totaltime);
    printf("Deallocating tree structure... \n\n");

    cleanup(troot);

    return;

} /* END function treecode */

