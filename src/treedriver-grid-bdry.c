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

void treecode_grid_bdry(double *xS, double *yS, double *zS, double *qS, 
                        double *xyzminmax, int *xyzdim, double zyx, int dir, int numparsS,
                        double *tEn, double *tpeng, int order, double theta, 
                        int maxparnode, double *timetree, int pot_type, double kappa)
{

    /* local variables */
    struct tnode *troot = NULL;
    int level;
    int numparsT;
    int xyzind[4];

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

    numparsT = xyzdim[0] * xyzdim[1];

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    time1 = MPI_Wtime();
    setup_grid_bdry(xyzminmax, xyzdim, xyzind, order, theta);

    cp_create_tree_n0_grid_bdry(&troot, maxparnode, xyzminmax, xyzdim, xyzind, level);

    time2 = MPI_Wtime();
    timetree[0] = time2-time1;

/*
    printf("Tree created.\n\n");
    printf("Tree information: \n\n");

    printf("                      numpar: %d\n", troot->numpar);
    printf("                       x_mid: %e\n", troot->x_mid);
    printf("                       y_mid: %e\n", troot->y_mid);
    printf("                  z pos, dir: %e %d\n", zyx, dir);
    printf("                      radius: %f\n\n", troot->radius);
    printf("                       x_len: %e\n", troot->x_max - troot->x_min);
    printf("                       y_len: %e\n", troot->y_max - troot->y_min);
    printf("                      torder: %d\n", torder);
    printf("                       theta: %f\n", theta);
    printf("        maxparnode (iflag 0): %d\n", maxparnode);
    printf("               tree maxlevel: %d\n", maxlevel);
    printf("               tree minlevel: %d\n", minlevel);
    printf("                tree maxpars: %d\n", maxpars);
    printf("                tree minpars: %d\n", minpars);
    printf("            number of leaves: %d\n", numleaves);
    printf("        number of x,y divs: %d, %d\n", xdiv, ydiv);
*/

    time1 = MPI_Wtime();

    cp_treecode_grid_bdry(troot, zyx, dir, xS, yS, zS, qS, tpeng, tEn,
                     numparsS, numparsT, &timetree[1]);

    time2 = MPI_Wtime();
    timetree[3] = time2-time1 + timetree[0];

    //printf("       Tree building time (s): %f\n", *timetree - totaltime);
    //printf("    Tree computation time (s): %f\n\n", totaltime);

    cleanup(troot);

    return;

} /* END function treecode */

