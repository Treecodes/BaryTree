#include <time.h>
//#include <sys/time.h>
#include <stdio.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"
#include "tree-cp.h"

#include "treedriver.h"
/* definition of primary treecode driver */

void treecode(double *xS, double *yS, double *zS, double *qS, 
              double *xT, double *yT, double *zT,
              int numparsS, int numparsT, double *tEn, double *tpeng, 
              int order, double theta, int shrinkS, int shrinkT, 
              int maxparnodeS, int maxparnodeT, double *timetree,
              int treelevelS, int treelevelT, int iflagS, int iflagT,
              int pot_type, double kappa)
{

        /* local variables */
        struct tnode *trootS, *trootT;
        int i, level, err;
        double *xyzminmax;

        /* date and time */
        time_t time1, time2;
        double totaltime;


        make_vector(xyzminmax, 6);

        /* call setup to allocate arrays for Taylor expansions and setup global vars */
        if (pot_type == 0)
                setup(xT, yT, zT, numparsT, order, theta, xyzminmax);
        else if (pot_type == 1)
                setup_yuk(xT, yT, zT, numparsT, order, theta, xyzminmax);


        time1 = time(NULL);

        /* set global variables to track tree levels during construction */
        level = 0;
        minlevel = 50000;

        printf("Creating tree... \n\n");


        if (iflagT == 0) {
                maxlevel = 0;
                create_tree_n0(&trootT, 1, numparsT, xT, yT, zT,
                               shrinkT, maxparnodeT, xyzminmax,
                               level, numparsT);
        } else {
                maxlevel = treelevelT;
                create_tree_lv(&trootT, 1, numparsT, xT, yT, zT,
                               shrinkT, treelevelT, xyzminmax,
                               level, numparsT);
        }


        time2 = time(NULL);
        totaltime = difftime(time2, time1);
        *timetree = totaltime;


        printf("Tree created.\n\n");
        printf("Tree information: \n\n");

        printf("       numpar: %d\n", trootT->numpar);
        printf("        x_mid: %e\n", trootT->x_mid);
        printf("        y_mid: %e\n", trootT->y_mid);
        printf("        z_mid: %e\n\n", trootT->z_mid);
        printf("       radius: %f\n", trootT->radius);
        printf("       torder: %d\n", torder);
        printf("        theta: %f\n", theta);
        printf("       shrink: %d\n", shrinkT);
        printf("   maxparnode: %d\n", maxparnodeT);
        printf("        iflag: %d\n", iflagT);
        printf("tree maxlevel: %d\n\n", treelevelT);


        time1 = time(NULL);


        if (pot_type == 0) {
                cp_treecode(trootT, xS, yS, zS, qS, xT, yT, zT,
                            tpeng, tEn, numparsS, numparsT);

        } else if (pot_type == 1) {
                cp_treecode_yuk(trootT, xS, yS, zS, qS, xT, yT, zT,
                                tpeng, tEn, numparsS, numparsT,
                                kappa);
        }

        time2 = time(NULL);
        totaltime = difftime(time2, time1);
        *timetree = *timetree + totaltime;

        printf("Deallocating tree structure... \n\n");

        cleanup(trootT);

        return;

} /* END function cp_treecode */

