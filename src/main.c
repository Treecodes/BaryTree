#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "sort.h"


/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* runtime parameters */
    int numparsS, numparsT, order;
    int maxparnode, treelevel;
    int iflag, pot_type, tree_type;
    int pflag, sflag, dflag, gflag;

    double theta, temp;
    double kappa;

    /* source particles */
    double *xS = NULL;
    double *yS = NULL;
    double *zS = NULL;
    double *qS = NULL;
    
    /* target particles */
    double *xT = NULL;
    double *yT = NULL;
    double *zT = NULL;
    int *iT = NULL;

    /* exact energy, treecode energy */
    double *denergyglob = NULL;
    double *tenergyglob = NULL;
    double *tenergy = NULL;

    /* for potential energy calculation */
    double tpeng = 0;
    double dpengglob = 0;
    double tpengglob = 0;

    /* insert variables for date-time calculation? */
    double time_direct, time_tree[4], time_preproc;
    double time_tree_glob[3][4];
    double time1, time2;

    /* input and output files */
    char *sampin1 = NULL;
    char *sampin2 = NULL;
    char *sampin3 = NULL;
    char *sampout = NULL;
    FILE *fp;

    /* variables for error calculations */
    double inferr, relinferr, n2err, reln2err;

    /* local variables */
    int i, j;
    int numparsTloc, maxparsTloc, numparsSloc;
    double buf[5];
    
    int *displs = NULL;
    int *scounts = NULL;
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;


    /* Executable statements begin here */

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       sampin1:  sources input file \n");               // "S10000.txt"
            printf("       sampin2:  targets input file \n");               // "T1000000.txt"
            printf("       sampin3:  direct calc potential input file \n"); // "ex_s4_t6.txt"
            printf("       sampout:  tree calc potential output file \n");  // "out.txt"
            printf("      numparsS:  number of sources \n");                // 10000
            printf("      numparsT:  number of targets \n");                // 1000000
            printf("         theta:  multipole acceptance criterion \n");   // 0.75
            printf("         order:  order of treecode Taylor expansion \n");        // 20
            printf("     tree_type:  0--cluster-particle, 1--particle-cluster \n");  // 0
            printf("    maxparnode:  maximum particles in leaf \n");                 // 500
            printf("     treelevel:  maximum tree levels \n");                       // 5
            printf("         iflag:  0--use maxparnode, 1--use treelevel \n");       // 0
            printf("         kappa:  screened Coulomb parameter \n");                // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");           // 1
            printf("         pflag:  distribute 0--targets, 1--sources \n");         // 0
            printf("         sflag:  on distr 0--sort, 1--no sort \n");              // 0
            printf("         dflag:  if sorted, direction 0--x, 1--y, 2--z \n");     // 0
        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampin3 = argv[3];
    sampout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    theta = atof(argv[7]);
    order = atoi(argv[8]);
    tree_type = atoi(argv[9]);
    maxparnode = atoi(argv[10]);
    treelevel = atoi(argv[11]);
    iflag = atoi(argv[12]);
    kappa = atof(argv[13]);
    pot_type = atoi(argv[14]);
    pflag = atoi(argv[15]);
    sflag = atoi(argv[16]);
    dflag = atoi(argv[17]);
    gflag = atoi(argv[18]);

    numparsTloc = numparsT;
    numparsSloc = numparsS;
    
    time1 = MPI_Wtime();

    if (pflag == 0) {

        numparsTloc = numparsT/p;
        maxparsTloc = numparsTloc + (numparsS - (numparsS / p) * p);
    
        make_vector(xS,numparsS);
        make_vector(yS,numparsS);
        make_vector(zS,numparsS);
        make_vector(qS,numparsS);
    
        if (rank == 0) {
            make_vector(xT,numparsT);
            make_vector(yT,numparsT);
            make_vector(zT,numparsT);
            make_vector(iT,numparsT);
        
            make_vector(tenergy,numparsT);
            make_vector(tenergyglob,numparsT);
            make_vector(denergyglob,numparsT);
        
            make_vector(displs,p);
            make_vector(scounts,p);
        } else {
            make_vector(xT,numparsTloc);
            make_vector(yT,numparsTloc);
            make_vector(zT,numparsTloc);
        
            make_vector(tenergy,numparsTloc);
        }

    
        if (rank == 0) {

            MPI_File_open(MPI_COMM_SELF, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
            MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
            for (i = 0; i < numparsT; i++) {
                MPI_File_read(fpmpi, buf, 3, MPI_DOUBLE, &status);
                xT[i] = buf[0];
                yT[i] = buf[1];
                zT[i] = buf[2];
                iT[i] = i;
            }
            MPI_File_close(&fpmpi);
        
            if (sflag == 0) sortTargets(xT, yT, zT, iT, numparsT, dflag);
            if (sflag == 1 && gflag == 1) interleaveGridTargets(xT, yT, zT, iT, numparsT, p);

            MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
            MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
            MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);

            MPI_File_read(fpmpi, denergyglob, numparsT, MPI_DOUBLE, &status);

            MPI_File_close(&fpmpi);
    
            scounts[0] = 0;
            displs[0] = 0;
    
            for (i=1; i < p; i++) {
                scounts[i] = numparsTloc;
                displs[i] = (i-1) * numparsTloc;
            }
        }
    
        MPI_Scatterv(&xT[maxparsTloc], scounts, displs, MPI_DOUBLE,
                     xT, numparsTloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&yT[maxparsTloc], scounts, displs, MPI_DOUBLE,
                     yT, numparsTloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(&zT[maxparsTloc], scounts, displs, MPI_DOUBLE,
                     zT, numparsTloc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

        /* Reading in coordinates and charges for the source particles*/
        MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        for (i = 0; i < numparsS; i++) {
            MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
            xS[i] = buf[0];
            yS[i] = buf[1];
            zS[i] = buf[2];
            qS[i] = buf[3];
        }
        MPI_File_close(&fpmpi);

        if (rank == 0) numparsTloc = maxparsTloc;

    } else if (pflag == 1) {

        numparsSloc = numparsS / p;
        if (rank == p-1) numparsSloc += (numparsS - (numparsS / p) * p);
    
        make_vector(xS,numparsSloc);
        make_vector(yS,numparsSloc);
        make_vector(zS,numparsSloc);
        make_vector(qS,numparsSloc);
    
        make_vector(xT,numparsT);
        make_vector(yT,numparsT);
        make_vector(zT,numparsT);
        make_vector(iT,numparsT);
        make_vector(tenergy,numparsT);
        
        /* Reading in coordinates and charges for the source particles*/
        MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
        MPI_File_seek(fpmpi, (MPI_Offset)(rank*(numparsS/p)*4*sizeof(double)), MPI_SEEK_SET);
        for (i = 0; i < numparsSloc; i++) {
            MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
            xS[i] = buf[0];
            yS[i] = buf[1];
            zS[i] = buf[2];
            qS[i] = buf[3];
        }
        MPI_File_close(&fpmpi);

        /* Reading in coordinates for target particles*/
        MPI_File_open(MPI_COMM_SELF, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        for (i = 0; i < numparsT; i++) {
            MPI_File_read(fpmpi, buf, 3, MPI_DOUBLE, &status);
            xT[i] = buf[0];
            yT[i] = buf[1];
            zT[i] = buf[2];
            iT[i] = i;
        }
        MPI_File_close(&fpmpi);

        if (rank == 0) {
            make_vector(tenergyglob,numparsT);
            make_vector(denergyglob,numparsT);

            MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
            MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
            MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
            MPI_File_read(fpmpi, denergyglob, numparsT, MPI_DOUBLE, &status);
            MPI_File_close(&fpmpi);
        }

    }

    time2 = MPI_Wtime();
    time_preproc = time2 - time1;


    /* Calling main treecode subroutine to calculate approximate energy */
    treecode(xS, yS, zS, qS, xT, yT, zT, numparsSloc, numparsTloc,
             tenergy, &tpeng, order, theta, 1, maxparnode, time_tree,
             treelevel, iflag, pot_type, kappa, tree_type);

    
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (rank == 0) dpengglob = sum(denergyglob, numparsT);
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    
    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n\n", time_direct);
        printf("              Pre-process time (s):  %f\n", time_preproc);
        printf("      Min, Max tree setup time (s):  %f, %f\n", time_tree_glob[0][0],
                                                                time_tree_glob[1][0]);
        if (tree_type == 0) {
            printf("             Min, Max cp1 time (s):  %f, %f\n", time_tree_glob[0][1],
                                                                    time_tree_glob[1][1]);
            printf("             Min, Max cp2 time (s):  %f, %f\n", time_tree_glob[0][2],
                                                                    time_tree_glob[1][2]);
        }
        
        printf("      Min, Max total tree time (s):  %f, %f\n\n", time_tree_glob[0][3],
                                                                  time_tree_glob[1][3]);
        printf(" Preproc + Max total tree time (s):  %f \n\n", time_tree_glob[1][3] + time_preproc);
        
        //printf("                 Avg tree time (s):  %f\n\n", time_tree_tot/(double)p);
        //printf("         Direct : Tree on %d procs:  %f\n\n",
        //       p, time_direct/(time_tree_max*(double)p));

    
        /* Printing error in potential energy and potential energy */
        printf("           Direct potential energy:  %f\n", dpengglob);
        printf("             Tree potential energy:  %f\n\n", tpengglob);
    
        printf("Absolute error for total potential:  %e\n",
               fabs(tpengglob-dpengglob));
        printf("Relative error for total potential:  %e\n\n",
               fabs((tpengglob-dpengglob)/dpengglob));
    }
    
    
    /* Computing pointwise potential errors */
    if (pflag == 0) { 

        if (rank == 0) {
            scounts[0] = maxparsTloc;
            displs[0] = 0;
            for (i=1; i < p; i++) {
                scounts[i] = numparsTloc;
                displs[i] = maxparsTloc + (i-1) * numparsTloc;
            }
        }

        MPI_Gatherv(tenergy, numparsTloc, MPI_DOUBLE, tenergyglob,
                    scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    } else if (pflag == 1) {
        MPI_Reduce(tenergy, tenergyglob, numparsT, 
                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    if (rank == 0)
    {
        inferr = 0.0;
        relinferr = 0.0;
        n2err = 0.0;
        reln2err = 0.0;

        for (j = 0; j < numparsT; j++) {
            temp = fabs(denergyglob[iT[j]] - tenergyglob[j]);
        
            if (temp >= inferr)
                inferr = temp;

            if (fabs(denergyglob[j]) >= relinferr)
                relinferr = fabs(denergyglob[j]);

            n2err = n2err + pow(denergyglob[iT[j]] - tenergyglob[j], 2.0);
            reln2err = reln2err + pow(denergyglob[j], 2.0);
        }

        relinferr = inferr / relinferr;
        reln2err = sqrt(n2err / reln2err);
        n2err = sqrt(n2err);

        printf("Absolute inf norm error in potential:  %e \n", inferr);
        printf("Relative inf norm error in potential:  %e \n\n", relinferr);
        printf("  Absolute 2 norm error in potential:  %e \n", n2err);
        printf("  Relative 2 norm error in potential:  %e \n\n", reln2err);
    }
    
    
    if (rank == 0) {
        fp = fopen(sampout, "a");
        fprintf(fp, "%s \t %s \t %s \t %d \t %d \t %f \t %d \t %d \t %d \t"
                "%d \t %d \t %f \t %d \t %d \t %d \t"
                "%d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t"
                "%f \t %f \t %f \t %f \t %f \t %f \t %f \t"
                "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",
                sampin1, sampin2, sampin3, numparsS, numparsT,
                theta, order, tree_type, maxparnode, treelevel, iflag,
                kappa, pot_type, sflag, pflag, //2 ends
                p, time_preproc,
                time_tree_glob[0][0], time_tree_glob[1][0],
                time_tree_glob[2][0]/(double)p,
                time_tree_glob[0][1], time_tree_glob[1][1],
                time_tree_glob[2][1]/(double)p, //3 ends
                time_tree_glob[0][2], time_tree_glob[1][2],
                time_tree_glob[2][2]/(double)p,
                time_tree_glob[0][3], time_tree_glob[1][3],
                time_tree_glob[2][3]/(double)p,
                time_tree_glob[1][3] + time_preproc, //4 ends
                dpengglob, tpengglob, fabs(tpengglob-dpengglob),
                fabs((tpengglob-dpengglob)/dpengglob),
                inferr, relinferr, n2err, reln2err); //5 ends
        fclose(fp);
    }
    
    
    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    
    free_vector(xT);
    free_vector(yT);
    free_vector(zT);
    free_vector(iT);
    free_vector(tenergy);

    if (rank == 0) {
        free_vector(denergyglob);
        free_vector(tenergyglob);
        if (pflag == 0) {
            free_vector(displs);
            free_vector(scounts);
        }
    }
    
    MPI_Finalize();
    return 0;
    
}
