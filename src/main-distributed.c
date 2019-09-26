#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "particles.h"
#include "sort.h"


/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
//	printf("Entering main.c\n");
    int rank, numProcs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    /* runtime parameters */
    int numparsS, numparsT, order;
    int maxparnode;
    int pot_type, tree_type;
    int pflag, sflag, dflag, gflag = 0;
    int batch_size;
    int numDevices;
    int numThreads;

    double theta, temp;
    double kappa;

    /* source particles */
    struct particles *sources = NULL;
    
    /* target particles */
    struct particles *targets = NULL;

    /* exact energy, treecode energy */
//    double *denergyglob = NULL;
//    double *tenergyglob = NULL;
    double *tenergy = NULL;

    /* for potential energy calculation */
    double dpeng = 0;
    double tpeng = 0;
    double dpengglob = 0;
    double tpengglob = 0;

    /* insert variables for date-time calculation? */
    double time_direct, time_tree[4], time_preproc, time_treedriver;
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
    int numparsTloc, maxparsTloc, numparsSloc, maxparsSloc;
    double buf[5];
    
    double *xS = NULL;
	double *yS = NULL;
	double *zS = NULL;
	double *qS = NULL;
	double *wS = NULL;

	/* target particles */
	double *xT = NULL;
	double *yT = NULL;
	double *zT = NULL;
	double *qT = NULL;

        /* exact energy */
        double *denergy = NULL;

    int *displs = NULL;
    int *scounts = NULL;
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;


    /* Executable statements begin here */
//    printf("Reading in arguments.c\n");
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
            printf("         kappa:  screened Coulomb parameter \n");                // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");           // 1
            printf("         pflag:  distribute 0--targets, 1--sources \n");         // 0
            printf("         sflag:  on distr 0--sort, 1--no sort \n");              // 0
            printf("         dflag:  if sorted, direction 0--x, 1--y, 2--z \n");     // 0
            printf("         batch_size:  size of target batches \n");     // 0
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
    kappa = atof(argv[11]);
    pot_type = atoi(argv[12]);
    pflag = atoi(argv[13]);
    sflag = atoi(argv[14]);
    dflag = atoi(argv[15]);
    batch_size = atoi(argv[16]);
    numDevices = atoi(argv[17]);
    numThreads = atoi(argv[18]);

//    printf("Read in arguments.c\n");


    numparsTloc = (int)floor((double)numparsT/(double)numProcs);
	maxparsTloc = numparsTloc + (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);

	numparsSloc = (int)floor((double)numparsS/(double)numProcs);
	maxparsSloc = numparsSloc + (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
	if (rank==0) printf("Num local T and local S: %i, %i.\n", numparsSloc, numparsTloc);

	if (rank == 0){
		numparsTloc = maxparsTloc;
		numparsSloc = maxparsSloc;
	}
	int globparsTloc, globparsSloc;
	globparsTloc = maxparsTloc + numparsTloc * (rank-1);
	globparsSloc = maxparsSloc + numparsSloc * (rank-1);


    make_vector(denergy,numparsTloc);
    
    time1 = MPI_Wtime();
    
    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));

    sources->num = numparsSloc;
	make_vector(sources->x, numparsSloc);
	make_vector(sources->y, numparsSloc);
	make_vector(sources->z, numparsSloc);
	make_vector(sources->q, numparsSloc);
	make_vector(sources->w, numparsSloc);
//	printf("Made source vectors.c\n");

	targets->num = numparsTloc;
	make_vector(targets->x, numparsTloc);
	make_vector(targets->y, numparsTloc);
	make_vector(targets->z, numparsTloc);
	make_vector(targets->q, numparsTloc);
	make_vector(targets->order, numparsTloc);
	make_vector(tenergy, numparsTloc);
//	printf("Made target vectors.c\n");

	/* Reading in coordinates and charges for the source particles*/
	int offset = rank*(numparsSloc);
//	printf("Proc %i reading in from %i to %i.\n", rank, offset, offset+numparsSloc);
	MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset) offset*5*sizeof(double), MPI_SEEK_SET);
	for (i = 0; i < numparsSloc; i++) {
		MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
		sources->x[i] = buf[0];
		sources->y[i] = buf[1];
		sources->z[i] = buf[2];
		sources->q[i] = buf[3];
		sources->w[i] = buf[4];
	}
	MPI_File_close(&fpmpi);
//	printf("Read in sources.\n");

	/* Reading in coordinates for target particles*/
	MPI_File_open(MPI_COMM_SELF, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset)(rank*(numparsTloc)*4*sizeof(double)), MPI_SEEK_SET);
	for (i = 0; i < numparsTloc; i++) {
		MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
		targets->x[i] = buf[0];
		targets->y[i] = buf[1];
		targets->z[i] = buf[2];
		targets->q[i] = buf[3];
		targets->order[i] = i;
	}
	MPI_File_close(&fpmpi);
//	printf("Read in targets.\n");

	make_vector(tenergy, numparsTloc);
	make_vector(denergy, numparsTloc);

	// reading in file containing direct sum results.
	MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset)(rank*(numparsTloc)*1*sizeof(double)), MPI_SEEK_SET);
	MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
	MPI_File_read(fpmpi, denergy, numparsTloc, MPI_DOUBLE, &status);
	MPI_File_close(&fpmpi);
//	printf("Did MPI file stuff.\n");



//    printf("Allocated sources and targets particles.\n");
//    printf("PFLAG = %i\n",pflag);
//
//    int ierr;
//    /* Reading in coordinates and charges for the source particles*/
//	ierr = MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
//	if (ierr != MPI_SUCCESS) {
//		fprintf(stderr,"FILE COULD NOT OPEN\n");
//		exit(1);
//	}
//
//	MPI_File_seek(fpmpi, (MPI_Offset)globparsSloc*5*sizeof(double), MPI_SEEK_SET);
//	for (i = 0; i < numparsSloc; i++) {
//		MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
//		xS[i] = buf[0];
//		yS[i] = buf[1];
//		zS[i] = buf[2];
//		qS[i] = buf[3];
//		wS[i] = buf[4];
//	}
//	MPI_File_close(&fpmpi);

//    MPI_Barrier(MPI_COMM_WORLD);

	/* Reading in coordinates for the targets */
//	ierr = MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
//	if (ierr != MPI_SUCCESS) {
//		fprintf(stderr,"FILE COULD NOT OPEN\n");
//		exit(1);
//	}
//
//	MPI_File_seek(fpmpi, (MPI_Offset)globparsTloc*4*sizeof(double), MPI_SEEK_SET);
//	for (i = 0; i < numparsTloc; i++) {
//		MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
//		xT[i] = buf[0];
//		yT[i] = buf[1];
//		zT[i] = buf[2];
//		qT[i] = buf[3];
//	}
//	MPI_File_close(&fpmpi);
        
        
    // Initialize all GPUs
//	if (numDevices>0){
//		#pragma omp parallel num_threads(numDevices)
//			{
//			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
//
//
//
//			acc_init(acc_get_device_type());
//			}
//	}


	#pragma acc set device_num(rank%numDevices) device_type(acc_device_nvidia)
	#pragma acc init device_type(acc_device_nvidia)



    time2 = MPI_Wtime();
    time_preproc = time2 - time1;
//    printf("Setup complete, calling treedriver...\n");
//    printf("numThreads: %i\n", numThreads);
//    fflush(stdout);
    /* Calling main treecode subroutine to calculate approximate energy */

    MPI_Barrier(MPI_COMM_WORLD);
    time1 = MPI_Wtime();
    treedriver(sources, targets,
               order, theta, maxparnode, batch_size,
               pot_type, kappa, tree_type,
               tenergy, &tpeng, time_tree, numDevices, numThreads);

    MPI_Barrier(MPI_COMM_WORLD);
    time2 = MPI_Wtime();
    time_treedriver = time2 - time1;

    
//    printf("About to do reductions.\n");
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    printf("Completed time reductions.\n");
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    printf("Completed reductions.\n");

    dpeng = sum(denergy, numparsTloc);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    
    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n\n", time_direct);
        printf("              Pre-process time (s):  %f\n", time_preproc);
        printf("              Treedriver time (s):  %f\n", time_treedriver);
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
//    if (pflag == 0) {
//
//        if (rank == 0) {
//            scounts[0] = maxparsTloc;
//            displs[0] = 0;
//            for (i=1; i < numProcs; i++) {
//                scounts[i] = numparsTloc;
//                displs[i] = maxparsTloc + (i-1) * numparsTloc;
//            }
//        }
//
//        MPI_Gatherv(tenergy, numparsTloc, MPI_DOUBLE, tenergyglob,
//                    scounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//    } else if (pflag == 1) {
//        MPI_Reduce(tenergy, tenergyglob, numparsT,
//                   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
//    }



    // have each proc compute its L2 error.  Then reduce and print with rank 0
//    printf("Rank %i computing its errors.\n", rank);
    double glob_reln2_err;
    double glob_relinf_err;

    inferr = 0.0;
	relinferr = 0.0;
	n2err = 0.0;
	reln2err = 0.0;
	double x,y,z;

	for (j = 0; j < numparsTloc; j++) {
		temp = fabs(denergy[targets->order[j]] - tenergy[j]);

		if (temp >= inferr){
			inferr = temp;
			x = targets->x[targets->order[j]];
			y = targets->y[targets->order[j]];
			z = targets->z[targets->order[j]];
		}


		if (fabs(denergy[j]) >= relinferr)
			relinferr = fabs(denergy[j]);

		n2err = n2err + pow(denergy[targets->order[j]]
						  - tenergy[j], 2.0)*sources->w[j];
		reln2err = reln2err + pow(denergy[j], 2.0)*sources->w[j];
	}

	relinferr = inferr / relinferr;
	reln2err = sqrt(n2err / reln2err);
	n2err = sqrt(n2err);


//	printf("Computed errors on each rank.  Now reducing.\n");
	MPI_Reduce(&reln2err, &glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&relinferr, &glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {


//        printf("Absolute inf norm error in potential:  %e \n", inferr);
        printf("Relative inf norm error in potential:  %e \n\n", glob_relinf_err);
//        printf("  Absolute 2 norm error in potential:  %e \n", n2err);
        printf("  Relative 2 norm error in potential:  %e \n\n", glob_reln2_err);

//        printf("inf error occurring at %f, %f, %f \n\n", x, y,z);
    }
    
    
//    if (rank == 0) {
//        fp = fopen(sampout, "a");
//        fprintf(fp, "%s \t %s \t %s \t %d \t %d \t %f \t %d \t %d \t %d \t"
//                "%f \t %d \t %d \t %d \t"
//                "%d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t"
//                "%f \t %f \t %f \t %f \t %f \t %f \t %f \t"
//                "%e \t %e \t %e \t %e \t %e \t %e \t %e \t %e \n",
//                sampin1, sampin2, sampin3, numparsS, numparsT,
//                theta, order, tree_type, maxparnode,
//                kappa, pot_type, sflag, pflag, //2 ends
//                p, time_preproc,The
//                time_tree_glob[0][0], time_tree_glob[1][0],
//                time_tree_glob[2][0]/(double)p,
//                time_tree_glob[0][1], time_tree_glob[1][1],
//                time_tree_glob[2][1]/(double)p, //3 ends
//                time_tree_glob[0][2], time_tree_glob[1][2],
//                time_tree_glob[2][2]/(double)p,
//                time_tree_glob[0][3], time_tree_glob[1][3],
//                time_tree_glob[2][3]/(double)p,
//                time_tree_glob[1][3] + time_preproc, //4 ends
//                dpengglob, tpengglob, fabs(tpengglob-dpengglob),
//                fabs((tpengglob-dpengglob)/dpengglob),
//                inferr, relinferr, n2err, reln2err); //5 ends
//        fclose(fp);
//    }
    if (rank == 0){
            fp = fopen(sampout, "a");
            fprintf(fp, "%s,%s,%s,%d,%d,%f,%d,%d,%d,%d,%f,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%e,%e,%e,%e,%e,%e,%e,%e,%d\n",
                    sampin1, sampin2, sampin3, numparsS, numparsT,
                    theta, order, tree_type, maxparnode,batch_size,
                    kappa, pot_type, sflag, pflag, //2 ends
                    numProcs, time_preproc,
                    time_tree_glob[0][0], time_tree_glob[1][0],
                    time_tree_glob[2][0]/(double)numProcs,
                    time_tree_glob[0][1], time_tree_glob[1][1],
                    time_tree_glob[2][1]/(double)numProcs, //3 ends
                    time_tree_glob[0][2], time_tree_glob[1][2],
                    time_tree_glob[2][2]/(double)numProcs,
                    time_tree_glob[0][3], time_tree_glob[1][3],
                    time_tree_glob[2][3]/(double)numProcs,
                    time_tree_glob[1][3] + time_preproc, //4 ends
                    dpengglob, tpengglob, fabs(tpengglob-dpengglob),
                    fabs((tpengglob-dpengglob)/dpengglob),
                    inferr, relinferr, n2err, reln2err,numDevices); //5 ends
            fclose(fp);
        }
//    printf("Wrote to output file.\n");
    
    
    free_vector(sources->x);
    free_vector(sources->y);
    free_vector(sources->z);
    free_vector(sources->q);
    free_vector(sources->w);
//    free_vector(sources->order);

//    printf("Freed sources.\n");
    free_vector(targets->x);
    free_vector(targets->y);
    free_vector(targets->z);
    free_vector(targets->q);
//    free_vector(targets->w);
    free_vector(targets->order);
//    printf("Freed targets.\n");
    
    free(sources);
    free(targets);
    free_vector(tenergy);

//    if (rank == 0) {
//        free_vector(denergyglob);
//        free_vector(tenergyglob);
//        if (pflag == 0) {
//            free_vector(displs);
//            free_vector(scounts);
//        }
//    }


//    printf("Freed other vectors.\n");

    MPI_Finalize();
//    printf("Final.\n");
    return 0;
    
}
