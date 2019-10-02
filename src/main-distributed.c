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
    int maxparnode, batch_size, pot_type;
    double theta, kappa;

    struct particles *sources = NULL;
    struct particles *targets = NULL;

    double *tenergy = NULL;
    double *denergy = NULL;

    /* for potential energy calculation */
    double dpeng = 0;
    double tpeng = 0;
    double dpengglob = 0;
    double tpengglob = 0;

    /* variables for date-time calculation */
    double time_direct, time_tree[4], time_preproc, time_treedriver;
    double time_tree_glob[3][4];
    double time1, time2;

    /* input and output files */
    char *sampin1 = NULL;
    char *sampin2 = NULL;
    char *sampin3 = NULL;
    char *offset1 = NULL;
    char *offset2 = NULL;
    char *sampout = NULL;
    FILE *fp;

    double buf[5];
    
    double *xS = NULL;
	double *yS = NULL;
	double *zS = NULL;
	double *qS = NULL;
	double *wS = NULL;

	double *xT = NULL;
	double *yT = NULL;
	double *zT = NULL;
	double *qT = NULL;

    /* MPI Variables */
    int *displs = NULL;
    int *scounts = NULL;
    int mpi_err;
    MPI_File fpmpi;
    MPI_Status status;


    /* Executable statements begin here */
    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0) {
        if (rank == 0) {
            printf("Input arguments... \n");
            printf("        infile 1:  sources input file \n");
            printf("        infile 2:  targets input file \n");
            printf("        infile 3:  sources offset file \n");
            printf("        infile 4:  targets offset file \n");
            printf("        infile 5:  direct calc potential input file \n");
            printf("      csv output:  results summary to CSV file \n");
            printf("        numparsS:  number of sources \n");
            printf("        numparsT:  number of targets \n");
            printf("           theta:  multipole acceptance criterion \n");
            printf("           order:  number of Chebyshev interp. pts per Cartesian direction \n");
            printf("      maxparnode:  maximum particles in leaf \n");
            printf("      batch size:  maximum size of target batch \n");
            printf(" pot/approx type:  0--Coulomb, Lagrange approx.\n");
            printf("                   1--screened Coulomb/Yukawa, Lagrange approx.\n");
            printf("                   4--Coulomb, Hermite approx.\n");
            printf("                   5--screened Coulomb/Yukawa, Hermite approx.\n");
            printf("           kappa:  screened Coulomb parameter \n");
        }
        return 0;
    }
    
    sampin2 = argv[2];
    offset1 = argv[3];
    offset2 = argv[4];
    sampin3 = argv[5];
    sampout = argv[6];
    numparsS = atoi(argv[7]);
    numparsT = atoi(argv[8]);
    theta = atof(argv[9]);
    order = atoi(argv[10]);
    maxparnode = atoi(argv[11]);
    batch_size = atoi(argv[12]);
    pot_type = atoi(argv[13]);
    kappa = atof(argv[14]);
    
    int numparsSloc, numparsTloc, local_sources_offset, local_targets_offset;
    
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, offset1,
        MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        numparsSloc = (int)floor((double)numparsS/(double)numProcs);
        local_sources_offset = numparsSloc * rank;
        if (rank == numProcs-1)
            numparsSloc += (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
    } else {
        int *sources_offset;
        make_vector(sources_offset, numProcs);
        MPI_File_read(fpmpi, sources_offset, numProcs, MPI_INT, &status);
        if (rank == numProcs-1)
            numparsSloc = sources_offset[rank+1] - sources_offset[rank];
        else
            numparsSloc = numparsS - sources_offset[rank];
        local_sources_offset = sources_offset[rank];
        MPI_File_close(&fpmpi);
        free_vector(sources_offset);
    }
    
    if ((mpi_err = MPI_File_open(MPI_COMM_WORLD, offset1,
        MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi)) != MPI_SUCCESS) {
        numparsTloc = (int)floor((double)numparsT/(double)numProcs);
        local_targets_offset = numparsTloc * rank;
        if (rank == numProcs-1)
            numparsTloc += (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);
    } else {
        int *targets_offset;
        make_vector(targets_offset, numProcs);
        MPI_File_read(fpmpi, targets_offset, numProcs, MPI_INT, &status);
        if (rank == numProcs-1)
            numparsTloc = targets_offset[rank+1] - targets_offset[rank];
        else
            numparsTloc = numparsT - targets_offset[rank];
        local_targets_offset = targets_offset[rank];
        MPI_File_close(&fpmpi);
        free_vector(targets_offset);
    }

    
	int globparsTloc = numparsT;
	int globparsSloc = numparsS;
    
    time1 = MPI_Wtime();
    
    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));

    sources->num = numparsSloc;
	make_vector(sources->x, numparsSloc);
	make_vector(sources->y, numparsSloc);
	make_vector(sources->z, numparsSloc);
	make_vector(sources->q, numparsSloc);
	make_vector(sources->w, numparsSloc);

	targets->num = numparsTloc;
	make_vector(targets->x, numparsTloc);
	make_vector(targets->y, numparsTloc);
	make_vector(targets->z, numparsTloc);
	make_vector(targets->q, numparsTloc);
	make_vector(targets->order, numparsTloc);
	make_vector(tenergy, numparsTloc);

	/* Reading in coordinates and charges for the source particles*/
	MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset) (local_sources_offset * 5 * sizeof(double)), MPI_SEEK_SET);
	for (int i = 0; i < numparsSloc; ++i) {
		MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
		sources->x[i] = buf[0];
		sources->y[i] = buf[1];
		sources->z[i] = buf[2];
		sources->q[i] = buf[3];
		sources->w[i] = buf[4];
	}
	MPI_File_close(&fpmpi);

	/* Reading in coordinates for target particles*/
	MPI_File_open(MPI_COMM_SELF, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset) (local_targets_offset * 4 * sizeof(double)), MPI_SEEK_SET);
	for (int i = 0; i < numparsTloc; ++i) {
		MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
		targets->x[i] = buf[0];
		targets->y[i] = buf[1];
		targets->z[i] = buf[2];
		targets->q[i] = buf[3];
		targets->order[i] = i;
	}
	MPI_File_close(&fpmpi);

	make_vector(tenergy, numparsTloc);
	make_vector(denergy, numparsTloc);

	// reading in file containing direct sum results.
	MPI_File_open(MPI_COMM_SELF, sampin3, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
	MPI_File_seek(fpmpi, (MPI_Offset)(local_targets_offset * sizeof(double)), MPI_SEEK_SET);
	MPI_File_read(fpmpi, &time_direct, 1, MPI_DOUBLE, &status);
	MPI_File_read(fpmpi, denergy, numparsTloc, MPI_DOUBLE, &status);
	MPI_File_close(&fpmpi);


	#pragma acc set device_num(rank) device_type(acc_device_nvidia)
	#pragma acc init device_type(acc_device_nvidia)

    time2 = MPI_Wtime();
    time_preproc = time2 - time1;


    /* Calling main treecode subroutine to calculate approximate energy */

    MPI_Barrier(MPI_COMM_WORLD);
    time1 = MPI_Wtime();
    treedriver(sources, targets, order, theta, maxparnode, batch_size,
               pot_type, kappa, 1, tenergy, &tpeng, time_tree);

    MPI_Barrier(MPI_COMM_WORLD);
    time2 = MPI_Wtime();
    time_treedriver = time2 - time1;

    
    
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&tpeng, &tpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

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
        
        printf("      Min, Max total tree time (s):  %f, %f\n\n", time_tree_glob[0][3],
                                                                  time_tree_glob[1][3]);
        printf(" Preproc + Max total tree time (s):  %f \n\n", time_tree_glob[1][3] + time_preproc);
        

        printf("           Direct potential energy:  %f\n", dpengglob);
        printf("             Tree potential energy:  %f\n\n", tpengglob);
    
        printf("Absolute error for total potential:  %e\n",
               fabs(tpengglob-dpengglob));
        printf("Relative error for total potential:  %e\n\n",
               fabs((tpengglob-dpengglob)/dpengglob));
    }
    
    
    double glob_reln2_err, glob_relinf_err;
    double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
	double x, y, z, temp;

	for (int j = 0; j < numparsTloc; ++j) {
		temp = fabs(denergy[targets->order[j]] - tenergy[j]);

		if (temp >= inferr) {
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


	MPI_Reduce(&reln2err, &glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&relinferr, &glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("Relative inf norm error in potential:  %e \n\n", glob_relinf_err);
        printf("  Relative 2 norm error in potential:  %e \n\n", glob_reln2_err);
    }
    
    if (rank == 0) {
        fp = fopen(sampout, "a");
        fprintf(fp, "%s,%s,%s,%d,%d,%f,%d,%d,%d,%f,%d,%d,%f,"
                    "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%e,%e,%e,%e,%e,%e,%e,%e,%d\n",
            sampin1, sampin2, sampin3, numparsS, numparsT,
            theta, order, maxparnode, batch_size, kappa, pot_type, //2 ends
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
            inferr, relinferr, n2err, reln2err); //5 ends
        fclose(fp);
    }
    
    free_vector(sources->x);
    free_vector(sources->y);
    free_vector(sources->z);
    free_vector(sources->q);
    free_vector(sources->w);

    free_vector(targets->x);
    free_vector(targets->y);
    free_vector(targets->z);
    free_vector(targets->q);
    free_vector(targets->order);
    
    free(sources);
    free(targets);
    free_vector(tenergy);

    MPI_Finalize();

    return 0;
}
