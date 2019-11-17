#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>

#include "array.h"
#include "tools.h"
#include "kernels/kernels_direct.h"


static void directSummation(double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
                double *target_x, double *target_y, double *target_z, double *target_charge,
                int number_of_sources, int number_of_targets, double *potential, double *total_potential, double kernel_parameter,
                char *kernel_name);


int main(int argc, char **argv)
{
    int rank, numProcs, provided;
    
    MPI_Init(&argc, &argv);
    int ierr;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    /* runtime parameters */
    int numparsS, numparsT;

    /* arrays for coordinates, charges, energy of target particles */
    /* source particles */
    double *source_x = NULL; 
    double *source_y = NULL;
    double *source_z = NULL;
    double *source_charge = NULL;
    double *source_weight = NULL;

    /* target particles */
    double *target_x = NULL;
    double *target_y = NULL;
    double *target_z = NULL;
    double *target_charge = NULL;

    /* exact energy */
    double *potential = NULL;

    /* for potential energy calculation */
    double total_potential_local, total_potential_global;

    double kernel_parameter;

    // insert variables for date-time calculation?
    double time1, time2, time3, time4, time_direct, time_direct_tot, timeCommunicate, timeInteract;
    double time_direct_max, time_direct_min;
    double total_time_start, total_time_stop;
    
    /* input and output files */
    char *sampin1, *sampin2, *sampout, *sampdatout, *kernelName;
    FILE *fp;

    //local variables
    int i;
    int number_of_targets_local, maxparsTloc, location_in_global_target_array;
    int number_of_sources_local, maxparsSloc, location_in_global_source_array;
    double buf[5];
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       infile 1:  sources input file \n");                // "S10000.bin"
            printf("       infile 2:  targets input file \n");                // "T1000000.bin"
            printf("    pot outfile:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("    csv outfile:  human readable data output file \n");   // "out.txt"
            printf("       numparsS:  number of sources \n");                 // 10000
            printf("       numparsT:  number of targets \n");                 // 1000000
            printf("       pot type:  0--Coulomb\n");
            printf("                  1--screened Coulomb/Yukawa\n");
            printf("          kappa:  screened Coulomb parameter \n");        // 0.00


        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampout = argv[3];
    sampdatout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    kernel_parameter = atof(argv[7]);
    kernelName = argv[8];
    



    number_of_targets_local = (int)floor((double)numparsT/(double)numProcs);
    maxparsTloc = number_of_targets_local + (numparsT - (int)floor((double)numparsT/(double)numProcs) * numProcs);

	number_of_sources_local = (int)floor((double)numparsS/(double)numProcs);
	maxparsSloc = number_of_sources_local + (numparsS - (int)floor((double)numparsS/(double)numProcs) * numProcs);
    if (rank==0) printf("Num local T and local S: %i, %i.\n", number_of_sources_local, number_of_targets_local);
    
    if (rank == 0){
    	number_of_targets_local = maxparsTloc;
    	number_of_sources_local = maxparsSloc;
    }

    location_in_global_target_array = maxparsTloc + number_of_targets_local * (rank-1);
    location_in_global_source_array = maxparsSloc + number_of_sources_local * (rank-1);

    if (rank==0) printf("Max local T and local S: %i, %i.\n", maxparsTloc, maxparsSloc);

    double *S_local;
    make_vector(S_local, 5*maxparsSloc);
    source_x = &S_local[0*maxparsSloc];
    source_y = &S_local[1*maxparsSloc];
    source_z = &S_local[2*maxparsSloc];
    source_charge = &S_local[3*maxparsSloc];
    source_weight = &S_local[4*maxparsSloc];


    for (int i = 0; i < 5*maxparsSloc; i++) {
    	S_local[i] = 0.0;
    }


    double *T_local;
    make_vector(T_local, 4*number_of_targets_local);
    target_x = &T_local[0*number_of_targets_local];
    target_y = &T_local[1*number_of_targets_local];
    target_z = &T_local[2*number_of_targets_local];
    target_charge = &T_local[3*number_of_targets_local];


    make_vector(potential,number_of_targets_local);

    // Allocate and zero out receive buffers.
	int sendTo, recvFrom;
	double *xS_foreign1;
	double *yS_foreign1;
	double *zS_foreign1;
	double *qS_foreign1;
	double *wS_foreign1;
	double *S_foreign1;
	double *xS_foreign2;
	double *yS_foreign2;
	double *zS_foreign2;
	double *qS_foreign2;
	double *wS_foreign2;
	double *S_foreign2;

	make_vector(S_foreign1, 5*maxparsSloc);
//	memset(S_foreign1,0,5*maxparsSloc);
	xS_foreign1 = &S_foreign1[0*maxparsSloc];
	yS_foreign1 = &S_foreign1[1*maxparsSloc];
	zS_foreign1 = &S_foreign1[2*maxparsSloc];
	qS_foreign1 = &S_foreign1[3*maxparsSloc];
	wS_foreign1 = &S_foreign1[4*maxparsSloc];

	make_vector(S_foreign2, 5*maxparsSloc);
	xS_foreign2 = &S_foreign2[0*maxparsSloc];
	yS_foreign2 = &S_foreign2[1*maxparsSloc];
	zS_foreign2 = &S_foreign2[2*maxparsSloc];
	qS_foreign2 = &S_foreign2[3*maxparsSloc];
	wS_foreign2 = &S_foreign2[4*maxparsSloc];


    /* Reading in coordinates and charges for the source particles*/
    ierr = MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (ierr != MPI_SUCCESS) {
        fprintf(stderr,"FILE COULD NOT OPEN\n");
        exit(1);
    }

    MPI_File_seek(fpmpi, (MPI_Offset)location_in_global_source_array*5*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < number_of_sources_local; i++) {
        MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
        source_x[i] = buf[0];
        source_y[i] = buf[1];
        source_z[i] = buf[2];
        source_charge[i] = buf[3];
        source_weight[i] = buf[4];
    }
    MPI_File_close(&fpmpi);
    
//    MPI_Barrier(MPI_COMM_WORLD);

    /* Reading in coordinates for the targets */
    ierr = MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    if (ierr != MPI_SUCCESS) {
        fprintf(stderr,"FILE COULD NOT OPEN\n");
        exit(1);
    }

    MPI_File_seek(fpmpi, (MPI_Offset)location_in_global_target_array*4*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < number_of_targets_local; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        target_x[i] = buf[0];
        target_y[i] = buf[1];
        target_z[i] = buf[2];
        target_charge[i] = buf[3];
    }
    MPI_File_close(&fpmpi);

    total_potential_local = 0.0;
//    for (i = 0; i < maxparsTloc; i++) {
//		denergy[i] = 0.0;
//	}
 
#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    // Set up kernel
	if       (strcmp(kernelName,"coulomb")==0){
		if (rank==0) printf("Set kernel to coulombKernel.\n");
		for (i = 0; i < maxparsTloc; i++) {
			potential[i] = 0.0;
		}

	}else if (strcmp(kernelName,"yukawa")==0){
		if (rank==0) printf("Set kernel to yukawaKernel.\n");
		for (i = 0; i < maxparsTloc; i++) {
			potential[i] = 0.0;
		}

	}else if (strcmp(kernelName,"coulomb_SS")==0){
		if (rank==0) printf("Set kernel to coulombKernel_SS.\n");
		for (i = 0; i < maxparsTloc; i++) {
			potential[i] = 2.0*M_PI*kernel_parameter*kernel_parameter*target_charge[i];
		}

	}else if (strcmp(kernelName,"yukawa_SS")==0){
		if (rank==0) printf("Set kernel to yukawaKernel_SS.\n");
		for (i = 0; i < maxparsTloc; i++) {
			potential[i] = 4.0*M_PI*target_charge[i]/kernel_parameter/kernel_parameter;
		}

	}else{
		if (rank==0) printf("kernelName = %s.\n", kernelName);
		if (rank==0) printf("Invalid command line argument for kernelName... aborting.\n");
		return 1;
	}

    /* Interact with self */
    time1 = MPI_Wtime();
    total_time_start = MPI_Wtime();


	MPI_Request request1s;
	MPI_Request request1r;

	for (i=0; i < 5*maxparsSloc; i++) {
		S_foreign1[i]=0.0;
		S_foreign2[i]=0.0;
	}


	if (numProcs == 1) { // one 1 proc, won't enter into the round robin below.  So just interact with self here.
		directSummation(source_x, source_y, source_z, source_charge, source_weight, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local,
						potential, &total_potential_local, kernel_parameter, kernelName);
	} else {

		for (int procID = 1; procID < numProcs; procID++) {
			timeCommunicate = MPI_Wtime();

			// Send and receive source particles
			sendTo = (rank+procID)%numProcs;
			recvFrom = (numProcs+rank-procID)%numProcs;

			if (procID%2==1) { // communicate w/ foreign1, compute on foreign2
				MPI_Isend(S_local, 5*maxparsSloc, MPI_DOUBLE, sendTo, 1, MPI_COMM_WORLD, &request1s);
				MPI_Irecv(S_foreign1, 5*maxparsSloc, MPI_DOUBLE, recvFrom, 1, MPI_COMM_WORLD, &request1r);


				if (procID==1) { // in first iteration of loop, interact with self.
					directSummation(source_x, source_y, source_z, source_charge, source_weight, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local,
											potential, &total_potential_local, kernel_parameter, kernelName);
				} else { // in subsequent iterations, interact with foreign
					directSummation(xS_foreign2, yS_foreign2, zS_foreign2, qS_foreign2, wS_foreign2, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local, potential, &total_potential_local, kernel_parameter,kernelName);
					for (i=0;i<5*maxparsSloc;i++){
							S_foreign2[i]=0.0;
						}
				}
			} else { // procId%2=0, communicate foreign2 and compute with foreign1
				MPI_Isend(S_local, 5*maxparsSloc, MPI_DOUBLE, sendTo, 1, MPI_COMM_WORLD, &request1s);
				MPI_Irecv(S_foreign2, 5*maxparsSloc, MPI_DOUBLE, recvFrom, 1, MPI_COMM_WORLD, &request1r);

				directSummation(xS_foreign1, yS_foreign1, zS_foreign1, qS_foreign1, wS_foreign1, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local, potential, &total_potential_local, kernel_parameter, kernelName);
				for (i=0;i<5*maxparsSloc;i++){
						S_foreign1[i]=0.0;
					}
			}

			time4 = MPI_Wtime();
			MPI_Wait(&request1r, &status);
			MPI_Wait(&request1s, &status);



		} // end round robin

		if ((numProcs-1)%2==1) { // in final loop, S_foreign1 was received but not yet computed with
			directSummation(xS_foreign1, yS_foreign1, zS_foreign1, qS_foreign1, wS_foreign1, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local,
											potential, &total_potential_local, kernel_parameter, kernelName);
		} else { // S_foreign2 is the one that needs to be computed with
			directSummation(xS_foreign2, yS_foreign2, zS_foreign2, qS_foreign2, wS_foreign2, target_x, target_y, target_z, target_charge, maxparsSloc, number_of_targets_local,
											potential, &total_potential_local, kernel_parameter, kernelName);
		}
	}

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    
    /* Reducing values to root process */
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&total_potential_local, &total_potential_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
//        } // end OMP PARALLEL REGION
    total_time_stop = MPI_Wtime();


    MPI_File_open(MPI_COMM_WORLD, sampout, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0) {
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        MPI_File_write(fpmpi, &time_direct_max, 1, MPI_DOUBLE, &status);
    }
    MPI_File_seek(fpmpi, (MPI_Offset)(location_in_global_target_array+1)*sizeof(double), MPI_SEEK_SET);
    MPI_File_write(fpmpi, potential, number_of_targets_local, MPI_DOUBLE, &status);
    MPI_File_close(&fpmpi);

    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                 Max direct time (s):  %f\n", time_direct_max);
        printf("                 Min direct time (s):  %f\n", time_direct_min);
        printf("                 Avg direct time (s):  %f\n", time_direct_tot/(double)numProcs);
        printf("   Total direct time (s) on %d procs:  %f\n\n", numProcs, time_direct_tot);
    
        /* Calculating value dpeng by summing all values in denergy */
        printf("		          Wallclock time (s): %f\n", total_time_stop-total_time_start);
        printf("		         Cumulative time (s): %f\n", (total_time_stop-total_time_start)*numProcs);
        printf("             Direct potential energy:  %.15f\n\n\n", total_potential_global);
    }

    
    if (rank == 0) {
        fp = fopen(sampdatout, "a");
        fprintf(fp, "%s, %s, %s, %d, %d, %f, %s, %d,"
                "%f, %f, %f, %e \n",
                sampin1, sampin2, sampout, numparsS, numparsT,
                kernel_parameter, kernelName, numProcs, time_direct_max, time_direct_min,
                time_direct_tot/(double)numProcs, total_potential_global);
        fclose(fp);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) printf("Finished writing files.\n");
    free_vector(T_local);
    free_vector(S_local);


    free_vector(potential);

	free_vector(S_foreign1);
	free_vector(S_foreign2);
    MPI_Barrier(MPI_COMM_WORLD);
	if (rank==0) printf("Finished freeing memory.\n");
    MPI_Finalize();
    return 0;
}




static void directSummation(double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
                double *target_x, double *target_y, double *target_z, double *target_charge,
                int number_of_sources, int number_of_targets, double *potential, double *total_potential, double kernel_parameter,
                char *kernel_name)
{

#ifdef OPENACC_ENABLED
    #pragma acc data copyin (source_x[0:number_of_sources], source_y[0:number_of_sources], source_z[0:number_of_sources], \
                             source_charge[0:number_of_sources], source_weight[0:number_of_sources], source_x[0:number_of_targets], \
                             source_y[0:number_of_targets], source_z[0:number_of_targets], source_charge[0:number_of_targets])
    {
#endif



/***************************************/
/************* Coulomb *****************/
/***************************************/

    if (strcmp(kernel_name, "coulomb") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += coulombKernel(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/*************** Yukawa ****************/
/***************************************/

    } else if (strcmp(kernel_name, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += yukawaKernel(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Coulomb with SS ************/
/***************************************/

    } else if (strcmp(kernel_name, "coulomb_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += coulombKernel_SS(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Yukawa with SS *************/
/***************************************/

    } else if (strcmp(kernel_name, "yukawa_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += yukawaKernel_SS(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif

    } // end kernel selection


#ifdef OPENACC_ENABLED
    } // end acc data region
#endif

    *total_potential = sum(potential, number_of_targets);

    return;
}





