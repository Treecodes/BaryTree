#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <omp.h>
#include <float.h>


#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numDevices, int numThreads);

/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    printf("rank = %i\n",rank);
    
    /* runtime parameters */
    int numparsS, numparsT, numDevices, numThreads;
    int pot_type;

    /* arrays for coordinates, charges, energy of target particles */
    /* source particles */
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

    /* for potential energy calculation */
    double dpeng, dpengglob;

    double kappa;

    // insert variables for date-time calculation?
    double time1, time2, time_direct, time_direct_tot;
    double time_direct_max, time_direct_min;
    
    /* input and output files */
    char *sampin1, *sampin2, *sampout, *sampdatout;
    FILE *fp;

    //local variables
    int i;
    int numparsTloc, maxparsTloc, globparsTloc;
    double buf[5];
    
    /* MPI Variables */
    MPI_File fpmpi;
    MPI_Status status;
    printf("Setup MPI file variables.\n");
//    omp_set_num_threads(4);

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
    {
        if (rank == 0)
        {
            printf("Input arguments: \n");
            printf("       sampin1:  sources input file \n");                // "S10000.bin"
            printf("       sampin2:  targets input file \n");                // "T1000000.bin"
            printf("       sampout:  direct calc potential output file \n"); // "ex_s4_t6.bin"
            printf("    sampdatout:  human readable data output file \n");   // "out.txt"
            printf("      numparsS:  number of sources \n");                 // 10000
            printf("      numparsT:  number of targets \n");                 // 1000000
            printf("         kappa:  screened Coulomb parameter \n");        // 0.00
            printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");   // 1
            printf("   num devices:  number of GPUs available \n");   // 1
            printf("   num threads:  number of OpenMP threads \n");   // 1
        }
        return 0;
    }
    
    sampin2 = argv[2];
    sampout = argv[3];
    sampdatout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    kappa = atof(argv[7]);
    pot_type = atoi(argv[8]);
    numDevices = atoi(argv[9]);
    numThreads = atoi(argv[10]);
    
    numparsTloc = (int)floor((double)numparsT/(double)p);
    maxparsTloc = numparsTloc + (numparsT - (int)floor((double)numparsT/(double)p) * p);
    
    if (rank == 0) numparsTloc = maxparsTloc;
    globparsTloc = maxparsTloc + numparsTloc * (rank-1);

    make_vector(xS,numparsS);
    make_vector(yS,numparsS);
    make_vector(zS,numparsS);
    make_vector(qS,numparsS);
    make_vector(wS,numparsS);

    make_vector(xT,numparsTloc);
    make_vector(yT,numparsTloc);
    make_vector(zT,numparsTloc);
    make_vector(qT,numparsTloc);

    make_vector(denergy,numparsTloc);

    /* Reading in coordinates and charges for the source particles*/
    MPI_File_open(MPI_COMM_WORLD, sampin1, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
    for (i = 0; i < numparsS; i++) {
        MPI_File_read(fpmpi, buf, 5, MPI_DOUBLE, &status);
        xS[i] = buf[0];
        yS[i] = buf[1];
        zS[i] = buf[2];
        qS[i] = buf[3];
        wS[i] = buf[4];
    }
    MPI_File_close(&fpmpi);
    
    /* Reading in coordinates for the targets */
    MPI_File_open(MPI_COMM_WORLD, sampin2, MPI_MODE_RDONLY, MPI_INFO_NULL, &fpmpi);
    MPI_File_seek(fpmpi, (MPI_Offset)globparsTloc*4*sizeof(double), MPI_SEEK_SET);
    for (i = 0; i < numparsTloc; i++) {
        MPI_File_read(fpmpi, buf, 4, MPI_DOUBLE, &status);
        xT[i] = buf[0];
        yT[i] = buf[1];
        zT[i] = buf[2];
        qT[i] = buf[3];
    }
    MPI_File_close(&fpmpi);

    printf("Allocated sources and targets particles.\n");
	double *originalWeights, *rho;
    make_vector(originalWeights, numparsS);
	make_vector(rho, numparsS);
    for (int i = 0; i < numparsT; ++i) {
		originalWeights[i] = wS[i];
		rho[i] = qS[i];
    }
	printf("originalWeights filled.  Starting timer.\n");

	// Initialize all GPUs
	if (numDevices > 0) {
		#pragma omp parallel num_threads(numDevices)
        {
			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
			acc_init(acc_get_device_type());
        }
	}

    /* Calling main treecode subroutine to calculate approximate energy */
    time1 = MPI_Wtime();
    direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, numparsS, numparsTloc,
                denergy, &dpeng, pot_type, kappa, numDevices, numThreads);

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    double hartreeEnergy = 0.0;
    for (int i = 0; i < numparsTloc; ++i) {
    	hartreeEnergy += denergy[i] * rho[i] * originalWeights[i] / 2.0;  // sum ( V_H * RHO * W )
    }
    printf("Hartree energy = %1.12f\n", hartreeEnergy);

    
    /* Reducing values to root process */
    MPI_Reduce(&time_direct, &time_direct_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_tot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&dpeng, &dpengglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_File_open(MPI_COMM_WORLD, sampout, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fpmpi);
    if (rank == 0)
    {
        MPI_File_seek(fpmpi, (MPI_Offset)0, MPI_SEEK_SET);
        MPI_File_write(fpmpi, &time_direct_max, 1, MPI_DOUBLE, &status);
    }
    MPI_File_seek(fpmpi, (MPI_Offset)(globparsTloc+1)*sizeof(double), MPI_SEEK_SET);
    MPI_File_write(fpmpi, denergy, numparsTloc, MPI_DOUBLE, &status);
    MPI_File_close(&fpmpi);

    if (rank == 0)
    {
        /* Printing direct and treecode time calculations: */
        printf("                 Max direct time (s):  %f\n", time_direct_max);
        printf("                 Min direct time (s):  %f\n", time_direct_min);
        printf("                 Avg direct time (s):  %f\n", time_direct_tot/(double)p);
        printf("   Total direct time (s) on %d procs:  %f\n\n", p, time_direct_tot);
    
        /* Calculating value dpeng by summing all values in denergy */
        printf("             Direct potential energy:  %.15f\n\n\n", dpengglob);
        printf("			 Hartree energy = %1.12f\n", hartreeEnergy);

    }

    if (rank == 0) {
        fp = fopen(sampdatout, "a");
        fprintf(fp, "%s \t %s \t %s \t %d \t %d \t %f \t %d \t %d \t"
                "%f \t %f \t %f \t %e \n",
                sampin1, sampin2, sampout, numparsS, numparsT,
                kappa, pot_type, p, time_direct_max, time_direct_min,  
                time_direct_tot/(double)p, dpengglob);
        fclose(fp);
    }

    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    free_vector(xT);
    free_vector(yT);
    free_vector(zT);
    free_vector(denergy);

    MPI_Finalize();
    return 0;
}


void direct_eng( double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numDevices, int numThreads)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, qi, teng, rad;

    if (numDevices > 0) {
	    #pragma omp parallel num_threads(numDevices)
        {
            acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
	        #pragma acc data copyin ( xS [ 0 : numparsS ] , yS [ 0 : numparsS ] , zS [ 0 : numparsS ] , qS [ 0 : numparsS ] , wS [ 0 : numparsS ] , \
			xT [ 0 : numparsT ] , yT [ 0 : numparsT ] , zT [ 0 : numparsT ] , qT [ 0 : numparsT ] )
			{

			int this_thread = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
			if (this_thread == 0) printf("numDevices: %i\n", numDevices);
			if (this_thread == 0) printf("num_threads: %i\n", num_threads);

			int targetStart, targetEnd;
			#pragma omp for schedule(guided)
			for (int deviceNumber = 0; deviceNumber < num_threads; ++deviceNumber) {

				targetStart = deviceNumber*numparsT/num_threads;
				targetEnd = (deviceNumber+1)*numparsT/num_threads;
				printf("Device number: %i\n",deviceNumber);
				printf("Target start: %i\n", targetStart);
				printf("Target end: %i\n", targetEnd);
                
			    if (pot_type == 0 || pot_type == 4) { // Coulomb with singularity skipping.
					#pragma acc kernels
					{
					#pragma acc loop independent
					for (int i = targetStart; i < targetEnd; ++i) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        teng = 0.0;
                        #pragma acc loop independent
                        for (int j = 0; j < numparsS; ++j) {
                            tx = xi - xS[j];
                            ty = yi - yS[j];
                            tz = zi - zS[j];
                            rad = sqrt(tx*tx + ty*ty + tz*tz);
                            if (rad > 1e-14) {
                                teng = teng + qS[j]*wS[j] / rad;
                            }
                        }
                        denergy[i] =  teng;
					}
				    }

			    } else if (pot_type == 1 || pot_type == 5) {
	                #pragma acc kernels
				    {
	                #pragma acc loop independent
                    for (i = targetStart; i < targetEnd; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        teng = 0.0;

                        for (j = 0; j < numparsS; j++) {
                            tx = xi - xS[j];
                            ty = yi - yS[j];
                            tz = zi - zS[j];
                            rad = sqrt(tx*tx + ty*ty + tz*tz);
                            if (rad > 1e-14) {
                                teng += qS[j] * wS[j] * exp(-kappa * rad) / rad;
                            }
                        }
                        denergy[i] = teng;
                    }
			        }

			    } else if (pot_type == 2 || pot_type == 6) {
				    double kappaSq = kappa * kappa;
	                #pragma acc kernels
				    {
	                #pragma acc loop independent
					for (i = targetStart; i < targetEnd; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        qi = qT[i];
                        teng = 2 * M_PI * kappaSq * qi;  // 2pi alpha^2*f_t for SS scheme exp(-r^2/alpha^2)

                        for (j = 0; j < numparsS; j++) {
                            tx = xi - xS[j];
                            ty = yi - yS[j];
                            tz = zi - zS[j];
                            rad = sqrt(tx*tx + ty*ty + tz*tz);
                            if (rad > 1e-14) {
                                teng += (qS[j] - qi * exp(-rad*rad/kappaSq)) * wS[j] / rad;
                            }
                        }
                        denergy[i] = teng;
					}
			        }

			    } else if (pot_type == 3 || pot_type == 7) {
	                #pragma acc kernels
				    {
	                #pragma acc loop independent
				    for (i = targetStart; i < targetEnd; i++) {
						xi = xT[i];
						yi = yT[i];
						zi = zT[i];
						qi = qT[i];
						teng = 4 * M_PI * qi / kappa / kappa;  // 4pi*f_t/k^2

						for (j = 0; j < numparsS; j++) {
                            tx = xi - xS[j];
                            ty = yi - yS[j];
                            tz = zi - zS[j];
                            rad = sqrt(tx*tx + ty*ty + tz*tz);
                            if (rad>1e-14){
                                teng += ( qS[j] - qi) *wS[j] * exp(-kappa * rad) / rad;
                            }
						}
						denergy[i] = teng;
				    }
			        }

			    } // end pot
			} // end loop over device number
			} // end acc parallel
        } // end omp parallel
        
    } else { //numDevices=0, parallelize with openMP
        printf("Using OpenMP only, no GPU devices.\n");
        #pragma omp parallel num_threads(numThreads)
        {
		    int this_thread = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
		    if (this_thread == 0) printf("numDevices: %i\n", numDevices);
		    if (this_thread == 0) printf("num_threads: %i\n", num_threads);

            if (pot_type == 0 || pot_type == 4) { // Coulomb with singularity skipping.
        	    printf("Entering pot_type==0 or 4 section.\n");
			    #pragma omp for private(j,xi,yi,zi,tx,ty,tz,rad,teng) schedule(guided)
			    for (i = 0; i < numparsT; i++) {
                    xi = xT[i];
                    yi = yT[i];
                    zi = zT[i];
                    teng = 0.0;
                    for (j = 0; j < numparsS; j++) {
                        tx = xi - xS[j];
                        ty = yi - yS[j];
                        tz = zi - zS[j];
                        rad = sqrt(tx*tx + ty*ty + tz*tz);
                        if (rad > 1e-14) {
                            teng = teng + qS[j] * wS[j] / rad;
                        }
                    }
                    denergy[i] =  teng;
				}

            } else if (pot_type == 1 || pot_type == 5) {
				#pragma omp for private(j,xi,yi,zi,tx,ty,tz,rad,teng)
        		for (i = 0; i < numparsT; i++) {
                    xi = xT[i];
                    yi = yT[i];
                    zi = zT[i];
                    teng = 0.0;

                    for (j = 0; j < numparsS; j++) {
                        tx = xi - xS[j];
                        ty = yi - yS[j];
                        tz = zi - zS[j];
                        rad = sqrt(tx*tx + ty*ty + tz*tz);
                        if (rad > 1e-14) {
                            teng += qS[j] * wS[j] * exp(-kappa * rad) / rad;
                        }
                    }
                    denergy[i] = teng;
                }

            } else if (pot_type == 2 || pot_type == 6) {
        	    double kappaSq = kappa*kappa;
			    #pragma omp for private(j,xi,yi,zi,qi,tx,ty,tz,rad,teng) schedule(guided)
			    for (i = 0; i < numparsT; i++) {
				    xi = xT[i];
				    yi = yT[i];
				    zi = zT[i];
				    qi = qT[i];
				    teng = 2.0 * M_PI * kappaSq * qi;  // 2pi alpha^2*f_t for SS scheme exp(-r^2/alpha^2)

				    for (j = 0; j < numparsS; j++) {
                        tx = xi - xS[j];
                        ty = yi - yS[j];
                        tz = zi - zS[j];
                        rad = sqrt(tx*tx + ty*ty + tz*tz);
                        if (rad > DBL_MIN) {
                            teng += (qS[j] - qi * exp(-rad*rad/kappaSq)) * wS[j] / rad;
                        }
				    }
				    denergy[i] = teng;
                }

            } else if (pot_type == 3 || pot_type == 7) {
			    #pragma omp for private(j,xi,yi,zi,tx,ty,tz,rad,teng)
			    for (i = 0; i < numparsT; i++) {
                    xi = xT[i];
                    yi = yT[i];
                    zi = zT[i];
                    qi = qT[i];
                    teng = 4 * M_PI * qi / kappa / kappa;  // 4pi*f_t/k^2

                    for (j = 0; j < numparsS; j++) {
                        tx = xi - xS[j];
                        ty = yi - yS[j];
                        tz = zi - zS[j];
                        rad = sqrt(tx*tx + ty*ty + tz*tz);
                        if (rad > 1e-14) {
                            teng += (qS[j] - qi) * wS[j] * exp(-kappa * rad) / rad;
                        }
                    }
                    denergy[i] = teng;
			    }
            } // end pot
        } // end omp parallel
	} // end numDevices=0 region

    *dpeng = sum(denergy, numparsT);

    return;
}
