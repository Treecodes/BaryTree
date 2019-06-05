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
                int pot_type, double kappa);

/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{

    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* runtime parameters */
    int numparsS, numparsT;
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


    /* Calling main treecode subroutine to calculate approximate energy */

    time1 = MPI_Wtime();
    direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, numparsS, numparsTloc,
                denergy, &dpeng, pot_type, kappa);

    time2 = MPI_Wtime();
    time_direct = time2-time1;

    
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


//void direct_eng( double *xS,  double *yS,  double *zS,  double *qS,  double *wS,
//		 double *xT,  double *yT,  double *zT,  double *qT,
//                int numparsS, int numparsT, double *denergy, double *dpeng,
//                int pot_type, double kappa)
//{
//        /* local variables */
//        int i, j;
//        double tx, ty, tz, xi, yi, zi, qi, teng, rad;
//
//
//
//        if (pot_type == 0) {
////#pragma omp parallel for private(xi,yi,zi,teng,j,rad,tx,ty,tz)
//
//        	for (i = 0; i < numparsT; i++) {
//                        xi = xT[i];
//                        yi = yT[i];
//                        zi = zT[i];
//                        teng = 0.0;
//
//                        for (j = 0; j < numparsS; j++) {
//                                tx = xi - xS[j];
//                                ty = yi - yS[j];
//                                tz = zi - zS[j];
//                                rad = sqrt(tx*tx + ty*ty + tz*tz);
//                                if (rad>DBL_MIN){
//                                	teng = teng + qS[j]*wS[j] / rad;
//                                }
//                        }
//                        denergy[i] = teng;
//
//        	}
//
//        } else if (pot_type == 1) {
////#pragma omp parallel for private(xi,yi,zi,teng,j,rad,tx,ty,tz)
//                for (i = 0; i < numparsT; i++) {
//                        xi = xT[i];
//                        yi = yT[i];
//                        zi = zT[i];
//                        teng = 0.0;
//
//                        for (j = 0; j < numparsS; j++) {
//                                tx = xi - xS[j];
//                                ty = yi - yS[j];
//                                tz = zi - zS[j];
//                                rad = sqrt(tx*tx + ty*ty + tz*tz);
//                                if (rad>DBL_MIN){
//                                	teng = teng + qS[j]*wS[j] * exp(-kappa * rad) / rad;
//                                }
//                        }
//                        denergy[i] = teng;
//                }
//        } else if (pot_type == 3) {
//        	//#pragma omp parallel for private(xi,yi,zi,teng,j,rad,tx,ty,tz)
//						for (i = 0; i < numparsT; i++) {
//								xi = xT[i];
//								yi = yT[i];
//								zi = zT[i];
//								qi = qT[i];
//								teng = 4*M_PI*qi/kappa/kappa;  // 4pi*f_t/k^2
//
//								for (j = 0; j < numparsS; j++) {
//										tx = xi - xS[j];
//										ty = yi - yS[j];
//										tz = zi - zS[j];
//										rad = sqrt(tx*tx + ty*ty + tz*tz);
//										if (rad>DBL_MIN){
//											teng = teng + ( qS[j] - qi) *wS[j] * exp(-kappa * rad) / rad;
//										}
//								}
//								denergy[i] = teng;
//						}
//				}
//
//
//
//
//        *dpeng = sum(denergy, numparsT);
//
//        return;

//void direct_eng(__restrict__ double *xS, __restrict__ double *yS, __restrict__ double *zS, __restrict__ double *qS, __restrict__ double *wS,
//		__restrict__ double *xT, __restrict__ double *yT, __restrict__ double *zT, __restrict__ double *qT,
//                int numparsS, int numparsT, double *denergy, double *dpeng,
//                int pot_type, double kappa)
void direct_eng( double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa)
{
        /* local variables */
        int i, j;
        double tx, ty, tz, xi, yi, zi, qi, teng, rad;
//        int numDevices = acc_get_num_devices(acc_device_nvidia);

//        int idevtype = acc_get_device_type();
//        void acc_init ( idevtype );
#pragma omp parallel num_threads(acc_get_num_devices(acc_get_device_type()))
        {
        acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
//        int queue = 1;

#pragma acc data copyin ( xS [ 0 : numparsS ] , yS [ 0 : numparsS ] , zS [ 0 : numparsS ] , qS [ 0 : numparsS ] , wS [ 0 : numparsS ] , \
		xT [ 0 : numparsT ] , yT [ 0 : numparsT ] , zT [ 0 : numparsT ] , qT [ 0 : numparsT ] )
        {


		int numDevices = acc_get_num_devices(acc_get_device_type());
		printf("numDevices: %i\n", numDevices);
		int this_thread = omp_get_thread_num(), num_threads = omp_get_num_threads();
		printf("this_thread: %i\n", this_thread);
		printf("num_threads: %i\n", num_threads);

//		numDevices=1;
//		num_threads=1;

		int targetStart, targetEnd;
        if (pot_type == 0) {
			#pragma omp for schedule(static)
        	for (int deviceNumber=0;deviceNumber<num_threads;deviceNumber++){

        		targetStart = deviceNumber*numparsT/num_threads;
        		targetEnd = (deviceNumber+1)*numparsT/num_threads;

				#pragma acc kernels
				{
				#pragma acc loop independent
				for (i = targetStart; i < targetEnd; i++) {
							xi = xT[i];
							yi = yT[i];
							zi = zT[i];
							teng = 0.0;
							#pragma acc loop independent
							for (j = 0; j < numparsS; j++) {
									tx = xi - xS[j];
									ty = yi - yS[j];
									tz = zi - zS[j];
									rad = sqrt(tx*tx + ty*ty + tz*tz);
									if (rad>1e-14){
										teng = teng + qS[j]*wS[j] / rad;
									}
							}
							denergy[i] = teng;
					}
				}

        	}
//        	queue = (queue%2)+1;


        } else if (pot_type == 1) {
#pragma acc kernels
        	{
#pragma acc loop independent
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
                                if (rad>1e-14){
                                	teng = teng + qS[j]*wS[j] * exp(-kappa * rad) / rad;
                                }
                        }
                        denergy[i] = teng;
                }
        }

        } else if (pot_type == 2) {
        	double kappaSq = kappa*kappa;
#pragma acc kernels
        	{
#pragma acc loop independent
        		for (i = 0; i < numparsT; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        qi = qT[i];
						teng = 2*M_PI*kappaSq*qi;  // 2pi alpha^2*f_t for SS scheme exp(-r^2/alpha^2)

                        for (j = 0; j < numparsS; j++) {
                                tx = xi - xS[j];
                                ty = yi - yS[j];
                                tz = zi - zS[j];
                                rad = sqrt(tx*tx + ty*ty + tz*tz);
                                if (rad>1e-14){
                                	teng = teng + ( qS[j] - qi* exp(-rad*rad/kappaSq)) *wS[j]/ rad;
                                }
                        }
                        denergy[i] = teng;
                }
        }

        } else if (pot_type == 3) {
#pragma acc kernels
        	{
#pragma acc loop independent
			for (i = 0; i < numparsT; i++) {
					xi = xT[i];
					yi = yT[i];
					zi = zT[i];
					qi = qT[i];
					teng = 4*M_PI*qi/kappa/kappa;  // 4pi*f_t/k^2

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

        }
        }
        // Instead of just summing the final values, use their quadrature weights (assuming targets=sources)
//		  for (i=0;i<numparsT;i++){
//			  denergy[i] *= wS[i];
//		  }

        *dpeng = sum(denergy, numparsT);


} // end pragma omp parallel


        return;

}
