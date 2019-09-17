#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>
#include <float.h>


#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numThreads);

/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
    /* runtime parameters */
    int numparsS, numparsT, pot_type, numThreads;
    double kappa;

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
    double dpeng;

    // insert variables for date-time calculation?
    double time1, time2, time_direct;
    
    /* input and output files */
    char *sampin1, *sampin2, *sampout, *sampdatout;
    FILE *fp;

    //local variables
    double buf[5];
    

    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
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
        printf("   num threads:  number of OpenMP threads (or GPUs) \n");   // 1
        return 0;
    }
    
    sampin2 = argv[2];
    sampout = argv[3];
    sampdatout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    kappa = atof(argv[7]);
    pot_type = atoi(argv[8]);
    numThreads = atoi(argv[9]);

	// Initialize all GPUs
#ifdef OPENACC_ENABLED
	if (numThreads > 0) {
		#pragma omp parallel num_threads(numThreads)
        {
			acc_set_device_num(omp_get_thread_num(), acc_get_device_type());
			acc_init(acc_get_device_type());
        }
	} else {
        printf("Error! At least one GPU must be present for GPU version!\n");
        return 1;
    }
#endif
    
    make_vector(xS,numparsS);
    make_vector(yS,numparsS);
    make_vector(zS,numparsS);
    make_vector(qS,numparsS);
    make_vector(wS,numparsS);

    make_vector(xT,numparsT);
    make_vector(yT,numparsT);
    make_vector(zT,numparsT);
    make_vector(qT,numparsT);

    make_vector(denergy,numparsT);

    /* Reading in coordinates and charges for the source particles*/
    fp = fopen(sampin1, "rb");
    for (int i = 0; i < numparsS; i++) {
        fread(buf, sizeof(double), 5, fp);
        xS[i] = buf[0];
        yS[i] = buf[1];
        zS[i] = buf[2];
        qS[i] = buf[3];
        wS[i] = buf[4];
    }
    fclose(fp);
    
    /* Reading in coordinates for the targets */
    fp = fopen(sampin2, "rb");
    for (int i = 0; i < numparsT; i++) {
        fread(buf, sizeof(double), 4, fp);
        xT[i] = buf[0];
        yT[i] = buf[1];
        zT[i] = buf[2];
        qT[i] = buf[3];
    }
    fclose(fp);

    printf("Allocated sources and targets particles.\n");

    /* Calling main treecode subroutine to calculate approximate energy */
    time1 = omp_get_wtime();
    direct_eng(xS, yS, zS, qS, wS, xT, yT, zT, qT, numparsS, numparsT,
                denergy, &dpeng, pot_type, kappa, numThreads);

    time2 = omp_get_wtime();
    time_direct = time2-time1;

    
    fp = fopen(sampout, "wb");
    fwrite(&time_direct, sizeof(double), 1, fp);
    fwrite(denergy, sizeof(double), numparsT, fp);
    fclose(fp);

    /* Printing direct and treecode time calculations: */
    printf("                     Direct time (s):  %f\n\n", time_direct);
    
    /* Calculating value dpeng by summing all values in denergy */
    printf("             Direct potential energy:  %.15f\n\n\n", dpeng);


    fp = fopen(sampdatout, "a");
    fprintf(fp, "%s, %s, %s, %d, %d, %f, %d, %f, %e \n",
            sampin1, sampin2, sampout, numparsS, numparsT,
            kappa, pot_type, time_direct, dpeng);
    fclose(fp);


    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);
    free_vector(xT);
    free_vector(yT);
    free_vector(zT);
    free_vector(denergy);

    return 0;
}


void direct_eng( double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa, int numThreads)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, qi, teng, rad;

#ifdef OPENACC_ENABLED
	    #pragma omp parallel num_threads(numThreads)
        {
            acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
	        #pragma acc data copyin ( xS [ 0 : numparsS ] , yS [ 0 : numparsS ] , zS [ 0 : numparsS ] , qS [ 0 : numparsS ] , wS [ 0 : numparsS ] , \
			xT [ 0 : numparsT ] , yT [ 0 : numparsT ] , zT [ 0 : numparsT ] , qT [ 0 : numparsT ] )
			{

			int this_thread = omp_get_thread_num();
            int num_threads = omp_get_num_threads();

			int targetStart, targetEnd;
			#pragma omp for schedule(guided)
			for (int deviceNumber = 0; deviceNumber < num_threads; ++deviceNumber) {

				targetStart = deviceNumber*numparsT/num_threads;
				targetEnd = (deviceNumber+1)*numparsT/num_threads;
                
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
        
#else
        printf("Using OpenMP only, no GPU devices.\n");
        #pragma omp parallel num_threads(numThreads)
        {
		    int this_thread = omp_get_thread_num();
            int num_threads = omp_get_num_threads();
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
#endif

    *dpeng = sum(denergy, numparsT);

    return;
}
