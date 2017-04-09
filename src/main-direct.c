#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

#include "array.h"
#include "tools.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xT, double *yT, double *zT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa);

/* The treedriver routine in Fortran */
int main()
{

        /* runtime parameters */
        int numparsS, numparsT;
        int pot_type;

        double temp;

        /* arrays for coordinates, charges, energy of target particles */
        double *xS, *yS, *zS, *qS;  /* source particles */
        double *xT, *yT, *zT;       /* target particles */
        double *denergy;  /* exact energy, treecode energy */

        /* for potential energy calculation */
        double dpeng;

        double kappa;


        // insert variables for date-time calculation?
        double time_direct;
        time_t time1, time2;


        /* input and output files */
        char sampin1[20], sampin2[20], sampin3[20], sampout[20];
        //char *sampin1, *sampin2, *sampout;
        FILE *fp;


        //local variables (these decs may be unnecessary in C)
        int i, j, err, ii;


        printf("Enter name of input file 1: \n");
        scanf("%s", sampin1);

        //sampin1 = "S10000.txt";

        printf("Enter name of input file 2: \n");
        scanf("%s", sampin2);

        //sampin2 = "T10000.txt";

        printf("Enter name of output file: \n");
        scanf("%s", sampout);

        //sampout = "ex_s4_t4.txt";

        printf("Enter particle number for source and target: \n");
        scanf("%d %d", &numparsS, &numparsT);

        //numparsS = 10000;
        //numparsT = 10000;

        printf("Enter kappa value: \n");
        scanf("%lf", &kappa);

        //kappa = 0.1;

        printf("Enter potential type (0 Coulomb; 1 Yukawa) \n");
        scanf("%d", &pot_type);

        //pot_type = 1;


        make_vector(xS,numparsS);
        make_vector(yS,numparsS);
        make_vector(zS,numparsS);
        make_vector(qS,numparsS);
        
        printf("xS, yS, zS, qS allocated\n");

        make_vector(xT,numparsT);
        make_vector(yT,numparsT);
        make_vector(zT,numparsT);

        printf("xT, yT, zT allocated\n");

        make_vector(denergy,numparsT);

        printf("denergy allocated\n\n");


        /* Reading in coordinates and charges for the source particles*/
        fp = fopen(sampin1, "r");
        if (!fp) {
                perror("File opening failed");
                return EXIT_FAILURE;
        }

        for (int i = 0; i < numparsS; i++)
        {
                fscanf(fp, " %lf %lf %lf %lf", &xS[i], &yS[i], &zS[i], &qS[i]); 
        }

        fclose(fp);

        /* Reading in coordinates for the targets */
        
        fp = fopen(sampin2, "r");
        if (!fp) {
                perror("File opening failed");
                return EXIT_FAILURE;
        }

        for (int i = 0; i < numparsT; i++)
                fscanf(fp, " %lf %lf %lf", &xT[i], &yT[i], &zT[i]);

        fclose(fp);


        /* Calling main treecode subroutine to calculate approximate energy */

        fp = fopen(sampout, "a");

        time1 = time(NULL);

        direct_eng(xS, yS, zS, qS, xT, yT, zT, numparsS, numparsT,
                   denergy, &dpeng, pot_type, kappa);

        time2 = time(NULL);
        time_direct = difftime(time2, time1);


        fprintf(fp, " %lf\n", time_direct);
        for (int i = 0; i < numparsT; i++) 
                fprintf(fp, " %d %lf\n", i, denergy[i]); 
        fclose(fp); 


        /* Printing direct and treecode time calculations: */
        printf("          Direct time (s):  %f\n", time_direct);


        /* Calculating value dpeng by summing all values in denergy */ 
        printf("  Direct potential energy:  %f\n", dpeng);


        return 0;

}


void direct_eng(double *xS, double *yS, double *zS, double *qS, 
                double *xT, double *yT, double *zT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa)
{
        /* local variables */
        int i, j;
        double tx, ty, tz, xi, yi, zi, teng, rad;

        *dpeng = 0.0;

        if (pot_type == 0) {
                for (i = 0; i < numparsT; i++) {
                        xi = xT[i];
                        yi = yT[i];
                        zi = zT[i];
                        teng = 0.0;

                        for (j = 0; j < numparsS; j++) {
                                tx = xi - xS[j];
                                ty = yi - yS[j];
                                tz = zi - zS[j];
                                teng = teng + qS[j] / sqrt(tx*tx + ty*ty + tz*tz);
                        }
                        denergy[i] = teng;
                }

        } else if (pot_type == 1) {
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
                                teng = teng + qS[j] * exp(-kappa * rad) / rad;
                        }
                        denergy[i] = teng;
                }
        }


        *dpeng = sum(denergy, numparsT);


        return;

}
