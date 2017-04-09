#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"

/* The treedriver routine in Fortran */
int main()
{

        /* runtime parameters */
        int numparsS, numparsT, order, maxparnodeS, maxparnodeT;
        int shrinkS, shrinkT, treelevelS, treelevelT, iflagS, iflagT;
        int pot_type;

        double theta, temp;
        double kappa;

        /* arrays for coordinates, charges, energy of target particles */
        double *xS, *yS, *zS, *qS;  /* source particles */
        double *xT, *yT, *zT;       /* target particles */
        double *denergy, *tenergy;  /* exact energy, treecode energy */

        /* for potential energy calculation */
        double tpeng, dpeng;


        // insert variables for date-time calculation?
        double total_time, time_direct, time_tree;


        /* input and output files */
        char sampin1[20], sampin2[20], sampin3[20], *sampout;
        //char *sampin1, *sampin2, *sampin3, *sampout;
        FILE *fp;

        /* variables for error calculations */
        double inferr, relinferr, n2err, reln2err;

        //local variables (these decs may be unnecessary in C)
        int i, j, err, ii;


        printf("Enter name of input file 1: \n");
        scanf("%s", sampin1);

        //sampin1 = "S10000.txt";

        printf("Enter name of input file 2: \n");
        scanf("%s", sampin2);

        //sampin2 = "T1000000.txt";

        printf("Enter name of input file 3: \n");
        scanf("%s", sampin3);

        //sampin3 = "ex_s4_t6.txt";

        //printf("Enter name of output file: \n");
        //scanf("%s", sampout);

        sampout = "out.txt";

        printf("Enter particle number for source and target: \n");
        scanf("%d %d", &numparsS, &numparsT);

        //numparsS = 10000;
        //numparsT = 10000;

        //printf("Enter theta: \n");
        //scanf("%lf", &theta);
        
        theta = 0.75;

        //printf("Enter order: \n");
        //scanf("%d", &order);
        
        order = 20;

        //printf("Enter maxparnode for source and target: \n");
        //scanf("%d %d", &maxparnodeS, &maxparnodeT);
        
        maxparnodeS = 500;
        maxparnodeT = 500;

        //printf("Enter treelevel for source and target: \n");
        //scanf("%d %d", &treelevelS, &treelevelT);
        
        treelevelS = 5;
        treelevelT = 5;

        //printf("Enter shrink for source and target (0-No 1-Yes): \n");
        //scanf("%d %d", &shrinkS, &shrinkT);
        
        shrinkS = 1;
        shrinkT = 1;

        //printf("Enter IFLAG for source and target, where options are... \n");
        //printf("(0-divide tree till number of particles in leaf is less"
        //       "or equal to maxparnode) \n");
        //printf("(1-divide tree till maxlevel attained): \n");
        //scanf("%d %d", &iflagS, &iflagT);
        
        iflagS = 0;
        iflagT = 0;


        printf("Enter kappa value: \n");
        scanf("%lf", &kappa);

        //kappa = 0.00;


        printf("Enter potential type (0 Coulomb; 1 Yukawa): \n");
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

        make_vector(tenergy,numparsT);
        make_vector(denergy,numparsT);

        printf("tenergy, denergy allocated\n\n");


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

        /* Reading in the values for exact energy at each target */
        
        fp = fopen(sampin3, "r");
        if (!fp) {
                perror("File opening failed");
                return EXIT_FAILURE;
        }

        fscanf(fp, " %lf", &time_direct);
        for (int i = 0; i < numparsT; i++)
                fscanf(fp, " %d %lf", &ii, &denergy[i]);

        fclose(fp);

        /* Calling main treecode subroutine to calculate approximate energy */

        fp = fopen(sampout, "a");

        treecode(xS, yS, zS, qS, xT, yT, zT, numparsS, numparsT,
                 tenergy, &tpeng, order, theta, shrinkS, shrinkT,
                 maxparnodeS, maxparnodeT, &time_tree,
                 treelevelS, treelevelT, iflagS, iflagT,
                 pot_type, kappa);

        /* Printing direct and treecode time calculations: */
        printf("                   Direct time (s):  %f\n", time_direct);
        printf("                     Tree time (s):  %f\n\n", time_tree);


        /* Calculating value dpeng by summing all values in denergy */ 

        dpeng = sum(denergy, numparsT);

        printf("           Direct potential energy:  %f\n", dpeng);
        printf("             Tree potential energy:  %f\n\n", tpeng);


        /* Calculating error in potential energy */

        printf("Absolute error for total potential:  %e\n", fabs(tpeng-dpeng));
        printf("Relative error for total potential:  %e\n\n", fabs((tpeng-dpeng)/dpeng));

        /* Computing pointwise potential errors */

        inferr = 0.0;
        relinferr = 0.0;
        n2err = 0.0;
        reln2err = 0.0;

        for (j = 0; j < numparsT; j++) {
                temp = fabs(denergy[j] - tenergy[j]);

                if (temp >= inferr) 
                        inferr = temp;

                if (fabs(denergy[j]) >= relinferr)
                        relinferr = fabs(denergy[j]);

                n2err = n2err + pow(denergy[j] - tenergy[j], 2.0);
                reln2err = reln2err + pow(denergy[j], 2.0);
        }

        relinferr = inferr / relinferr;
        reln2err = sqrt(n2err / reln2err);
        n2err = sqrt(n2err);

        printf("Absolute inf norm error in potential:  %e \n", inferr);
        printf("Relative inf norm error in potential:  %e \n\n", relinferr);
        printf("  Absolute 2 norm error in potential:  %e \n", n2err);
        printf("  Relative 2 norm error in potential:  %e \n\n", reln2err);

        return 0;

}
