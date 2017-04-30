#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "array.h"
#include "sort.h"

int main()
{

    /* input and output files */
    char sampin1[20], sampout[20];
    FILE *fp;

    //local variables (these decs may be unnecessary in C)
    int i, numparsS, intype;
    char c1[10], c2[10], c3[10], c4[10], c5[10], c6[10], c7[10];
    double a1, a2, a3, a4, a5;
	double *xS, *yS, *zS, *qS, *mS;

    printf("Enter name of input file: \n");
    scanf("%s", sampin1);
    //sampin1 = "S10000.txt";

    printf("Enter name of output file: \n");
    scanf("%s", sampout);
    //sampout = "S10000.bin";
    
    printf("Enter type of input file (0--sources, 1--targets): \n");
    scanf("%d", &intype);
    //intype = 1;

    printf("Enter particle number: \n");
    scanf("%d", &numparsS);

    make_vector(xS,numparsS);
    make_vector(yS,numparsS);
    make_vector(zS,numparsS);
    make_vector(qS,numparsS);

    printf("xS, yS, zS, qS, mS allocated\n");

    printf("Reading in sources...\n");

    i = 0;
    fp = fopen(sampin1, "r");
    if (intype == 0)
    {
        while (fscanf(fp, " %lf %lf %lf %lf",
                  &a1, &a2, &a3, &a4) != EOF) {
            xS[i] = a1;
            yS[i] = a2;
            zS[i] = a3;
            qS[i] = a4;
            i++;
        }
    } else {
        while (fscanf(fp, " %lf %lf %lf",
                      &a1, &a2, &a3) != EOF) {
            xS[i] = a1;
            yS[i] = a2;
            zS[i] = a3;
            i++;
        }
    }
	fclose(fp);

    sortTargets(xS, yS, zS, numparsS);

    printf("Printing sources...\n");

    fp = fopen(sampout, "wb");
    for (i = 0; i < numparsS; i++) {
        //fprintf(fp, "%e  %e  %e  %e  %e  \n", 
        //        xS[i], yS[i], zS[i], qS[i], mS[i]);
        fwrite(&xS[i], sizeof(double), 1, fp);
        fwrite(&yS[i], sizeof(double), 1, fp);
        fwrite(&zS[i], sizeof(double), 1, fp);
        if (intype == 0)
            fwrite(&qS[i], sizeof(double), 1, fp);
    }
    fclose(fp);

    free_vector(xS);
    free_vector(yS);
    free_vector(zS);
    free_vector(qS);

	return 0;
}
