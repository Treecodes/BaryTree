#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <omp.h>
#include <string.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "particles.h"
#include "sort.h"


/* The treedriver routine in Fortran */
int main(int argc, char **argv)
{
	printf("Welcome to BaryTree.\n");

    /* runtime parameters */
    int numparsS, numparsT;
    int order, maxparnode, batch_size;
    int pot_type, tree_type;
    int numDevices, numThreads;
    double theta, temp, kappa;

    /* particles */
    struct particles *sources = NULL, *targets = NULL;

    /* exact energy, treecode energy */
    double *denergy = NULL, *tenergy = NULL;

    /* for potential energy calculation */
    double tpeng = 0, dpeng = 0;

    /* insert variables for date-time calculation? */
    double time_direct, time_tree[4], time_preproc, time_treedriver;
    double time1, time2;

    /* input and output files */
    char *sampin1 = NULL;
    char *sampin2 = NULL;
    char *sampin3 = NULL;
    char *sampout = NULL;
    FILE *fp;

    /* local variables */
    double buf[5];
    

    /* Executable statements begin here */
    sampin1 = argv[1];
    if (strcmp(sampin1,"--help") == 0)
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
//        printf("     tree_type:  0--cluster-particle, 1--particle-cluster \n");  // 0
        printf("    maxparnode:  maximum particles in leaf \n");                 // 500
        printf("         kappa:  screened Coulomb parameter \n");                // 0.00
        printf("      pot_type:  0--Coulomb, 1--screened Coulomb \n");           // 1
        printf("    batch_size:  size of target batches \n");     // 0
        printf("   num devices:  number of GPUs available \n");     // 0
        printf("   num threads:  number of OpenMP threads \n");     // 0

        return 0;
    }
    
    sampin2 = argv[2];
    sampin3 = argv[3];
    sampout = argv[4];
    numparsS = atoi(argv[5]);
    numparsT = atoi(argv[6]);
    theta = atof(argv[7]);
    order = atoi(argv[8]);
    tree_type = 1;
    maxparnode = atoi(argv[9]);
    kappa = atof(argv[10]);
    pot_type = atoi(argv[11]);
    batch_size = atoi(argv[12]);
    numDevices = atoi(argv[13]);
    numThreads = atoi(argv[14]);
    
    time1 = omp_get_wtime();
    
    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));
    
    sources->num = numparsS;
    make_vector(sources->x, numparsS);
    make_vector(sources->y, numparsS);
    make_vector(sources->z, numparsS);
    make_vector(sources->q, numparsS);
    make_vector(sources->w, numparsS);

    targets->num = numparsT;
    make_vector(targets->x, numparsT);
    make_vector(targets->y, numparsT);
    make_vector(targets->z, numparsT);
    make_vector(targets->q, numparsT);
    make_vector(targets->order, numparsT);
    
    make_vector(tenergy, numparsT);
    make_vector(denergy, numparsT);


    /* Reading in coordinates and charges for the source particles*/
    fp = fopen(sampin1, "rb");
    for (int i = 0; i < numparsS; ++i) {
        fread(buf, sizeof(double), 5, fp);
        sources->x[i] = buf[0];
        sources->y[i] = buf[1];
        sources->z[i] = buf[2];
        sources->q[i] = buf[3];
        sources->w[i] = buf[4];
    }
    fclose(fp);

    fp = fopen(sampin2, "rb");
    for (int i = 0; i < numparsT; ++i) {
        fread(buf, sizeof(double), 4, fp);
        targets->x[i] = buf[0];
        targets->y[i] = buf[1];
        targets->z[i] = buf[2];
        targets->q[i] = buf[3];
        targets->order[i] = i;
    }
    fclose(fp);

    fp = fopen(sampin3, "rb");
    fread(&time_direct, sizeof(double), 1, fp);
    fread(denergy, sizeof(double), numparsT, fp);
    fclose(fp);


    // Initialize all GPUs
	if (numDevices > 0) {
		#pragma omp parallel num_threads(numDevices)
        {
			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
			acc_init(acc_get_device_type());
        }
	}
    
    time2 = omp_get_wtime();
    time_preproc = time2 - time1;
    
    printf("Setup complete, calling treedriver...\n");
    
    /* Calling main treecode subroutine to calculate approximate energy */
    time1 = omp_get_wtime();
    treedriver(sources, targets,
               order, theta, maxparnode, batch_size,
               pot_type, kappa, tree_type,
               tenergy, &tpeng, time_tree, numDevices, numThreads);
    time2 = omp_get_wtime();
    time_treedriver = time2 - time1;


    dpeng = sum(denergy, numparsT);

    /* Printing direct and treecode time calculations: */
    printf("                   Direct time (s):  %f\n\n", time_direct);
    printf("              Pre-process time (s):  %f\n", time_preproc);
    printf("               Treedriver time (s):  %f\n", time_treedriver);
    printf("               Tree setup time (s):  %f, %f\n", time_tree[0]);
    
    if (tree_type == 0) {
        printf("             Min, Max cp1 time (s):  %f, %f\n", time_tree[1]);
        printf("             Min, Max cp2 time (s):  %f, %f\n", time_tree[2]);
    }
        
    printf("      Min, Max total tree time (s):  %f, %f\n\n", time_tree[3]);
    printf(" Preproc + Max total tree time (s):  %f \n\n", time_tree[3] + time_preproc);
    
    printf("           Direct potential energy:  %f\n", dpeng);
    printf("             Tree potential energy:  %f\n\n", tpeng);
    
    printf("Absolute error for total potential:  %e\n", fabs(tpeng-dpeng));
    printf("Relative error for total potential:  %e\n\n", fabs((tpeng-dpeng)/dpeng));
    
    
    /* Computing pointwise potential errors */

    double inferr = 0.0, relinferr = 0.0;
    double n2err = 0.0, reln2err = 0.0;
    double xinf, yinf, zinf;

    for (int j = 0; j < numparsT; ++j) {
        temp = fabs(denergy[targets->order[j]] - tenergy[j]);
        
        if (temp >= inferr) {
            inferr = temp;
            xinf = targets->x[targets->order[j]];
            yinf = targets->y[targets->order[j]];
            zinf = targets->z[targets->order[j]];
        }

        if (fabs(denergy[j]) >= relinferr)
            relinferr = fabs(denergy[j]);

        n2err += pow(denergy[targets->order[j]]
                    - tenergy[j], 2.0) * sources->w[j];
        reln2err += pow(denergy[j], 2.0) * sources->w[j];
    }

    relinferr = inferr / relinferr;
    reln2err = sqrt(n2err / reln2err);
    n2err = sqrt(n2err);

    printf("Absolute inf norm error in potential:  %e \n", inferr);
    printf("Relative inf norm error in potential:  %e \n\n", relinferr);
    printf("  Absolute 2 norm error in potential:  %e \n", n2err);
    printf("  Relative 2 norm error in potential:  %e \n\n", reln2err);
    printf("              Inf error occurring at: %f, %f, %f \n\n", xinf, yinf, zinf);
    

    fp = fopen(sampout, "a");
    fprintf(fp, "%s,%s,%s,%d,%d,%f,%d,%d,%d,%d,%f,%d,"
        "%f,%f,%f,%f,%f,%f,"
        "%e,%e,%e,%e,%e,%e,%e,%e,%d,%d\n",
        sampin1, sampin2, sampin3, numparsS, numparsT,
        theta, order, tree_type, maxparnode,batch_size,
        kappa, pot_type, //1 ends
        time_preproc, time_tree[0], time_tree[1], time_tree[2], time_tree[3],
        time_tree[3] + time_preproc, //2 ends
        dpeng, tpeng, fabs(tpeng-dpeng), fabs((tpeng-dpeng)/dpeng),
        inferr, relinferr, n2err, reln2err, numDevices, numThreads); //3 ends
    fclose(fp);
    
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

    free_vector(denergy);
    free_vector(tenergy);

    printf("BaryTree has finished.\n");
    
    return 0;
}
