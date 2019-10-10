#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include “zoltan.h”



void load_balance(double *x, double *y, double *z, double *q, double *w, int numpars_local)
{
//    /* local variables */
//	int rank, numProcs;
//
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//
//
//
////  Gather an array containing the number of points on each processor
//	int * pointsOnEachProc;
//	make_vector(pointsOnEachProc, numProcs);
//	MPI_Allgather(&numpars_local, 1, MPI_INT, pointsOnEachProc, 1, MPI_INT, MPI_COMM_WORLD);
//
////	Each processor computes its offset
//	int localOffset=0;
//	for (int i=0;i<rank;i++){
//		localOffset += pointsOnEachProc[i];
//	}


	/* Initialize Zoltan */
	rc = Zoltan_Initialize(argc, argv, &ver);
	if (rc != ZOLTAN_OK){
	printf("sorry...\n");
	free(Pts); free(Gids);
	exit(0);
	}

	/* Allocate and initialize memory for Zoltan structure */
	zz = Zoltan_Create(MPI_COMM_WORLD);
	/* Set general parameters */
	Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
	Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
	Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1");
	Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
	Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
	/* Set RCB parameters */
	Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
	Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
	Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1");
	/* Register call-back query functions. */
	Zoltan_Set_Num_Obj_Fn(zz, exGetNumberOfAssignedObjects, NULL);
	Zoltan_Set_Obj_List_Fn(zz, exGetObjectList, NULL);
	Zoltan_Set_Num_Geom_Fn(zz, exGetObjectSize, NULL);
	Zoltan_Set_Geom_Multi_Fn(zz, exGetObject, NULL);


	/* Perform partitioning */
	rc = Zoltan_LB_Partition(zz,
	&changes, /* Flag indicating whether partition changed */
	&numGidEntries, &numLidEntries,
	&numImport, /* objects to be imported to new part */
	&importGlobalGids, &importLocalGids,
	&importProcs, &importToPart,
	&numExport, /* # objects to be exported from old part */
	&exportGlobalGids, &exportLocalGids,
	&exportProcs, &exportToPart);

	/* Process partitioning results;
	** in this case, print information;
	** in a "real" application, migrate data here.
	*/
	if (!rc){
	exPrintGlobalResult("Recursive Coordinate Bisection",
	nprocs, me,
	MyNumPts, numImport, numExport, changes);
	}
	else{
	free(Pts);
	free(Gids);
	Zoltan_Destroy(&zz);
	MPI_Finalize();
	exit(0);
	}

	/* Free Zoltan memory allocated by Zoltan_LB_Partition. */
	Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids,
	&importProcs, &importToPart);
	Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids,
	&exportProcs, &exportToPart);
	/* Free Zoltan memory allocated by Zoltan_Create. */
	Zoltan_Destroy(&zz);
	/* Free Application memory */
	free(Pts); free(Gids);
	/**********************
	** all done ***********
	**********************/



//	char *lb_method;
//	int new, num_imp, num_exp, *imp_procs, *exp_procs;
//	int *imp_to_part, *exp_to_part;
//	int num_gid_entries, num_lid_entries;
//	ZOLTAN_ID_PTR imp_global_ids, exp_global_ids;
//	ZOLTAN_ID_PTR imp_local_ids, exp_local_ids;
//	/* Set load-balancing method. */
//	read_load_balancing_info_from_input_file(&lb_method);
//	Zoltan_Set_Param(zz, "LB_METHOD", lb_method);
//
//	/* Reset some load-balancing parameters. */
//	Zoltan_Set_Param(zz, "RCB_Reuse", "TRUE");
//
//	/* Perform computations */
//	...
//	/* Perform load balancing */
//	Zoltan_LB_Partition(zz,&new,&num_gid_entries,&num_lid_entries,
//	    &num_imp,&imp_global_ids,&imp_local_ids,&imp_procs,&imp_to_part,
//	    &num_exp,&exp_global_ids,&exp_local_ids,&exp_procs,&exp_to_part);
//	if (new)
//	  perform_data_migration(...);
//
//	/* Free memory allocated for load-balancing results by Zoltan library */
//	Zoltan_LB_Free_Part(&imp_global_ids, &imp_local_ids, &imp_procs, &imp_to_part);
//	Zoltan_LB_Free_Part(&exp_global_ids, &exp_local_ids, &exp_procs, &exp_to_part);


    return;

} /* END of function load_balance */


