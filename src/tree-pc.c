/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"


void pc_create_tree_n0(struct tnode **p, int ibeg, int iend,
                       double *x, double *y, double *z, double *q,
                       int maxparnode, double *xyzmm,
                       int level)
{
    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild, nump;
    
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];


    for (i = 0; i < 8; i++) {
        for (j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }
                        

    (*p) = malloc(sizeof(struct tnode));


    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;

    nump = iend - ibeg + 1;

    (*p)->x_min = minval(x + ibeg - 1, nump);
    (*p)->x_max = maxval(x + ibeg - 1, nump);
    (*p)->y_min = minval(y + ibeg - 1, nump);
    (*p)->y_max = maxval(y + ibeg - 1, nump);
    (*p)->z_min = minval(z + ibeg - 1, nump);
    (*p)->z_max = maxval(z + ibeg - 1, nump);


    /*compute aspect ratio*/
    xl = (*p)->x_max - (*p)->x_min;
    yl = (*p)->y_max - (*p)->y_min;
    zl = (*p)->z_max - (*p)->z_min;
        
    lmax = max3(xl, yl, zl);
    t1 = lmax;
    t2 = min3(xl, yl, zl);


    if (t2 != 0.0)
        (*p)->aspect = t1/t2;
    else
        (*p)->aspect = 0.0;


    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    t1 = (*p)->x_max - (*p)->x_mid;
    t2 = (*p)->y_max - (*p)->y_mid;
    t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    /*set particle limits, tree level of node, and nullify child pointers*/
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;


    if (maxlevel < level) maxlevel = level;

    (*p)->num_children = 0;
    for (i = 0; i < 8; i++)
        (*p)->child[i] = NULL;

    
    if ((*p)->numpar > maxparnode) {

    /*
     * set IND array to 0, and then call PARTITION_8 routine.
     * IND array holds indices of the eight new subregions.
     * Also, setup XYZMMS array in the case that SHRINK = 1.
     */
        xyzmms[0][0] = (*p)->x_min;
        xyzmms[1][0] = (*p)->x_max;
        xyzmms[2][0] = (*p)->y_min;
        xyzmms[3][0] = (*p)->y_max;
        xyzmms[4][0] = (*p)->z_min;
        xyzmms[5][0] = (*p)->z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        x_mid = (*p)->x_mid;
        y_mid = (*p)->y_mid;
        z_mid = (*p)->z_mid;

        pc_partition_8(x, y, z, q, xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                pc_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                                  ind[i][0], ind[i][1], x, y, z, q,
                                  maxparnode, lxyzmm, loclev);
            }
        }

    } else {

        if (level < minlevel) minlevel = level;
        if (minpars > (*p)->numpar) minpars = (*p)->numpar;
        if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
        numleaves++;
    }

    return;

} /* END of function create_tree_n0 */



void pc_partition_8(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    /* local variables */
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {

        pc_partition(x, y, z, q, orderarr, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];

        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;

    }

    if (yl >= critlen) {

        for (i = 0; i < *numposchild; i++) {
            pc_partition(y, x, z, q, orderarr, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);
                        
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    if (zl >= critlen) {

        for (i = 0; i < *numposchild; i++) {
            pc_partition(z, x, y, q, orderarr, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);
                        
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function partition_8 */




void pc_treecode(struct tnode *p, double *xS, double *yS, double *zS,
                 double *qS, double *xT, double *yT, double *zT,
                 double *tpeng, double *EnP, int numparsS, int numparsT,
                 int **batch_index, double **batch_center, double *batch_radius,
                 int batch_num, int *batch_reorder)
{
    /* local variables */
    int i, j;
    double penglocal, peng;


    for (i = 0; i < numparsT; i++)
        EnP[i] = 0.0;

//    for (i = 0; i < numparsT; i++) {
//        peng = 0.0;
//        tarpos[0] = xT[i];
//        tarpos[1] = yT[i];
//        tarpos[2] = zT[i];
//
//        /* Copy target coordinates to GPU */
//        #pragma acc data copyin(tarpos[3])
//
//        for (j = 0; j < p->num_children; j++) {
//            compute_pc(p->child[j], &penglocal, xS, yS, zS, qS, numparsS);
//            peng += penglocal;
//        }
//
//        EnP[i] = peng;
//    }
    
    for (i = 0; i < batch_num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc(p->child[j], EnP, xS, yS, zS, qS, xT, yT, zT,
            batch_index[i], batch_center[i], batch_radius[i], batch_reorder);
        }
    }

    *tpeng = sum(EnP, numparsT);

    return;

} /* END of function pc_treecode */




void compute_pc(struct tnode *p, double *EnP,
                double *x, double *y, double *z, double *q,
                double *xT, double *yT, double *zT,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_reorder)
{
    /* local variables */
    double tx, ty, tz, distsq, dist, penglocal;
    int i, j, k, kk, ii, nn;

    //printf("Inside compute_cp1... 1\n");

    /* determine DISTSQ for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    distsq = tx*tx + ty*ty + tz*tz;
    dist = sqrt(distsq);

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */

        if (p->exist_ms == 0) {
            make_vector(p->ms, torderflat);

            for (i = 0; i < torderflat; i++)
                p->ms[i] = 0.0;

            pc_comp_ms(p, x, y, z, q);
            p->exist_ms = 1;
        }
        
        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
//            nn = batch_reorder[ii] - 1;
            tx = xT[ii] - p->x_mid;
            ty = yT[ii] - p->y_mid;
            tz = zT[ii] - p->z_mid;

            comp_tcoeff(tx, ty, tz);
        
            kk = -1;
            
            for (k = 0; k < torder + 1; k++) {
                for (j = 0; j < torder - k + 1; j++) {
                    for (i = 0; i < torder - k - j + 1; i++) {
                        EnP[ii] += b1[i][j][k] * p->ms[++kk];
                    }
                }
            }
        }

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            pc_comp_direct(EnP, p->ibeg, p->iend, x, y, z, q,
                           batch_ind[0], batch_ind[1], xT, yT, zT, batch_reorder);
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_pc(p->child[i], EnP, x, y, z, q, xT, yT, zT,
                           batch_ind, batch_mid, batch_rad, batch_reorder);
            }
        }
    }

    return;

} /* END of function compute_pc */



/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct(double *EnP, int ibeg, int iend,
                    double *restrict x, double *restrict y, double *restrict z, double *restrict q,
                    int batch_ibeg, int batch_iend, double *restrict xT, double *restrict yT, double *restrict zT,
                    int *restrict batch_reorder)
{
    /* local variables */
    int i, ii, nn;
    double tx, ty, tz;

/*
//    *peng = 0.0;
    double d_peng=0.0;
//    #pragma acc kernels loop worker(16) reduction(+:d_peng)
	#pragma acc data present(x[numparsS],y[numparsS],z[numparsS],q[numparsS],tarpos[3])
	#pragma acc kernels
    { 															// begin acc kernels region
//	#pragma acc kernels loop gang(32), vector(128)
    for (i = ibeg - 1; i < iend; i++) {
        tx = x[i] - tarpos[0];
        ty = y[i] - tarpos[1];
        tz = z[i] - tarpos[2];
        
        d_peng += q[i] / sqrt(tx*tx + ty*ty + tz*tz);
    }
    } 															// end acc kernels region
    *peng = d_peng;
*/

    double d_peng=0.0;
//	#pragma acc data present(x,y,z,q,xT,yT,zT)
	#pragma acc data present(x,y,z,q)
	#pragma acc region
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
//        nn = batch_reorder[ii] - 1;
        d_peng=0.0;
        for (i = ibeg - 1; i < iend; i++) {
            tx = x[i] - xT[ii];
            ty = y[i] - yT[ii];
            tz = z[i] - zT[ii];
            
            d_peng += q[i] / sqrt(tx*tx + ty*ty + tz*tz);
        }
        EnP[ii] += d_peng;
        
        //printf("%d; EnP %lf \n", ii, EnP[ii]);
    }

    return;

} /* END function pc_comp_direct */




/*
 * cp_comp_ms computes the moments for node p needed in the Taylor approximation
 */
void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q)
{

    int i, k1, k2, k3, kk;
    double dx, dy, dz, tx, ty, tz, qloc;
    
    for (i = p->ibeg-1; i < p->iend; i++) {
        dx = x[i] - p->x_mid;
        dy = y[i] - p->y_mid;
        dz = z[i] - p->z_mid;
        qloc = q[i];
        
        kk = -1;
        tz = 1.0;
        for (k3 = 0; k3 < torder + 1; k3++) {
            ty = 1.0;
            for (k2 = 0; k2 < torder - k3 + 1; k2++) {
                tx = 1.0;
                for (k1 = 0; k1 < torder - k3 - k2 + 1; k1++) {
                    p->ms[++kk] += qloc * tx*ty*tz;
                    tx *= dx;
                }
                ty *= dy;
            }
            tz *= dz;
        }
    }
    
    return;
    
} /* END function cp_comp_ms */
