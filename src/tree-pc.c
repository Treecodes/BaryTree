/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"


void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level)
{

    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild;
    
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


    /* increment number of nodes */
    (*p)->node_index = numnodes;
    numnodes++;

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;
    
    (*p)->x_min = minval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(sources->z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(sources->z + ibeg - 1, (*p)->numpar);
    
    //(*p)->x_min = xyzmm[0];
    //(*p)->x_max = xyzmm[1];
    //(*p)->y_min = xyzmm[2];
    //(*p)->y_max = xyzmm[3];
    //(*p)->z_min = xyzmm[4];
    //(*p)->z_max = xyzmm[5];
    

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

        pc_partition_8(sources->x, sources->y, sources->z, sources->q,
                       xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                pc_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                                  sources, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev);
            }
        }
        
    } else {
        
        if (level < minlevel) minlevel = level;
        if (minpars > (*p)->numpar) minpars = (*p)->numpar;
        if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
        /* increment number of leaves */
        numleaves++;
    }

    return;

} /* END of function create_tree_n0 */



void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array)
{
    int i;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    tree_array->x_mid[p->node_index] = p->x_mid;
    tree_array->y_mid[p->node_index] = p->y_mid;
    tree_array->z_mid[p->node_index] = p->z_mid;

    tree_array->ibeg[p->node_index] = p->ibeg;
    tree_array->iend[p->node_index] = p->iend;

    for (i = 0; i < p->num_children; i++) {
        pc_create_tree_array(p->child[i], tree_array);
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




void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets,
                 double *tpeng, double *EnP)
{
    /* local variables */
    int i, j;

    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;
    
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q,
                targets->x, targets->y, targets->z, EnP);
        }
    }

    *tpeng = sum(EnP, targets->num);

    return;

} /* END of function pc_treecode */




void compute_pc(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS,
                double *xT, double *yT, double *zT, double *EnP)
{
    /* local variables */
    double dxt, dyt, dzt, dist;
    double tx, ty, tz;
    int i, j, k, kk, ii;
    double *temp_i, *temp_j, *temp_k;

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);

    
    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than n0 particles
     * in the box, use the expansion for the approximation.
     */
     
        make_vector(temp_i, torderlim);
        make_vector(temp_j, torderlim);
        make_vector(temp_k, torderlim);

        if (p->exist_ms == 0) {
            make_vector(p->ms, (torderlim)*(torderlim)*(torderlim));
            make_vector(p->tx, torderlim);
            make_vector(p->ty, torderlim);
            make_vector(p->tz, torderlim);
            

            for (i = 0; i < (torderlim)*(torderlim)*(torderlim); i++)
                p->ms[i] = 0.0;

            pc_comp_ms(p, xS, yS, zS, qS);
            p->exist_ms = 1;
        }
        
        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
        
            for (i = 0; i < torderlim; i++) {
                dxt = xT[ii] - p->tx[i];
                dyt = yT[ii] - p->ty[i];
                dzt = zT[ii] - p->tz[i];
                temp_i[i] = dxt * dxt;
                temp_j[i] = dyt * dyt;
                temp_k[i] = dzt * dzt;
            }
        
            kk = -1;
            for (k = 0; k < torderlim; k++) {
                for (j = 0; j < torderlim; j++) {
                    for (i = 0; i < torderlim; i++) {
                        kk++;
                        EnP[ii] += p->ms[kk] / sqrt(temp_i[i] + temp_j[j] + temp_k[k]);
                    }
                }
            }
        }
        
        free_vector(temp_i);
        free_vector(temp_j);
        free_vector(temp_k);
        
    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            pc_comp_direct(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, xT, yT, zT, EnP);
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_pc(p->child[i], batch_ind, batch_mid, batch_rad,
                           xS, yS, zS, qS, xT, yT, zT, EnP);
            }
        }
    }

    return;

} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *restrict xS, double *restrict yS, double *restrict zS, double *restrict qS,
                    double *restrict xT, double *restrict yT, double *restrict zT, double *restrict EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;
	
    double d_peng;

    #pragma acc data present(xS, yS, zS, qS)
    #pragma acc kernels loop
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
        d_peng = 0.0;
        for (i = ibeg - 1; i < iend; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            
            d_peng += qS[i] / sqrt(tx*tx + ty*ty + tz*tz);
        }
        EnP[ii] += d_peng;
    }

    return;

} /* END function pc_comp_direct */




/*
 * cp_comp_ms computes the moments for node p needed in the Taylor approximation
 */
void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q)
{

    int i, j, k1, k2, k3, kk;
    int a1exactind, a2exactind, a3exactind;
    double dx, dy, dz, tx, ty, tz, qloc;
    double x0, x1, y0, y1, z0, z1;
    double sumA1, sumA2, sumA3;
    double xx, yy, zz;
    double *xibeg, *yibeg, *zibeg, *qibeg;
    
    double *w1i, *w2j, *w3k, *dj, *Dd;
    double **a1i, **a2j, **a3k;
    
    xibeg = &(x[p->ibeg-1]);
    yibeg = &(y[p->ibeg-1]);
    zibeg = &(z[p->ibeg-1]);
    qibeg = &(q[p->ibeg-1]);
    
    x0 = p->x_min;
    x1 = p->x_max;
    y0 = p->y_min;
    y1 = p->y_max;
    z0 = p->z_min;
    z1 = p->z_max;
    
    for (i = 0; i < torderlim; i++) {
        p->tx[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        p->ty[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        p->tz[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    
    make_vector(w1i, torderlim);
    make_vector(w2j, torderlim);
    make_vector(w3k, torderlim);
    make_vector(dj, torderlim);
    
    make_matrix(a1i, torderlim, p->numpar);
    make_matrix(a2j, torderlim, p->numpar);
    make_matrix(a3k, torderlim, p->numpar);
    
    make_vector(Dd, p->numpar);
    
    dj[0] = 0.5;
    dj[torder] = 0.5;
    for (j = 1; j < torder; j++)
        dj[j] = 1.0;
    
    for (j = 0; j < torderlim; j++) {
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
        w2j[j] = w1i[j];
        w3k[j] = w1i[j];
    }
    
    sumA1 = 0.0;
    sumA2 = 0.0;
    sumA3 = 0.0;
    
    for (i = 0; i < p->numpar; i++) {
        xx = xibeg[i];
        yy = yibeg[i];
        zz = zibeg[i];
        
        a1exactind = -1;
        a2exactind = -1;
        a3exactind = -1;
        
        for (j = 0; j < torderlim; j++) {
            a1i[j][i] = w1i[j] / (xx - p->tx[j]);
            a2j[j][i] = w2j[j] / (yy - p->ty[j]);
            a3k[j][i] = w3k[j] / (zz - p->tz[j]);
            
            //if (isinf(a2j[j][i - p->ibeg + 1]))
            //    printf("a2j %d, %d: %f, %f, %f\n", i, j, a2j[j][i - p->ibeg + 1], yy, p->ty[j]);
            
            //if (isinf(a1i[j][i - p->ibeg + 1]))
            //    printf("a1i %d, %d: %f, %f, %f\n", i, j, a1i[j][i - p->ibeg + 1], xx, p->tx[j]);

            //if (isinf(a3k[j][i - p->ibeg + 1]))
            //    printf("a3k %d, %d: %f, %f, %f\n", i, j, a3k[j][i - p->ibeg + 1], zz, p->tz[j]);

            sumA1 += a1i[j][i];
            sumA2 += a2j[j][i];
            sumA3 += a3k[j][i];
            
            if (fabs(xx - p->tx[j]) < 1e-10) a1exactind = j;
            if (fabs(yy - p->ty[j]) < 1e-10) a2exactind = j;
            if (fabs(zz - p->tz[j]) < 1e-10) a3exactind = j;
        }
        
        if (a1exactind > -1) {
            sumA1 = 1.0;
            for (j = 0; j < torderlim; j++) a1i[j][i] = 0.0;
            a1i[a1exactind][i] = 1.0;
        }
        
        if (a2exactind > -1) {
            sumA2 = 1.0;
            for (j = 0; j < torderlim; j++) a2j[j][i] = 0.0;
            a2j[a2exactind][i] = 1.0;
        }
        
        if (a3exactind > -1) {
            sumA3 = 1.0;
            for (j = 0; j < torderlim; j++) a3k[j][i] = 0.0;
            a3k[a3exactind][i] = 1.0;
        }
        
        Dd[i] = 1.0 / (sumA1 * sumA2 * sumA3);
        
        sumA1 = 0.0;
        sumA2 = 0.0;
        sumA3 = 0.0;
    }
    
    kk = -1;
    for (k3 = 0; k3 < torderlim; k3++) {
        for (k2 = 0; k2 < torderlim; k2++) {
            for (k1 = 0; k1 < torderlim; k1++) {
                kk++;
                for (i = 0; i < p->numpar; i++) {
                    p->ms[kk] += a1i[k1][i] * a2j[k2][i] * a3k[k3][i]
                               * Dd[i] * qibeg[i];
                    
                    //if (p->ms[kk] != p->ms[kk])
                    //    printf("%d, %d, %d, %d: %f, %f, %f, %f, %f, %f\n", k1, k2, k3, i,
                    //    p->ms[kk], Dd[i], q[p->ibeg-1+i], a1i[k1][i], a2j[k2][i], a3k[k3][i]);
                }
            }
        }
    }
    
    
    
    free_vector(w1i);
    free_vector(w2j);
    free_vector(w3k);
    
    free_vector(dj);
    
    free_matrix(a1i);
    free_matrix(a2j);
    free_matrix(a3k);

    free_vector(Dd);
    
    return;
    
} /* END function cp_comp_ms */




void pc_make_interaction_list(struct tnode *p, struct batch *batches,
                              int **tree_inter_list, int **direct_inter_list)
{
    /* local variables */
    int i, j;
    int tree_index_counter;
    int direct_index_counter;

    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < numnodes; j++) {
            tree_inter_list[i][j] = -1;
        }
        for (j = 0; j < numleaves; j++) {
            direct_inter_list[i][j] = -1;
        }
    }
    
    for (i = 0; i < batches->num; i++) {
        tree_index_counter = 0;
        direct_index_counter = 0;
        
        pc_compute_interaction_list(p,
                batches->index[i], batches->center[i], batches->radius[i],
                tree_inter_list[i], direct_inter_list[i],
                &tree_index_counter, &direct_index_counter);
    }

    return;

} /* END of function pc_treecode */




void pc_compute_interaction_list(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list,
                int *tree_index_counter, int *direct_index_counter)
{
    /* local variables */
    double tx, ty, tz, dist;
    int i;

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */

        batch_tree_list[*tree_index_counter] = p->node_index;
        (*tree_index_counter)++;

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            batch_direct_list[*direct_index_counter] = p->node_index;
            (*direct_index_counter)++;

        } else {
            for (i = 0; i < p->num_children; i++) {
                pc_compute_interaction_list(p->child[i], batch_ind, batch_mid, batch_rad,
                           batch_tree_list, batch_direct_list,
                           tree_index_counter, direct_index_counter);
            }
        }
    }

    return;

} /* END of function compute_pc */
