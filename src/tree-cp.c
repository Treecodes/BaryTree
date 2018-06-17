/*
 *Procedures for Cluster-Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"


/* Definition of variables declared extern in globvars.h */
double *cf = NULL;
double *cf1 = NULL;
double *cf2 = NULL;
double ***b1 = NULL;

int *orderarr = NULL;

int torder, torderlim, torderflat;
int minlevel, maxlevel;
int maxpars, minpars;
int numleaves;

int numnodes;

double tarpos[3];
double thetasq, tarposq;

/* variables used by Yukawa treecode */
double *cf3 = NULL;
double ***a1 = NULL;


void setup(struct particles *particles, int order, double theta,
           double *xyzminmax)
{
    /* local variables */
    int i;
    double t1;

    /* changing values of our extern variables */
    torder = order;
    torderlim = torder + 1;
    thetasq = theta * theta;
    torderflat = torderlim * (torderlim + 1) * (torderlim + 2) / 6;

    /* allocating global Taylor expansion variables */
    make_vector(cf, torder+1);
    make_vector(cf1, torderlim);
    make_vector(cf2, torderlim);

    make_3array(b1, torderlim, torderlim, torderlim);


    /* initializing arrays for Taylor sums and coefficients */
    for (i = 0; i < torder + 1; i++)
        cf[i] = -i + 1.0;

    for (i = 0; i < torderlim; i++) {
        t1 = 1.0 / (i + 1.0);
        cf1[i] = 1.0 - (0.5 * t1);
        cf2[i] = 1.0 - t1;
    }

    /* find bounds of Cartesian box enclosing the particles */
    xyzminmax[0] = minval(particles->x, particles->num);
    xyzminmax[1] = maxval(particles->x, particles->num);
    xyzminmax[2] = minval(particles->y, particles->num);
    xyzminmax[3] = maxval(particles->y, particles->num);
    xyzminmax[4] = minval(particles->z, particles->num);
    xyzminmax[5] = maxval(particles->z, particles->num);

    make_vector(orderarr, particles->num);

    for (i = 0; i < particles->num; i++)
        orderarr[i] = i+1;

    return;
    
} /* END of function setup */




void cp_create_tree_n0(struct tnode **p, struct particles *targets,
                       int ibeg, int iend, int maxparnode,
                       double *xyzmm, int level)
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

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;

    (*p)->x_min = minval(targets->x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(targets->x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(targets->y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(targets->y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(targets->z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(targets->z + ibeg - 1, (*p)->numpar);

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

        cp_partition_8(targets->x, targets->y, targets->z,
                       xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                cp_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                                  targets, ind[i][0], ind[i][1],
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

} /* end of function create_tree_n0 */



void cp_partition_8(double *x, double *y, double *z, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{

    /* local variables */
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {
        cp_partition(x, y, z, orderarr, ind[0][0], ind[0][1],
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
            cp_partition(y, x, z, orderarr, ind[i][0], ind[i][1],
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
            cp_partition(z, x, y, orderarr, ind[i][0], ind[i][1],
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

} /* END of function cp_partition_8 */




void cp_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets,
                 double *tpeng, double *EnP, double *timetree)
{
    /* local variables */
    int i, j;
    double time1, time2, time3;

    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;

    time1 = MPI_Wtime();
    
    for (i = 0; i < sources->num; i++) {
        tarpos[0] = sources->x[i];
        tarpos[1] = sources->y[i];
        tarpos[2] = sources->z[i];
        tarposq = sources->q[i];

        for (j = 0; j < p->num_children; j++)
            compute_cp1(p->child[j], EnP, targets->x, targets->y, targets->z);
    }
    
    time2 = MPI_Wtime();

    compute_cp2(p, targets->x, targets->y, targets->z, EnP);
    
    time3 = MPI_Wtime();
    timetree[0] = time2 - time1;
    timetree[1] = time3 - time2;

    *tpeng = sum(EnP, targets->num);

    return;

} /* END of function cp_treecode */




void compute_cp1(struct tnode *p, double *EnP,
                 double *x, double *y, double *z)
{
    /* local variables */
    double tx, ty, tz, distsq;
    int i;

    /* determine DISTSQ for MAC test */
    tx = tarpos[0] - p->x_mid;
    ty = tarpos[1] - p->y_mid;
    tz = tarpos[2] - p->z_mid;
    distsq = tx*tx + ty*ty + tz*tz;

    if ((p->sqradius < distsq * thetasq) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */
        
        comp_tcoeff(tx, ty, tz);

        if (p->exist_ms == 0) {
            make_vector(p->ms, torderflat);

            for (i = 0; i < torderflat; i++)
                p->ms[i] = 0.0;
            
            p->exist_ms = 1;
        }

        cp_comp_ms(p);

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            cp_comp_direct(EnP, p->ibeg, p->iend, x, y, z);
        } else {
            for (i = 0; i < p->num_children; i++)
                compute_cp1(p->child[i], EnP, x, y, z);
        }
    }

    return;

} /* END of function compute_cp1 */




void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP)
{
    /* local variables */
    double tx, ty, peng;
    double xm, ym, zm, dx, dy, dz;
    int i, nn, j, k1, k2, k3, kk, porder, porder1;

    porder = torder;
    porder1 = porder - 1;

    if (ap->exist_ms == 1) {
        xm = ap->x_mid;
        ym = ap->y_mid;
        zm = ap->z_mid;

        for (i = ap->ibeg-1; i < ap->iend; i++) {
            nn = orderarr[i];
            dx = x[i] - xm;
            dy = y[i] - ym;
            dz = z[i] - zm;

            kk=0;
            peng = ap->ms[kk];

            for (k3 = porder1; k3 > -1; k3--) {
                ty = ap->ms[++kk];

                for (k2 = porder1 - k3; k2 > -1; k2--) {
                    tx = ap->ms[++kk];

                    for (k1 = porder1 - k3 - k2; k1 > -1; k1--) {
                        tx = dx*tx + ap->ms[++kk];
                    }

                    ty = dy * ty + tx;
                }

                peng = dz * peng + ty;
            }

            EnP[nn-1] = EnP[nn-1] + peng;
        }
    }

    for (j = 0; j < ap->num_children; j++)
        compute_cp2(ap->child[j], x, y, z, EnP);

    return;

} /* END function compute_cp2 */




void comp_tcoeff(double dx, double dy, double dz)
{
    /* local variables */
    double tdx, tdy, tdz, fac, sqfac;
    int i, j, k, i1, i2, j1, j2, k1, k2;

    /* setup variables */
    tdx = 2.0 * dx;
    tdy = 2.0 * dy;
    tdz = 2.0 * dz;
    fac = 1.0 / (dx*dx + dy*dy + dz*dz);
    sqfac = sqrt(fac);

    /* 0th coefficient or function val */
    b1[0][0][0] = sqfac;

    /* set of indices for which two of them are 0 */
    b1[1][0][0] = fac * dx * sqfac;
    b1[0][1][0] = fac * dy * sqfac;
    b1[0][0][1] = fac * dz * sqfac;

    for (i = 2; i < torderlim; i++) {
        i1 = i - 1;
        i2 = i - 2;
        
        b1[i][0][0] = fac * (tdx * cf1[i1] * b1[i1][0][0]
                                 - cf2[i1] * b1[i2][0][0]);

        b1[0][i][0] = fac * (tdy * cf1[i1] * b1[0][i1][0]
                                 - cf2[i1] * b1[0][i2][0]);

        b1[0][0][i] = fac * (tdz * cf1[i1] * b1[0][0][i1]
                                 - cf2[i1] * b1[0][0][i2]);
    }

    /* set of indices for which one is 0, one is 1, and other is >= 1 */
    b1[1][1][0] = fac * (dx * b1[0][1][0] + tdy * b1[1][0][0]);
    b1[1][0][1] = fac * (dx * b1[0][0][1] + tdz * b1[1][0][0]);
    b1[0][1][1] = fac * (dy * b1[0][0][1] + tdz * b1[0][1][0]);

    for (i = 2; i < torderlim - 1; i++) {

        i1 = i - 1;
        i2 = i - 2;

        b1[1][0][i] = fac * (dx * b1[0][0][i] + tdz * b1[1][0][i1]
                                - b1[1][0][i2]);

        b1[0][1][i] = fac * (dy * b1[0][0][i] + tdz * b1[0][1][i1]
                                - b1[0][1][i2]);

        b1[0][i][1] = fac * (dz * b1[0][i][0] + tdy * b1[0][i1][1]
                                - b1[0][i2][1]);

        b1[1][i][0] = fac * (dx * b1[0][i][0] + tdy * b1[1][i1][0]
                                - b1[1][i2][0]);

        b1[i][1][0] = fac * (dy * b1[i][0][0] + tdx * b1[i1][1][0]
                                - b1[i2][1][0]);

        b1[i][0][1] = fac * (dz * b1[i][0][0] + tdx * b1[i1][0][1]
                                - b1[i2][0][1]);
    }

    /* set of indices for which one is 0, others are >=2 */
    for (i = 2; i < torderlim - 2; i++) {
        i1 = i - 1;
        i2 = i - 2;

        for (j = 2; j < torderlim - i; j++) {
            j1 = j - 1;
            j2 = j - 2;

            b1[i][j][0] = fac * (tdx * cf1[i1] * b1[i1][j][0]
                                         + tdy * b1[i][j1][0]
                                     - cf2[i1] * b1[i2][j][0]
                                               - b1[i][j2][0]);

            b1[i][0][j] = fac * (tdx * cf1[i1] * b1[i1][0][j]
                                         + tdz * b1[i][0][j1]
                                     - cf2[i1] * b1[i2][0][j]
                                              - b1[i][0][j2]);

            b1[0][i][j] = fac * (tdy * cf1[i1] * b1[0][i1][j]
                                         + tdz * b1[0][i][j1]
                                     - cf2[i1] * b1[0][i2][j]
                                               - b1[0][i][j2]);
        }
    }

    /* set of indices for which two are 1, other is >= 1 */
    b1[1][1][1] = fac * (dx * b1[0][1][1] + tdy * b1[1][0][1]
                                          + tdz * b1[1][1][0]);

    for (i = 2; i < torderlim - 2; i++) {

        i1 = i - 1;
        i2 = i - 2;

        b1[1][1][i] = fac * (dx * b1[0][1][i] + tdy * b1[1][0][i]
                          + tdz * b1[1][1][i1]      - b1[1][1][i2]);

        b1[1][i][1] = fac * (dx * b1[0][i][1] + tdy * b1[1][i1][1]
                          + tdz * b1[1][i][0]       - b1[1][i2][1]);

        b1[i][1][1] = fac * (dy * b1[i][0][1] + tdx * b1[i1][1][1]
                        + tdz * b1[i][1][0]       - b1[i2][1][1]);
    }

    /* set of indices for which one is 1, others are >= 2 */
    for (i = 2; i < torderlim - 3; i++) {
        i1 = i - 1;
        i2 = i - 2;

        for (j = 2; j < torderlim - i - 1; j++) {
            j1 = j - 1;
            j2 = j - 2;

            b1[1][i][j] = fac * (dx * b1[0][i][j] + tdy * b1[1][i1][j]
                              + tdz * b1[1][i][j1]      - b1[1][i2][j]
                                    - b1[1][i][j2]);

            b1[i][1][j] = fac * (dy * b1[i][0][j] + tdx * b1[i1][1][j]
                              + tdz * b1[i][1][j1]      - b1[i2][1][j]
                                    - b1[i][1][j2]);

            b1[i][j][1] = fac * (dz * b1[i][j][0] + tdx * b1[i1][j][1]
                              + tdy * b1[i][j1][1]      - b1[i2][j][1]
                                    - b1[i][j2][1]);
        }
    }

    /* set of indices for which all are >= 2 */
  
    for (k = 2; k < torderlim - 4; k++) {
        k1 = k - 1;
        k2 = k - 2;
                
        for (j = 2; j < torderlim - k - 2; j++) {
            j1 = j - 1;
            j2 = j - 2;

            for (i = 2; i < torderlim - k - j - 1; i++) {
                i1 = i - 1;
                i2 = i - 2;

                b1[i][j][k] = fac * (tdx * cf1[i1] * b1[i1][j][k]
                                             + tdy * b1[i][j1][k]
                                             + tdz * b1[i][j][k1]
                                         - cf2[i1] * b1[i2][j][k]
                                    - b1[i][j2][k] - b1[i][j][k2]);
            }
        }
    }

    return;

} /* END function comp_tcoeff */




/*
 * comp_cms computes the moments for node p needed in the Taylor approximation
 */
void cp_comp_ms(struct tnode *p)
{
    int k1, k2, k3, kk=-1;

    for (k3 = torder; k3 > -1; k3--) {
        for (k2 = torder - k3; k2 > -1; k2--) {
            for (k1 = torder - k3 - k2; k1 > -1; k1--) {
                kk++;
                p->ms[kk] = p->ms[kk]
                          + tarposq * b1[k1][k2][k3];
            }
        }
    }

    return;

} /* END function comp_cms */



/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void cp_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z)
{

    /* local variables */
    int i, nn;
    double tx, ty, tz;

    for (i = ibeg - 1; i < iend; i++) {
        nn = orderarr[i];
        tx = x[i] - tarpos[0];
        ty = y[i] - tarpos[1];
        tz = z[i] - tarpos[2];
        EnP[nn-1] = EnP[nn-1] + tarposq / sqrt(tx*tx + ty*ty + tz*tz);
    }

    return;

} /* END function cp_comp_direct */




/*
 * cleanup deallocates allocated global variables and then calls
 * recursive function remove_node to delete the tree
 */
void cleanup(struct tnode *p)
{
    free_vector(cf);
    free_vector(cf1);
    free_vector(cf2);
    free_vector(cf3);
    
    free_3array(b1);
    free_3array(a1);

    free_vector(orderarr);

    remove_node(p);
    free(p);

    return;

} /* END function cleanup */




void remove_node(struct tnode *p)
{
    /* local variables */
    int i;

    if (p->exist_ms == 1)
        free(p->ms);

    if (p->num_children > 0) {
        for (i = 0; i < p->num_children; i++) {
            remove_node(p->child[i]);
            free(p->child[i]);
        }
    }

    return;

} /* END function remove_node */
