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
#include "tools.h"

#include "tree.h"


/* Global variables in globvars.h needed for virtual grid */
double dglobx, dgloby, dglobz;
int xglobdim, yglobdim, zglobdim;


void setup_grid(double *xyzminmax, int *xyzdim, int *xyzind, 
                int order, double theta)
{
    /* local variables */
    int i;
    double t1;

    /* changing values of our extern variables */
    torder = order;
    torderlim = torder + 1;
    thetasq = theta * theta;
    torderflat = torderlim * (torderlim + 1) * (torderlim + 2) / 6;

    xglobdim = xyzdim[0];
    yglobdim = xyzdim[1];
    zglobdim = xyzdim[2];

    dglobx = (xyzminmax[1]-xyzminmax[0]) / (xglobdim - 1);
    dgloby = (xyzminmax[3]-xyzminmax[2]) / (yglobdim - 1);
    dglobz = (xyzminmax[5]-xyzminmax[4]) / (zglobdim - 1);


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


    xyzind[0] = 0;
    xyzind[1] = xglobdim - 1;
    xyzind[2] = 0;
    xyzind[3] = yglobdim - 1;
    xyzind[4] = 0;
    xyzind[5] = zglobdim - 1;

} /* END of function setup */




void cp_create_tree_n0_grid(struct tnode **p, int maxparnode, 
                            double *xyzmm, int *xyzdim, int *xyzind, int level)
{
    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t2;
    int i, j, loclev, numposchild;
    
    double xyzmms[6][8];
    int xyzdims[3][8];
    int xyzinds[6][8];

    double lxyzmm[6];
    int lxyzdim[3];
    int lxyzind[6];


    for (j = 0; j < 8; j++) {
        for (i = 0; i < 6; i++) xyzmms[i][j] = 0.0;
        for (i = 0; i < 3; i++) xyzdims[i][j] = 0;
        for (i = 0; i < 6; i++) xyzinds[i][j] = 0;
    }

    for (i = 0; i < 6; i++) lxyzmm[i] = 0.0;
    for (i = 0; i < 3; i++) lxyzdim[i] = 0;
    for (i = 0; i < 6; i++) lxyzind[i] = 0;

    (*p) = malloc(sizeof(struct tnode));


    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = xyzdim[0] * xyzdim[1] * xyzdim[2];
    (*p)->exist_ms = 0;

    (*p)->x_min = xyzmm[0];
    (*p)->x_max = xyzmm[1];
    (*p)->y_min = xyzmm[2];
    (*p)->y_max = xyzmm[3];
    (*p)->z_min = xyzmm[4];
    (*p)->z_max = xyzmm[5];

    (*p)->xdim = xyzdim[0];
    (*p)->ydim = xyzdim[1];
    (*p)->zdim = xyzdim[2];

    (*p)->xlind = xyzind[0];
    (*p)->xhind = xyzind[1];
    (*p)->ylind = xyzind[2];
    (*p)->yhind = xyzind[3];
    (*p)->zlind = xyzind[4];
    (*p)->zhind = xyzind[5];

    /*compute aspect ratio*/
    xl = (*p)->x_max - (*p)->x_min;
    yl = (*p)->y_max - (*p)->y_min;
    zl = (*p)->z_max - (*p)->z_min;
        
    lmax = max3(xl, yl, zl);
    t2 = min3(xl, yl, zl);

    if (t2 != 0.0) (*p)->aspect = lmax/t2;
    else (*p)->aspect = 0.0;

/*
    printf("\n\nNode at level %d: \n", level);
    printf("    bounds: %f, %f, %f, %f, %f, %f \n", xyzmm[0], xyzmm[1], xyzmm[2],
                                                    xyzmm[3], xyzmm[4], xyzmm[5]);
    printf("    dimensions: %d, %d, %d\n", xyzdim[0], xyzdim[1], xyzdim[2]);
    printf("    indices: %d, %d, %d, %d, %d, %d \n", xyzind[0], xyzind[1], xyzind[2],
                                                     xyzind[3], xyzind[4], xyzind[5]);
    printf("\n\n");
*/
    
    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    (*p)->sqradius = (xl*xl + yl*yl + zl*zl)/4.0;
    (*p)->radius = sqrt((*p)->sqradius);

    
    /*set tree level of node, and nullify child pointers*/
    (*p)->level = level;
    if (maxlevel < level) maxlevel = level;

    (*p)->num_children = 0;
    for (i = 0; i < 8; i++) (*p)->child[i] = NULL;
    

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

        xyzdims[0][0] = (*p)->xdim;
        xyzdims[1][0] = (*p)->ydim;
        xyzdims[2][0] = (*p)->zdim;

        xyzinds[0][0] = (*p)->xlind;
        xyzinds[1][0] = (*p)->xhind;
        xyzinds[2][0] = (*p)->ylind;
        xyzinds[3][0] = (*p)->yhind;
        xyzinds[4][0] = (*p)->zlind;
        xyzinds[5][0] = (*p)->zhind;

        x_mid = (*p)->x_mid;
        y_mid = (*p)->y_mid;
        z_mid = (*p)->z_mid;

        cp_partition_8_grid(xyzmms, xyzdims, xyzinds, xl, yl, zl, lmax, 
                            &numposchild, x_mid, y_mid, z_mid);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (xyzinds[0][i] <= xyzinds[1][i] &&
                xyzinds[2][i] <= xyzinds[3][i] &&
                xyzinds[4][i] <= xyzinds[5][i]) {

                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++) lxyzmm[j] = xyzmms[j][i];
                for (j = 0; j < 3; j++) lxyzdim[j] = xyzdims[j][i];
                for (j = 0; j < 6; j++) lxyzind[j] = xyzinds[j][i];

                //cp_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                //                  ind[i][0], ind[i][1], x, y, z, shrink,
                //                  maxparnode, lxyzmm, loclev);

                cp_create_tree_n0_grid(&((*p)->child[(*p)->num_children - 1]), 
                                       maxparnode, lxyzmm, lxyzdim, lxyzind, loclev);
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




void cp_partition_8_grid(double xyzmms[6][8], int xyzdims[3][8], int xyzinds[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid)
{

    /* local variables */
    int i, j;
    double critlen;
    int xdim, ydim, zdim, xn, yn, zn;
    int xlowind, xhighind, ylowind, yhighind, zlowind, zhighind;
    double xlowmid, xhighmid, ylowmid, yhighmid, zlowmid, zhighmid;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    xdim = xyzdims[0][0];
    ydim = xyzdims[1][0];
    zdim = xyzdims[2][0];

    xn = xdim / 2;
    yn = ydim / 2;
    zn = zdim / 2;

    if (xl >= critlen) {
        xlowmid = xyzmms[0][0] + (xn-1)*dglobx;
        xhighmid = xyzmms[1][0] - (xdim-xn-1)*dglobx;

        xlowind = xyzinds[0][0] + (xn-1);
        xhighind = xyzinds[1][0] - (xdim-xn-1);

        for (i = 0; i < 6; i++) xyzmms[i][1] = xyzmms[i][0];
        for (i = 0; i < 6; i++) xyzinds[i][1] = xyzinds[i][0];
        
        xyzmms[1][0] = xlowmid;
        xyzmms[0][1] = xhighmid;

        xyzinds[1][0] = xlowind;
        xyzinds[0][1] = xhighind;
        
        *numposchild = 2 * *numposchild;
        xdiv++;
    }

    if (yl >= critlen) {

        ylowmid = xyzmms[2][0] + (yn-1)*dgloby;
        yhighmid = xyzmms[3][0] - (ydim-yn-1)*dgloby;

        ylowind = xyzinds[2][0] + (yn-1);
        yhighind = xyzinds[3][0] - (ydim-yn-1);

        for (i = 0; i < *numposchild; i++) {
                        
            for (j = 0; j < 6; j++) xyzmms[j][*numposchild + i] = xyzmms[j][i];
            for (j = 0; j < 6; j++) xyzinds[j][*numposchild + i] = xyzinds[j][i];

            xyzmms[3][i] = ylowmid;
            xyzmms[2][*numposchild + i] = yhighmid;

            xyzinds[3][i] = ylowind;
            xyzinds[2][*numposchild + i] = yhighind;
        }
        
        *numposchild = 2 * *numposchild;
        ydiv++;
    }

    if (zl >= critlen) {

        zlowmid = xyzmms[4][0] + (zn-1)*dglobz;
        zhighmid = xyzmms[5][0] - (zdim-zn-1)*dglobz;

        zlowind = xyzinds[4][0] + (zn-1);
        zhighind = xyzinds[5][0] - (zdim-zn-1);

        for (i = 0; i < *numposchild; i++) {
                        
            for (j = 0; j < 6; j++) xyzmms[j][*numposchild + i] = xyzmms[j][i];
            for (j = 0; j < 6; j++) xyzinds[j][*numposchild + i] = xyzinds[j][i];

            xyzmms[5][i] = zlowmid;
            xyzmms[4][*numposchild + i] = zhighmid;

            xyzinds[5][i] = zlowind;
            xyzinds[4][*numposchild + i] = zhighind;
        }
        
        *numposchild = 2 * *numposchild;
        zdiv++;
    }

    for (i = 0; i < *numposchild; i++) {

        xyzdims[0][i] = xyzinds[1][i] - xyzinds[0][i] + 1;
        xyzdims[1][i] = xyzinds[3][i] - xyzinds[2][i] + 1;
        xyzdims[2][i] = xyzinds[5][i] - xyzinds[4][i] + 1;
    }

    return;

} /* END of function cp_partition_8 */




void cp_treecode_grid(struct tnode *p, double *xS, double *yS, double *zS, double *qS, 
                      double *tpeng, double *EnP, int numparsS, int numparsT,
                      double *timetree)
{
    /* local variables */
    int i, j;
    double time1, time2, time3;

    for (i = 0; i < numparsT; i++)
        EnP[i] = 0.0;

    time1 = MPI_Wtime();
    
    for (i = 0; i < numparsS; i++) {
        tarpos[0] = xS[i];
        tarpos[1] = yS[i];
        tarpos[2] = zS[i];
        tarposq = qS[i];

        for (j = 0; j < p->num_children; j++) {
            compute_cp1_grid(p->child[j], EnP);
        }
    }
    
    time2 = MPI_Wtime();

    compute_cp2_grid(p, EnP);
    
    time3 = MPI_Wtime();
    timetree[0] = time2 - time1;
    timetree[1] = time3 - time2;

    *tpeng = sum(EnP, numparsT);

    return;

} /* END of function cp_treecode */




void compute_cp1_grid(struct tnode *p, double *EnP)
{
    /* local variables */
    double tx, ty, tz, distsq;
    int i, j, k;


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
            cp_comp_direct_grid(p, EnP);
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_cp1_grid(p->child[i], EnP);
            }
        }
    }

    return;

} /* END of function compute_cp1 */




void compute_cp2_grid(struct tnode *ap, double *EnP)
{
    /* local variables */
    double tx, ty, peng;
    double xm, ym, zm, dx, dy, dz, xl, yl, zl;
    int xlind, ylind, zlind, xhind, yhind, zhind;
    int i, nn, j, k, k1, k2, k3, kk, porder, porder1;

    porder = torder;
    porder1 = porder - 1;

    if (ap->exist_ms == 1) {
        xm = ap->x_mid;
        ym = ap->y_mid;
        zm = ap->z_mid;

        xl = ap->x_min;
        yl = ap->y_min;
        zl = ap->z_min;

        xlind = ap->xlind;
        ylind = ap->ylind;
        zlind = ap->zlind;

        xhind = ap->xhind;
        yhind = ap->yhind;
        zhind = ap->zhind;

        for (i = xlind; i <= xhind; i++) {
            dx = xl + (i-xlind)*dglobx  - xm;

            for (j = ylind; j <= yhind; j++) {
                dy = yl + (j-ylind)*dgloby  - ym;

                for (k = zlind; k <= zhind; k++) {

                    nn = (i * yglobdim * zglobdim) + (j * zglobdim) + k;
                    dz = zl + (k-zlind)*dglobz  - zm;

                    kk = 0;

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

                    EnP[nn] = EnP[nn] + peng;
                }
            }
        }
    }

    for (j = 0; j < ap->num_children; j++)
        compute_cp2_grid(ap->child[j], EnP);

    return;

} /* END function compute_cp2 */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void cp_comp_direct_grid(struct tnode *ap, double *EnP)
{

    /* local variables */
    int i, j, k, nn;
    double tx, ty, tz, xl, yl, zl;
    double xlind, ylind, zlind, xhind, yhind, zhind;

    xl = ap->x_min;
    yl = ap->y_min;
    zl = ap->z_min;

    xlind = ap->xlind;
    ylind = ap->ylind;
    zlind = ap->zlind;

    xhind = ap->xhind;
    yhind = ap->yhind;
    zhind = ap->zhind;

    for (i = xlind; i <= xhind; i++) {
        tx = xl + (i-xlind)*dglobx  - tarpos[0];

        for (j = ylind; j <= yhind; j++) {
            ty = yl + (j-ylind)*dgloby  - tarpos[1];

            for (k = zlind; k <= zhind; k++) {

                nn = (i * yglobdim * zglobdim) + (j * zglobdim) + k;
                tz = zl + (k-zlind)*dglobz  - tarpos[2];

                EnP[nn] = EnP[nn] + tarposq / sqrt(tx*tx + ty*ty + tz*tz);
            }
        }
    }

    return;

} /* END function cp_comp_direct */
