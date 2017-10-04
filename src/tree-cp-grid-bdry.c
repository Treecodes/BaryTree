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


void setup_grid_bdry(double *xyzminmax, int *xyzdim, int *xyzind,
                int order, double theta)
{
    /* local variables */
    int i;
    double t1;

    /* changing values of our extern variables */
    torder = order;
    torderlim = torder + 1;
    thetasq = theta * theta;

    xglobdim = xyzdim[0];
    yglobdim = xyzdim[1];

    dglobx = (xyzminmax[1]-xyzminmax[0]) / (xglobdim - 1);
    dgloby = (xyzminmax[3]-xyzminmax[2]) / (yglobdim - 1);


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

} /* END of function setup */




void cp_create_tree_n0_grid_bdry(struct tnode **p, int maxparnode,
                            double *xyzmm, int *xyzdim, int *xyzind, int level)
{
    /*local variables*/
    double x_mid, y_mid, xl, yl, lmax, t2;
    int i, j, loclev, numposchild;
    
    double xyzmms[4][4];
    int xyzdims[2][4];
    int xyzinds[4][4];

    double lxyzmm[4];
    int lxyzdim[2];
    int lxyzind[4];


    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) xyzmms[i][j] = 0.0;
        for (i = 0; i < 2; i++) xyzdims[i][j] = 0;
        for (i = 0; i < 4; i++) xyzinds[i][j] = 0;
    }

    for (i = 0; i < 4; i++) lxyzmm[i] = 0.0;
    for (i = 0; i < 2; i++) lxyzdim[i] = 0;
    for (i = 0; i < 4; i++) lxyzind[i] = 0;

    (*p) = malloc(sizeof(struct tnode));


    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = xyzdim[0] * xyzdim[1];
    (*p)->exist_ms = 0;

    (*p)->x_min = xyzmm[0];
    (*p)->x_max = xyzmm[1];
    (*p)->y_min = xyzmm[2];
    (*p)->y_max = xyzmm[3];

    (*p)->xdim = xyzdim[0];
    (*p)->ydim = xyzdim[1];

    (*p)->xlind = xyzind[0];
    (*p)->xhind = xyzind[1];
    (*p)->ylind = xyzind[2];
    (*p)->yhind = xyzind[3];

    /*compute aspect ratio*/
    xl = (*p)->x_max - (*p)->x_min;
    yl = (*p)->y_max - (*p)->y_min;
        
    lmax = fmax(xl, yl);
    t2 = fmin(xl, yl);

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

    (*p)->sqradius = (xl*xl + yl*yl)/4.0;
    (*p)->radius = sqrt((*p)->sqradius);

    
    /*set tree level of node, and nullify child pointers*/
    (*p)->level = level;
    if (maxlevel < level) maxlevel = level;

    (*p)->num_children = 0;
    for (i = 0; i < 4; i++) (*p)->child[i] = NULL;
    

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

        xyzdims[0][0] = (*p)->xdim;
        xyzdims[1][0] = (*p)->ydim;

        xyzinds[0][0] = (*p)->xlind;
        xyzinds[1][0] = (*p)->xhind;
        xyzinds[2][0] = (*p)->ylind;
        xyzinds[3][0] = (*p)->yhind;

        x_mid = (*p)->x_mid;
        y_mid = (*p)->y_mid;

        cp_partition_8_grid_bdry(xyzmms, xyzdims, xyzinds, xl, yl, lmax,
                            &numposchild, x_mid, y_mid);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (xyzinds[0][i] <= xyzinds[1][i] &&
                xyzinds[2][i] <= xyzinds[3][i]) {

                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 4; j++) lxyzmm[j] = xyzmms[j][i];
                for (j = 0; j < 2; j++) lxyzdim[j] = xyzdims[j][i];
                for (j = 0; j < 4; j++) lxyzind[j] = xyzinds[j][i];

                cp_create_tree_n0_grid_bdry(&((*p)->child[(*p)->num_children - 1]),
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




void cp_partition_8_grid_bdry(double xyzmms[4][4], int xyzdims[2][4], int xyzinds[4][4],
                    double xl, double yl, double lmax, int *numposchild,
                    double x_mid, double y_mid)
{

    /* local variables */
    int i, j;
    double critlen;
    int xdim, ydim, xn, yn;
    int xlowind, xhighind, ylowind, yhighind;
    double xlowmid, xhighmid, ylowmid, yhighmid;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    xdim = xyzdims[0][0];
    ydim = xyzdims[1][0];

    xn = xdim / 2;
    yn = ydim / 2;

    if (xl >= critlen) {
        xlowmid = xyzmms[0][0] + (xn-1)*dglobx;
        xhighmid = xyzmms[1][0] - (xdim-xn-1)*dglobx;

        xlowind = xyzinds[0][0] + (xn-1);
        xhighind = xyzinds[1][0] - (xdim-xn-1);

        for (i = 0; i < 4; i++) xyzmms[i][1] = xyzmms[i][0];
        for (i = 0; i < 4; i++) xyzinds[i][1] = xyzinds[i][0];
        
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
                        
            for (j = 0; j < 4; j++) xyzmms[j][*numposchild + i] = xyzmms[j][i];
            for (j = 0; j < 4; j++) xyzinds[j][*numposchild + i] = xyzinds[j][i];

            xyzmms[3][i] = ylowmid;
            xyzmms[2][*numposchild + i] = yhighmid;

            xyzinds[3][i] = ylowind;
            xyzinds[2][*numposchild + i] = yhighind;
        }
        
        *numposchild = 2 * *numposchild;
        ydiv++;
    }

    for (i = 0; i < *numposchild; i++) {
        xyzdims[0][i] = xyzinds[1][i] - xyzinds[0][i] + 1;
        xyzdims[1][i] = xyzinds[3][i] - xyzinds[2][i] + 1;
    }

    return;

} /* END of function cp_partition_8 */




void cp_treecode_grid_bdry(struct tnode *p, double zyx, int dir, double *xS, double *yS,
                           double *zS, double *qS, double *tpeng, double *EnP,
                           int numparsS, int numparsT, double *timetree)
{
    /* local variables */
    int i, j;
    double time1, time2, time3;

    for (i = 0; i < numparsT; i++)
        EnP[i] = 0.0;
    
    time1 = MPI_Wtime();
    
    if (p->num_children > 0) {
        for (i = 0; i < numparsS; i++) {
            tarpos[0] = xS[i];
            tarpos[1] = yS[i];
            tarpos[2] = zS[i];
            tarposq = qS[i];
            
            for (j = 0; j < p->num_children; j++) {
                compute_cp1_grid_bdry(p->child[j], zyx, dir, EnP);
            }
        }
    } else if (p->num_children == 0) {
        for (i = 0; i < numparsS; i++) {
            tarpos[0] = xS[i];
            tarpos[1] = yS[i];
            tarpos[2] = zS[i];
            tarposq = qS[i];
            
            cp_comp_direct_grid_bdry(p, zyx, dir, EnP);
        }
    }
    
    time2 = MPI_Wtime();

    compute_cp2_grid_bdry(p, zyx, dir, EnP);
    
    time3 = MPI_Wtime();
    timetree[0] = time2 - time1;
    timetree[1] = time3 - time2;

    *tpeng = sum(EnP, numparsT);

    return;

} /* END of function cp_treecode */




void compute_cp1_grid_bdry(struct tnode *p, double zyx, int dir, double *EnP)
{
    /* local variables */
    double tx=0, ty=0, tz=0, distsq;
    int i, j, k;


    /* determine DISTSQ for MAC test */
    if (dir == 0) {
        tx = tarpos[0] - p->x_mid;
        ty = tarpos[1] - p->y_mid;
        tz = tarpos[2] - zyx;
    } else if (dir == 1) {
        tx = tarpos[0] - p->x_mid;
        ty = tarpos[1] - zyx;
        tz = tarpos[2] - p->y_mid;
    } else if (dir == 2) {
        tx = tarpos[0] - zyx;
        ty = tarpos[1] - p->x_mid;
        tz = tarpos[2] - p->y_mid;
    }
    distsq = tx*tx + ty*ty + tz*tz;

    if ((p->sqradius < distsq * thetasq) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */
        for (i = 0; i < torderlim; i++) {
            for (j = 0; j < torderlim; j++) {
                for (k = 0; k < torderlim; k++) {
                    b1[i][j][k] = 0.0;
                }
            }
        }
        
        comp_tcoeff(tx, ty, tz);

        if (p->exist_ms == 0) {
            make_3array(p->ms, torder+1, torder+1, torder+1);

            for (i = 0; i < torder + 1; i++) {
                for (j = 0; j < torder + 1; j++) {
                    for (k = 0; k < torder + 1; k++) {
                        p->ms[i][j][k] = 0.0;
                    }
                }
            }
            
            p->exist_ms = 1;
        }

        cp_comp_ms(p);

    } else {
        
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            cp_comp_direct_grid_bdry(p, zyx, dir, EnP);
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_cp1_grid_bdry(p->child[i], zyx, dir, EnP);
            }
        }
    }

    return;

} /* END of function compute_cp1 */




void compute_cp2_grid_bdry(struct tnode *ap, double zyx, int dir, double *EnP)
{
    /* local variables */
    double tx=0, ty=0, peng=0;
    double xm=0, ym=0, zm=0, dx=0, dy=0, dz=0, xl=0, yl=0, zl=0;
    int xlind=0, ylind=0, zlind=0, xhind=0, yhind=0, zhind=0;
    int ydim=0, zdim=0;
    double ddx=0, ddy=0, ddz=0;
    int i, nn, j, k, k1, k2, k3, porder, porder1;

    porder = torder;
    porder1 = porder - 1;

    if (ap->exist_ms == 1) {
        if (dir == 0) {
            xm = ap->x_mid;  ym = ap->y_mid;  zm = zyx;
            xl = ap->x_min;  yl = ap->y_min;  zl = zyx;

            xlind = ap->xlind;  ylind = ap->ylind;  zlind = 0;
            xhind = ap->xhind;  yhind = ap->yhind;  zhind = 0;
            
            ydim = yglobdim;  zdim = 1;
            ddx = dglobx;  ddy = dgloby;  ddz = 0;
            
        } else if (dir == 1) {
            xm = ap->x_mid;  ym = zyx;  zm = ap->y_mid;
            xl = ap->x_min;  yl = zyx;  zl = ap->y_min;
            
            xlind = ap->xlind;  ylind = 0;  zlind = ap->ylind;
            xhind = ap->xhind;  yhind = 0;  zhind = ap->yhind;
            
            ydim = 1;  zdim = yglobdim;
            ddx = dglobx;  ddy = 0;  ddz = dgloby;
            
        } else if (dir == 2) {
            xm = zyx;  ym = ap->x_mid;  zm = ap->y_mid;
            xl = zyx;  yl = ap->x_min;  zl = ap->y_min;
            
            xlind = 0;  ylind = ap->xlind;  zlind = ap->ylind;
            xhind = 0;  yhind = ap->xhind;  zhind = ap->yhind;
            
            ydim = xglobdim;  zdim = yglobdim;
            ddx = 0;  ddy = dglobx;  ddz = dgloby;
        }

        for (i = xlind; i <= xhind; i++) {
            dx = xl + (i-xlind)*ddx  - xm;
            for (j = ylind; j <= yhind; j++) {
                dy = yl + (j-ylind)*ddy  - ym;
                for (k = zlind; k <= zhind; k++) {

                    nn = (i * ydim * zdim) + (j * zdim) + k;
                    dz = zl + (k-zlind)*ddz  - zm;

                    peng = ap->ms[0][0][porder];

                    for (k3 = porder1; k3 > -1; k3--) {
                        ty = ap->ms[0][porder - k3][k3];

                        for (k2 = porder1 - k3; k2 > -1; k2--) {
                            tx = ap->ms[porder - k3 - k2][k2][k3];

                            for (k1 = porder1 - k3 - k2; k1 > -1; k1--) {
                                tx = dx*tx + ap->ms[k1][k2][k3];
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
        compute_cp2_grid_bdry(ap->child[j], zyx, dir, EnP);

    return;

} /* END function compute_cp2 */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void cp_comp_direct_grid_bdry(struct tnode *ap, double zyx, int dir, double *EnP)
{

    /* local variables */
    int i, j, k, nn;
    double tx=0, ty=0, tz=0, xl=0, yl=0, zl=0;
    double xlind=0, ylind=0, zlind=0, xhind=0, yhind=0, zhind=0;
    int ydim=0, zdim=0;
    double ddx=0, ddy=0, ddz=0;

    
    if (dir == 0) {
        xl = ap->x_min;  yl = ap->y_min;  zl = zyx;
        
        xlind = ap->xlind;  ylind = ap->ylind;  zlind = 0;
        xhind = ap->xhind;  yhind = ap->yhind;  zhind = 0;
        
        ydim = yglobdim;  zdim = 1;
        ddx = dglobx;  ddy = dgloby;  ddz = 0;
        
    } else if (dir == 1) {
        xl = ap->x_min;  yl = zyx;  zl = ap->y_min;
        
        xlind = ap->xlind;  ylind = 0;  zlind = ap->ylind;
        xhind = ap->xhind;  yhind = 0;  zhind = ap->yhind;
        
        ydim = 1;  zdim = yglobdim;
        ddx = dglobx;  ddy = 0;  ddz = dgloby;
        
    } else if (dir == 2) {
        xl = zyx;  yl = ap->x_min;  zl = ap->y_min;
        
        xlind = 0;  ylind = ap->xlind;  zlind = ap->ylind;
        xhind = 0;  yhind = ap->xhind;  zhind = ap->yhind;
        
        ydim = xglobdim;  zdim = yglobdim;
        ddx = 0;  ddy = dglobx;  ddz = dgloby;
    }
    
    
    for (i = xlind; i <= xhind; i++) {
        tx = xl + (i-xlind)*ddx  - tarpos[0];
        for (j = ylind; j <= yhind; j++) {
            ty = yl + (j-ylind)*ddy  - tarpos[1];
            for (k = zlind; k <= zhind; k++) {

                nn = (i * ydim * zdim) + (j * zdim) + k;
                tz = zl + (k-zlind)*ddz  - tarpos[2];

                EnP[nn] = EnP[nn] + tarposq / sqrt(tx*tx + ty*ty + tz*tz);
            }
        }
    }

    return;

} /* END function cp_comp_direct */
