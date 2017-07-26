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

#include "partition.h"
#include "tree.h"


void setup_yuk(double *x, double *y, double *z,
               int numpars, int order, double theta, 
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
    make_vector(cf3, torderlim);

    make_3array(a1, torderlim+3, torderlim+3, torderlim+3);
    make_3array(b1, torderlim+3, torderlim+3, torderlim+3);


    /* initializing arrays for Taylor sums and coefficients */
    for (i = 0; i < torder + 1; i++)
        cf[i] = i + 1.0;

    for (i = 0; i < torderlim; i++) {
        t1 = 1.0 / (i + 1.0);
        cf1[i] = t1;
        cf2[i] = 1.0 - (0.5 * t1);
        cf3[i] = 1.0 - t1;
    }

    /* find bounds of Cartesian box enclosing the particles */
    xyzminmax[0] = minval(x, numpars);
    xyzminmax[1] = maxval(x, numpars);
    xyzminmax[2] = minval(y, numpars);
    xyzminmax[3] = maxval(y, numpars);
    xyzminmax[4] = minval(z, numpars);
    xyzminmax[5] = maxval(z, numpars);

    make_vector(orderarr, numpars);

    for (i = 0; i < numpars; i++)
        orderarr[i] = i+1;

    return;
    
} /* END of function setup_yuk */




void cp_treecode_yuk(struct tnode *p, double *xS, double *yS, double *zS,
                     double *qS, double *xT, double *yT, double *zT,
                     double *tpeng, double *EnP, int numparsS, int numparsT,
                     double kappa, double *timetree)
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

        for (j = 0; j < p->num_children; j++)
            compute_cp1_yuk(p->child[j], EnP, xT, yT, zT, kappa);
    }

    time2 = MPI_Wtime();
    
    compute_cp2(p, xT, yT, zT, EnP);
    
    time3 = MPI_Wtime();
    timetree[0] = time2 - time1;
    timetree[1] = time3 - time2;

    *tpeng = sum(EnP, numparsT);

    return;

} /* END of function cp_treecode */




void compute_cp1_yuk(struct tnode *p, double *EnP,
                     double *x, double *y, double *z, double kappa)
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
        for (i = 0; i < torderlim + 1; i++) {
            for (j = 0; j < torderlim + 1; j++) {
                for (k = 0; k < torderlim + 1; k++) {
                    b1[i][j][k] = 0.0;
                    a1[i][j][k] = 0.0;
                }
            }
        }

        comp_tcoeff_yuk(tx, ty, tz, kappa);

        if (p->exist_ms == 0) {
            make_vector(p->ms, torderflat);

            for (i = 0; i < torderflat; i++) {
                p->ms[i] = 0.0;
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
            cp_comp_direct_yuk(EnP, p->ibeg, p->iend, x, y, z, kappa);
        } else {
            for (i = 0; i < p->num_children; i++)
                compute_cp1_yuk(p->child[i], EnP, x, y, z, kappa);
        }
    }

    return;

} /* END of function compute_cp1 */




void comp_tcoeff_yuk(double dx, double dy, double dz, double kappa)
{
    /* local variables */
    double tdx, tdy, tdz, fac, dist, dist2;
    int i, j, k, i1, i2, j1, j2, k1, k2;

    double kappax, kappay, kappaz;

    /* setup variables */
    tdx = 2.0 * dx;
    tdy = 2.0 * dy;
    tdz = 2.0 * dz;

    kappax = kappa * dx;
    kappay = kappa * dy;
    kappaz = kappa * dz;

    dist2 = dx*dx + dy*dy + dz*dz;
    fac = 1.0 / (dist2);
    dist = sqrt(dist2);

    /* 0th coefficient or function val */
    a1[0][0][0] = exp(-kappa * dist);
    b1[0][0][0] = a1[0][0][0] / dist;

    /* set of indices for which two of them are 0 */
    a1[1][0][0] = kappax * b1[0][0][0];
    a1[0][1][0] = kappay * b1[0][0][0];
    a1[0][0][1] = kappaz * b1[0][0][0];

    b1[1][0][0] = fac * dx * (b1[0][0][0] + kappa*a1[0][0][0]);
    b1[0][1][0] = fac * dy * (b1[0][0][0] + kappa*a1[0][0][0]);
    b1[0][0][1] = fac * dz * (b1[0][0][0] + kappa*a1[0][0][0]);

    for (i = 2; i < torderlim+1; i++) {
        i1 = i - 1;
        i2 = i - 2;

        a1[i][0][0] = cf1[i1] * kappa * (dx * b1[i1][0][0] - b1[i2][0][0]);
        a1[0][i][0] = cf1[i1] * kappa * (dy * b1[0][i1][0] - b1[0][i2][0]);
        a1[0][0][i] = cf1[i1] * kappa * (dz * b1[0][0][i2] - b1[0][0][i2]);

        b1[i][0][0] = fac * (tdx * cf2[i1] * b1[i1][0][0]
                                 - cf3[i1] * b1[i2][0][0]
                   + cf1[i1] * kappa * (dx * a1[i1][0][0]
                                           - a1[i2][0][0]));

        b1[0][i][0] = fac * (tdy * cf2[i1] * b1[0][i1][0]
                                 - cf3[i1] * b1[0][i2][0]
                   + cf1[i1] * kappa * (dy * a1[0][i1][0]
                                           - a1[0][i2][0]));

        b1[0][0][i] = fac * (tdz * cf2[i1] * b1[0][0][i1]
                                 - cf3[i1] * b1[0][0][i2]
                   + cf1[i1] * kappa * (dz * a1[0][0][i1]
                                           - a1[0][0][i2]));
    }

    /* set of indices for which one is 0, one is 1, and other is >= 1 */
    a1[1][1][0] = kappax * b1[0][1][0];
    a1[1][0][1] = kappay * b1[0][0][1];
    a1[0][1][1] = kappaz * b1[0][0][1];


    b1[1][1][0] = fac * (dx * b1[0][1][0] + tdy * b1[1][0][0]
                                       + kappax * a1[0][1][0]);
    b1[1][0][1] = fac * (dx * b1[0][0][1] + tdz * b1[1][0][0]
                                       + kappax * a1[0][0][1]);
    b1[0][1][1] = fac * (dy * b1[0][0][1] + tdz * b1[0][1][0]
                                       + kappay * a1[0][0][1]);

    for (i = 2; i < torderlim; i++) {
        i1 = i - 1;
        i2 = i - 2;

        a1[1][0][i] = kappax * b1[0][0][i];
        a1[0][1][i] = kappay * b1[0][0][i];
        a1[0][i][1] = kappaz * b1[0][i][0];
        a1[1][i][0] = kappax * b1[0][i][0];
        a1[i][1][0] = kappay * b1[i][0][0];
        a1[i][0][1] = kappaz * b1[i][0][0];

        b1[1][0][i] = fac * (dx * b1[0][0][i] + tdz * b1[1][0][i1]
                                                    - b1[1][0][i2]
                                           + kappax * a1[0][0][i]);

        b1[0][1][i] = fac * (dy * b1[0][0][i] + tdz * b1[0][1][i1]
                                - b1[0][1][i2]
                                           + kappay * a1[0][0][i]);

        b1[0][i][1] = fac * (dz * b1[0][i][0] + tdy * b1[0][i1][1]
                                - b1[0][i2][1]
                                           + kappaz * a1[0][i][0]);

        b1[1][i][0] = fac * (dx * b1[0][i][0] + tdy * b1[1][i1][0]
                                - b1[1][i2][0]
                                           + kappax * a1[0][i][0]);

        b1[i][1][0] = fac * (dy * b1[i][0][0] + tdx * b1[i1][1][0]
                                - b1[i2][1][0]
                                           + kappay * a1[i][0][0]);

        b1[i][0][1] = fac * (dz * b1[i][0][0] + tdx * b1[i1][0][1]
                                - b1[i2][0][1]
                                           + kappaz * a1[i][0][0]);
    }

    /* set of indices for which one is 0, others are >=2 */
    for (i = 2; i < torderlim - 1; i++) {
        i1 = i - 1;
        i2 = i - 2;

        for (j = 2; j < torderlim - i + 1; j++) {
            j1 = j - 1;
            j2 = j - 2;

            a1[i][j][0] = cf1[i1] * kappa * (dx * b1[i1][j][0]
                                                - b1[i2][j][0]);

            a1[i][0][j] = cf1[i1] * kappa * (dx * b1[i1][0][j]
                                                - b1[i2][0][j]);

            a1[0][i][j] = cf1[i1] * kappa * (dy * b1[0][i1][j]
                                                - b1[0][i2][j]);
                        
            b1[i][j][0] = fac * (tdx * cf2[i1] * b1[i1][j][0]
                                         + tdy * b1[i][j1][0]
                                     - cf3[i1] * b1[i2][j][0]
                                               - b1[i][j2][0]
                       + cf1[i1] * kappa * (dx * a1[i1][j][0]
                                               - a1[i2][j][0]));

            b1[i][0][j] = fac * (tdx * cf2[i1] * b1[i1][0][j]
                                         + tdz * b1[i][0][j1]
                                     - cf3[i1] * b1[i2][0][j]
                                               - b1[i][0][j2]
                       + cf1[i1] * kappa * (dx * a1[i1][0][j]
                                               - a1[i2][0][j]));

            b1[0][i][j] = fac * (tdy * cf2[i1] * b1[0][i1][j]
                                         + tdz * b1[0][i][j1]
                                     - cf3[i1] * b1[0][i2][j]
                                               - b1[0][i][j2]
                       + cf1[i1] * kappa * (dy * a1[0][i1][j]
                                               - a1[0][i2][j]));
        }
    }

    /* set of indices for which two are 1, other is >= 1 */
    a1[1][1][1] = kappax * b1[0][1][1];
    b1[1][1][1] = fac * (dx * b1[0][1][1] + tdy * b1[1][0][1]
                                          + tdz * b1[1][1][0]
                                       + kappax * a1[0][1][1]);

    for (i = 2; i < torderlim - 1; i++) {
        i1 = i - 1;
        i2 = i - 2;

        a1[1][1][i] = kappax * b1[0][1][i];
        a1[1][i][1] = kappax * b1[0][i][1];
        a1[i][1][1] = kappay * b1[i][0][1];

        b1[1][1][i] = fac * (dx * b1[0][1][i] + tdy * b1[1][0][i]
                          + tdz * b1[1][1][i1]      - b1[1][1][i2]
                                          + kappax * a1[0][1][i]);

        b1[1][i][1] = fac * (dx * b1[0][i][1] + tdy * b1[1][i1][1]
                          + tdz * b1[1][i][0]       - b1[1][i2][1]
                                           + kappax * a1[0][i][1]);

        b1[i][1][1] = fac * (dy * b1[i][0][1] + tdx * b1[i1][1][1]
                          + tdz * b1[i][1][0]       - b1[i2][1][1]
                                           + kappay * a1[i][0][1]);
    }

    /* set of indices for which one is 1, others are >= 2 */
    for (i = 2; i < torderlim - 2; i++) {
        i1 = i - 1;
        i2 = i - 2;

        for (j = 2; j < torderlim - i + 1; j++) {
            j1 = j - 1;
            j2 = j - 2;

            a1[1][i][j] = kappax * b1[0][i][j];
            a1[i][1][j] = kappay * b1[i][0][j];
            a1[i][j][1] = kappaz * b1[i][j][0];

            b1[1][i][j] = fac * (dx * b1[0][i][j] + tdy * b1[1][i1][j]
                              + tdz * b1[1][i][j1]      - b1[1][i2][j]
                                    - b1[1][i][j2]
                                               + kappax * a1[0][i][j]);

            b1[i][1][j] = fac * (dy * b1[i][0][j] + tdx * b1[i1][1][j]
                              + tdz * b1[i][1][j1]      - b1[i2][1][j]
                                    - b1[i][1][j2]
                                               + kappay * a1[i][0][j]);

            b1[i][j][1] = fac * (dz * b1[i][j][0] + tdx * b1[i1][j][1]
                              + tdy * b1[i][j1][1]      - b1[i2][j][1]
                                    - b1[i][j2][1]
                                               + kappaz * a1[i][j][0]);
        }
    }

    /* set of indices for which all are >= 2 */
    for (k = 2; k < torderlim - 3; k++) {
        k1 = k - 1;
        k2 = k - 2;
                
        for (j = 2; j < torderlim - k - 1; j++) {
            j1 = j - 1;
            j2 = j - 2;

            for (i = 2; i < torderlim - k - j + 1; i++) {
                i1 = i - 1;
                i2 = i - 2;

                a1[i][j][k] = cf1[i1] * kappa * (dx * b1[i1][j][k]
                                                    - b1[i2][j][k]);

                b1[i][j][k] = fac * (tdx * cf2[i1] * b1[i1][j][k]
                                             + tdy * b1[i][j1][k]
                                             + tdz * b1[i][j][k1]
                                         - cf3[i1] * b1[i2][j][k]
                                    - b1[i][j2][k] - b1[i][j][k2]
                           + cf1[i1] * kappa * (dx * a1[i1][j][k]
                                                   - a1[i2][j][k]));
            }
        }
    }

    return;

} /* END function comp_tcoeff */




/*
 * comp_direct directly computes the potential on the targets in the current cluster due
 * to the current source, determined by the global variable TARPOS
 */
void cp_comp_direct_yuk(double *EnP, int ibeg, int iend,
                        double *x, double *y, double *z, double kappa)
{

    /* local variables */
    int i, nn;
    double tx, ty, tz, dist;

    for (i = ibeg - 1; i < iend; i++) {
        nn = orderarr[i];
        tx = x[i] - tarpos[0];
        ty = y[i] - tarpos[1];
        tz = z[i] - tarpos[2];
        dist = sqrt(tx*tx + ty*ty + tz*tz);

        EnP[nn-1] = EnP[nn-1] + tarposq * exp(-kappa * dist) / dist;
    }

    return;

} /* END function comp_direct */
