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


void pc_treecode_yuk(struct tnode *p, struct batch *batches,
                     double *xS, double *yS, double *zS,
                     double *qS, double *xT, double *yT, double *zT,
                     int numparsS, int numparsT, double kappa,
                     double *tpeng, double *EnP)
{
    /* local variables */
    int i, j;
    
    for (i = 0; i < numparsT; i++)
        EnP[i] = 0.0;
    
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc_yuk(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                xS, yS, zS, qS, xT, yT, zT, kappa, EnP);
        }
    }
    
    *tpeng = sum(EnP, numparsT);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_yuk(struct tnode *p,
                    int *batch_ind, double *batch_mid, double batch_rad,
                    double *xS, double *yS, double *zS, double *qS,
                    double *xT, double *yT, double *zT,
                    double kappa, double *EnP)
{
    /* local variables */
    double tx, ty, tz, dist;
    int i, j, k, kk, ii;

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
        
        if (p->exist_ms == 0) {
            make_vector(p->ms, torderflat);
            
            for (i = 0; i < torderflat; i++) {
                p->ms[i] = 0.0;
            }
            
            pc_comp_ms(p, xS, yS, zS, qS);
            p->exist_ms = 1;
        }
        
        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
            tx = xT[ii] - p->x_mid;
            ty = yT[ii] - p->y_mid;
            tz = zT[ii] - p->z_mid;
        
            comp_tcoeff_yuk(tx, ty, tz, kappa);
        
            kk = -1;
        
            for (k = torder; k > -1; k--) {
                for (j = torder - k; j > -1; j--) {
                    for (i = torder - k - j; i > -1; i--) {
                        EnP[ii] += a1[i+2][j+2][k+2] * p->ms[++kk];
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
            pc_comp_direct_yuk(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, xT, yT, zT, kappa, EnP);
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_pc_yuk(p->child[i], batch_ind, batch_mid, batch_rad,
                           xS, yS, zS, qS, xT, yT, zT, kappa, EnP);
            }
        }
    }
    
    return;
    
} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current cluster due
 * to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct_yuk(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *restrict xS, double *restrict yS, double *restrict zS, double *restrict qS,
                    double *restrict xT, double *restrict yT, double *restrict zT,
                    double kappa, double *EnP)
{
    
    /* local variables */
    int i, ii;
    double tx, ty, tz, dist;
    
    double d_peng;
    
    #pragma acc data present(xS, yS, zS, qS)
    #pragma acc region
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
        d_peng = 0.0;
        for (i = ibeg - 1; i < iend; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            
            dist = sqrt(tx*tx + ty*ty + tz*tz);
            d_peng += qS[i] * exp(-kappa * dist) / dist;
        }
        EnP[ii] += d_peng;
    }
    
    return;
    
} /* END function pc_comp_direct */
