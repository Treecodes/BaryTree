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


void pc_treecode_yuk(struct tnode *p, double *xS, double *yS, double *zS,
                     double *qS, double *xT, double *yT, double *zT,
                     double *tpeng, double *EnP,
                     int numparsS, int numparsT, double kappa)
{
    /* local variables */
    int i, j;
    double penglocal, peng;
    
    for (i = 0; i < numparsT; i++)
        EnP[i] = 0.0;
    
    for (i = 0; i < numparsT; i++) {
        peng = 0.0;
        tarpos[0] = xT[i];
        tarpos[1] = yT[i];
        tarpos[2] = zT[i];
        
        for (j = 0; j < p->num_children; j++) {
            compute_pc_yuk(p->child[j], &penglocal, xS, yS, zS, qS, kappa);
            peng += penglocal;
        }
        
        EnP[i] = peng;
    }
    
    *tpeng = sum(EnP, numparsT);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_yuk(struct tnode *p, double *peng,
                    double *x, double *y, double *z, double *q,
                    double kappa)
{
    /* local variables */
    double tx, ty, tz, distsq, penglocal;
    int i, j, k;

    
    /* determine DISTSQ for MAC test */
    tx = tarpos[0] - p->x_mid;
    ty = tarpos[1] - p->y_mid;
    tz = tarpos[2] - p->z_mid;
    distsq = tx*tx + ty*ty + tz*tz;
    
    *peng = 0.0;
    
    if ((p->sqradius < distsq * thetasq) && (p->sqradius != 0.00)) {
        /*
         * If MAC is accepted and there is more than 1 particle
         * in the box, use the expansion for the approximation.
         */
        
        if (p->exist_ms == 0) {
            make_3array(p->ms, torder+1, torder+1, torder+1);
            
            for (i = 0; i < torder + 1; i++) {
                for (j = 0; j < torder + 1; j++) {
                    for (k = 0; k < torder + 1; k++) {
                        p->ms[i][j][k] = 0.0;
                    }
                }
            }
            
            pc_comp_ms(p, x, y, z, q);
            p->exist_ms = 1;
        }
        
        comp_tcoeff_yuk(tx, ty, tz, kappa);
        
        for (k = 0; k < torder+1; k++) {
            for (j = 0; j < torder-k+1; j++) {
                for (i = 0; i < torder-k-j+1; i++) {
                    *peng += a1[i+2][j+2][k+2] * p->ms[i][j][k];
                }
            }
        }
        
    } else {
        
        /*
         * If MAC fails check to see if there are children. If not, perform direct
         * calculation. If there are children, call routine recursively for each.
         */
        
        if (p->num_children == 0) {
            pc_comp_direct_yuk(&penglocal, p->ibeg, p->iend, x, y, z, q, kappa);
            *peng = penglocal;
        } else {
            for (i = 0; i < p->num_children; i++) {
                compute_pc_yuk(p->child[i], &penglocal, x, y, z, q, kappa);
                *peng += penglocal;
            }
        }
    }
    
    return;
    
} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current cluster due
 * to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct_yuk(double *peng, int ibeg, int iend,
                        double *x, double *y, double *z, double *q, double kappa)
{
    
    /* local variables */
    int i;
    double tx, ty, tz, dist;
    
    *peng = 0.0;
    
    for (i = ibeg - 1; i < iend; i++) {
        tx = x[i] - tarpos[0];
        ty = y[i] - tarpos[1];
        tz = z[i] - tarpos[2];
        dist = sqrt(tx*tx + ty*ty + tz*tz);
        
        *peng += q[i] * exp(-kappa * dist) / dist;
    }
    
    return;
    
} /* END function pc_comp_direct */
