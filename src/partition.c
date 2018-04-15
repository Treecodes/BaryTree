#include <stdio.h>

#include "tools.h"
#include "partition.h"
/* 
 * definition of partition functions
 *
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */

void cp_partition(double *a, double *b, double *c, int *indarr,
                  int ibeg, int iend, double val, int *midind)
{
        /* local variables */
        double ta, tb, tc;
        int lower, upper, tind;

        if (ibeg < iend) {
        /* 
         * temporarily stores ibeg entries and set a(ibeg) = val
         * for the partitioning algorithm.
         */
         
                ta = a[ibeg - 1];
                tb = b[ibeg - 1];
                tc = c[ibeg - 1];
                tind = indarr[ibeg - 1];
                a[ibeg - 1] = val;
                upper = ibeg;
                lower = iend;

                while (upper != lower) {
                        while ((upper < lower) && (val < a[lower-1])) {
                                lower--;
                        }

                        if (upper != lower) {
                                a[upper - 1] = a[lower - 1];
                                b[upper - 1] = b[lower - 1];
                                c[upper - 1] = c[lower - 1];
                                indarr[upper - 1] = indarr[lower - 1];
                        }

                        while ((upper < lower) && (val >= a[upper - 1])) {
                                upper++;
                        }

                        if (upper != lower) {
                                a[lower - 1] = a[upper - 1];
                                b[lower - 1] = b[upper - 1];
                                c[lower - 1] = c[upper - 1];
                                indarr[lower - 1] = indarr[upper - 1];
                        }
                }

                *midind = upper;

                /* replace TA in position upper and change midind if ta > val */

                if (ta > val)
                        *midind = upper - 1;


                a[upper - 1] = ta;
                b[upper - 1] = tb;
                c[upper - 1] = tc;
                indarr[upper - 1] = tind;

        } else if (ibeg == iend) {

                if (a[ibeg - 1] < val) 
                        *midind = ibeg;
                else
                        *midind = ibeg - 1;

        } else {

                *midind = ibeg - 1;

        }

        return;

} /* END function cp_partition */



void pc_partition(double *a, double *b, double *c, double *d, int *indarr,
                  int ibeg, int iend, double val, int *midind)
{
    /* local variables */
    double ta, tb, tc, td;
    int lower, upper, tind;
    
    
    if (ibeg < iend) {
        /*
         * temporarily stores ibeg entries and set a(ibeg) = val
         * for the partitioning algorithm.
         */
        ta = a[ibeg - 1];
        tb = b[ibeg - 1];
        tc = c[ibeg - 1];
        td = d[ibeg - 1];
        tind = indarr[ibeg - 1];
        a[ibeg - 1] = val;
        upper = ibeg;
        lower = iend;
        
        while (upper != lower) {
            while ((upper < lower) && (val < a[lower-1])) {
                lower--;
            }
            
            if (upper != lower) {
                a[upper - 1] = a[lower - 1];
                b[upper - 1] = b[lower - 1];
                c[upper - 1] = c[lower - 1];
                d[upper - 1] = d[lower - 1];
                indarr[upper - 1] = indarr[lower - 1];
            }
            
            while ((upper < lower) && (val >= a[upper - 1])) {
                upper++;
            }
            
            if (upper != lower) {
                a[lower - 1] = a[upper - 1];
                b[lower - 1] = b[upper - 1];
                c[lower - 1] = c[upper - 1];
                d[lower - 1] = d[upper - 1];
                indarr[lower - 1] = indarr[upper - 1];
            }
        }
        
        *midind = upper;
        
        
        /* replace TA in position upper and change midind if ta > val */
        
        if (ta > val)
            *midind = upper - 1;
        
        
        a[upper - 1] = ta;
        b[upper - 1] = tb;
        c[upper - 1] = tc;
        d[upper - 1] = td;
        indarr[upper - 1] = tind;
        
    } else if (ibeg == iend) {
        
        if (a[ibeg - 1] < val)
            *midind = ibeg;
        else
            *midind = ibeg - 1;
        
    } else {
        
        *midind = ibeg - 1;
        
    }
    
    return;
    
} /* END function pc_partition */
