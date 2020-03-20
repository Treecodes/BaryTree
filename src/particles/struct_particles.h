#ifndef H_PARTICLES_H
#define H_PARTICLES_H

/* declaration of struct with tag particles */
struct Particles
{
        int num;
        double *x;
        double *y;
        double *z;
        double *q;
        // quadrature weights.  Set = 1 if interacting particles, not performing convolution integral.
        double *w;

		int *ibeg;
		int *iend;
  
        int *order;

};

#endif /* H_PARTICLES_H */
