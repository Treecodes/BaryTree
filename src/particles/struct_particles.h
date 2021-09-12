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

		int *ibeg;
		int *iend;
  
        int *order;

        double xmin;
        double xmax;
        
        double ymin;
        double ymax;
        
        double zmin;
        double zmax;

        int xdim;
        int ydim;
        int zdim;

        double xdd;
        double ydd;
        double zdd;
};

#endif /* H_PARTICLES_H */
