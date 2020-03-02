#ifndef H_OUTPUT_H
#define H_OUTPUT_H

/* declaration of struct with tag particles */
struct output
{
        int numberOfParticles;
        double *potential;
        double *forcesX;
        double *forcesY;
        double *forcesZ;
};

#endif /* H_OUTPUT_H */
