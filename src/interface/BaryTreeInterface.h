#ifndef H_BARYTREE_INTERFACE_H
#define H_BARYTREE_INTERFACE_H

    #ifndef H_BARYTREE_TYPES_H
    #define H_BARYTREE_TYPES_H
    
    typedef enum KERNEL
    {
        NO_KERNEL,
        COULOMB,
        YUKAWA,
        REGULARIZED_COULOMB,
        REGULARIZED_YUKAWA,
        ATAN,
        TCF,
        DCF,
        SIN_OVER_R,
        MQ,
        RBS_U,
        RBS_V,
        USER
    } KERNEL;
    
    
    typedef enum SINGULARITY
    {
        NO_SINGULARITY,
        SKIPPING,
        SUBTRACTION
    } SINGULARITY;
    
    
    typedef enum APPROXIMATION
    {
        NO_APPROX,
        LAGRANGE,
        HERMITE
    } APPROXIMATION;
    
    
    typedef enum COMPUTE_TYPE
    {
        NO_COMPUTE_TYPE,
        PARTICLE_CLUSTER,
        CLUSTER_PARTICLE,
        CLUSTER_CLUSTER,
    } COMPUTE_TYPE;
    
    
    #endif /* H_BARYTREE_TYPES_H */

void BaryTreeInterface(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray,
        KERNEL kernel, int numKernelParams, double *kernelParams,
        SINGULARITY singularity, APPROXIMATION approximation, COMPUTE_TYPE compute_type,
		double theta, int interpOrder, int maxPerSourceLeaf, int maxPerTargetLeaf,
        double sizeCheck, double beta, int verbosity);


#endif /* H_BARYTREE_INTERFACE_H */
