'''
'''
import os
import sys
import resource
import numpy as np
import mpi4py.MPI as MPI


sys.path.append(os.getcwd())
try:
    import BaryTreeInterface as BT
except ImportError:
    print('Unable to import BaryTreeInterface due to ImportError')
except OSError:
    print('Unable to import BaryTreeInterface due to OSError')


if __name__=="__main__":
    
    # set treecode parameters
    N = 5000
    maxPerSourceLeaf = 50
    maxPerTargetLeaf = 10
    GPUpresent = False
    theta = 0.8
    treecodeDegree = 4
    gaussianAlpha = 1.0
    verbosity = 0
    
    approximation = BT.Approximation.LAGRANGE
    singularity   = BT.Singularity.SUBTRACTION
    computeType   = BT.ComputeType.PARTICLE_CLUSTER
    
    kernel = BT.Kernel.YUKAWA
    numberOfKernelParameters = 1
    kernelParameters = np.array([0.5])


    # initialize some random data
    np.random.seed(1)
    RHO = np.random.rand(N)
    X = np.random.rand(N)
    Y = np.random.rand(N)
    Z = np.random.rand(N)
    W = np.ones(N)   # W stores quadrature weights for convolution integrals.  For particle simulations, set = ones.
    
    expectedOutput = 588.7432483318685  # using seed of 1, this is expected value of first element of output array.
    

    # call the treecode
        
    output = BT.callTreedriver(  N, N,
                                 X, Y, Z, RHO,
                                 np.copy(X), np.copy(Y), np.copy(Z), np.copy(RHO), np.copy(W),
                                 kernel, numberOfKernelParameters, kernelParameters,
                                 singularity, approximation, computeType,
                                 GPUpresent, verbosity, 
                                 theta=theta, degree=treecodeDegree, sourceLeafSize=maxPerSourceLeaf, targetLeafSize=maxPerTargetLeaf, sizeCheck=1.0)

    assert (abs(output[0]-expectedOutput) < 1e-14), "Error: didn't get the expected output using explicit theta/degree."
    
    
    
    
    
    
    beta = 0.1
    expectedOutput = 588.7445889051367  # this is expected value of first element of output array for beta = 0.1
    output = BT.callTreedriver(  N, N,
                                 X, Y, Z, RHO,
                                 np.copy(X), np.copy(Y), np.copy(Z), np.copy(RHO), np.copy(W),
                                 kernel, numberOfKernelParameters, kernelParameters,
                                 singularity, approximation, computeType,
                                 GPUpresent, verbosity, beta=beta,  sizeCheck=1.0)
    assert (abs(output[0]-expectedOutput) < 1e-14), "Error: didn't get the expected output using beta."
    
    
    print("If no errors printed, then the calls to the treecode wrapper worked (one using explicit theta/degree, one use beta)")



