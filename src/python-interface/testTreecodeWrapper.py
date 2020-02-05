'''
'''
import os
import sys
import resource
import numpy as np
import mpi4py.MPI as MPI


sys.path.insert(1, '/Users/nathanvaughn/Documents/GitHub/BaryTree/src/python-interface')

try:
    import treecodeWrappers
except ImportError:
    print('Unable to import treecodeWrapper due to ImportError')
except OSError:
    print('Unable to import treecodeWrapper due to OSError')
    import treecodeWrappers_distributed as treecodeWrappers


if __name__=="__main__":
    
    # set treecode parameters
    maxParNode=500
    batchSize=500
    GPUpresent=False
    theta=0.8
    treecodeOrder=7
    gaussianAlpha=1.0
    approximationName = "lagrange"
    singularityHandling = "subtraction"
    verbosity=0
    N=10000
    
    kernelName = "regularized-yukawa"
    numberOfKernelParameters=2
    kernelParameters=np.array([0.5, 0.1])

    
    # set number of iterations for the treecode wrapper calls
    n=10
    print('Initial memory usage before creating arrays: ', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss )
    # initialize some random data
    RHO = np.random.rand(N)
    X = np.random.rand(N)
    Y = np.random.rand(N)
    Z = np.random.rand(N)
    W = np.ones(N)
    
    
    initialMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    previousMemory=initialMemory
    print('Initial memory usage: ', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss )
    
    for i in range(10):
        
        output = treecodeWrappers.callTreedriver(  N, N, 
                                                   X, Y, Z, RHO, 
                                                   X, Y, Z, RHO, W,
                                                   kernelName, numberOfKernelParameters, kernelParameters, singularityHandling, approximationName,
                                                   treecodeOrder, theta, maxParNode, batchSize, GPUpresent, verbosity)
    
        newMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print("Memory growth this iteration: %i to %i " %(previousMemory,newMemory))
        previousMemory=newMemory
        
    finalMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("\nMemory growth over %i iterations: %i to %i" %(n,initialMemory,finalMemory))
        



