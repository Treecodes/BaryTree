'''
'''
import os
import sys
import resource
import numpy as np
import mpi4py.MPI as MPI



sys.path.append(os.getcwd())
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
    
    # initialize some random data
    RHO = np.random.rand(N)
    X = np.random.rand(N)
    Y = np.random.rand(N)
    Z = np.random.rand(N)
    W = np.ones(N)
    
    # measure memory usage
    initialMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    previousMemory=initialMemory
    print('Initial memory usage: ', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss )
    
    # call the treecode n times, measuring memory each time
    for i in range(n):
        
        output = treecodeWrappers.callTreedriver(  N, N, 
                                                   X, Y, Z, RHO, 
                                                   X, Y, Z, RHO, W,
                                                   kernelName, numberOfKernelParameters, kernelParameters, singularityHandling, approximationName,
                                                   treecodeOrder, theta, maxParNode, batchSize, GPUpresent, verbosity)
    
        
        newMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        print("Memory growth this iteration: %i to %i " %(previousMemory,newMemory))
        previousMemory=newMemory
        
        
    finalMemory = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    print("Memory growth over %i iterations: %i to %i, difference = %i" %(n,initialMemory,finalMemory,finalMemory-initialMemory))
    
    
    print("\n\nIf no errors occured, and memory didn't blow up, then the wrapper is working as expected!\n\n")
        



