'''
'''
import os
import sys
import resource
import numpy as np
import mpi4py.MPI as MPI
import time


sys.path.append(os.getcwd())
try:
    import BaryTreeInterface
except ImportError as e:
    print('Unable to import BaryTreeInterface due to ImportError')
    print(e)
#     exit(-1)
except OSError as e:
    print('Unable to import BaryTreeInterface due to OSError')
    print(e)
#     exit(-1)


if __name__=="__main__":
    
    # set treecode parameters 
    maxParNode=2000
    batchSize=2000 
    GPUpresent=True 
    theta=0.8
    treecodeOrder=4
    gaussianAlpha=1.0
    approximationName = "lagrange"
    singularityHandling = "subtraction"
    verbosity=0
    N=500000
    
    kernelName = "yukawa"
    numberOfKernelParameters=1
    kernelParameters=np.array([0.5])


    # initialize some random data
    np.random.seed(1)
    RHO = np.random.rand(N)
    X = np.random.rand(N)
    Y = np.random.rand(N)
    Z = np.random.rand(N)
    W = np.ones(N)   # W stores quadrature weights for convolution integrals.  For particle simulations, simply set = ones.
    
    expectedOutput=-977.407950538299  # using seed of 1, this is the expected value of the first element of the output array.
    

    # call the treecode
        
    start=time.time()
    output = BaryTreeInterface.callTreedriver(  N, N, 
                                                X, Y, Z, RHO, 
                                                X, Y, Z, RHO, W,
                                                kernelName, numberOfKernelParameters, kernelParameters,
                                                singularityHandling, approximationName,
                                                treecodeOrder, theta, maxParNode, batchSize, GPUpresent, verbosity)
    
    end=time.time()
    
#     assert (abs(output[0]-expectedOutput)<1e-14), "Error: didn't get the expected output."
    print("If no errors printed, then the call to the treecode wrapper worked!\nIt took %f seconds." %(end-start))



