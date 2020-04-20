import numpy as np
import ctypes
from mpi4py import MPI
import resource
from enum import IntEnum

class CEnum(IntEnum):
    @classmethod
    def from_param(cls, self):
        if not isinstance(self, cls):
            raise TypeError
        return int(self)
    
class Kernel(CEnum):
    NO_KERNEL = 0
    COULOMB = 1
    YUKAWA = 2
    REGULARIZED_COULOMB = 3
    REGULARIZED_YUKAWA = 4
    ATAN = 5
    TCF = 6
    DCF = 7
    SIN_OVER_R = 8
    MQ = 9
    
class Singularity(CEnum):
    NO_SINGULARITY = 0
    SKIPPING = 1
    SUBTRACTION = 2
    
class Approximation(CEnum):
    NO_APPROXIMATION = 0
    LAGRANGE = 1
    HERMITE = 2
    
class ComputeType(CEnum):
    NO_COMPUTE_TYPE = 0
    PARTICLE_CLUSTER = 1
    CLUSTER_PARTICLE = 2
    CLUSTER_CLUSTER = 3
    

""" LOAD TREECODE LIBRARY """
# tries to load .so shared libraries for Linux.  If this fails, tries to load .dylib library for Mac.
# paths may need to be adjusted depending on the install location of the treecode
try: 
    _cpu_treecodeRoutines = ctypes.CDLL('libBaryTree_cpu.so')
except OSError:
        _cpu_treecodeRoutines = ctypes.CDLL('libBaryTree_cpu.dylib')
        
try: 
    _gpu_treecodeRoutines = ctypes.CDLL('libBaryTree_gpu.so')
except OSError:
    try:
        _gpu_treecodeRoutines = ctypes.CDLL('libBaryTree_gpu.dylib')
    except OSError:
        print("Warning: Could not load GPU BaryTree library.  Ignore if not using GPUs.") 
        
    
""" Set argtypes of the wrappers. """
try:
    _gpu_treecodeRoutines.BaryTreeInterface.argtypes = ( ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), Kernel, ctypes.c_int, ctypes.POINTER(ctypes.c_double), Singularity, Approximation, ComputeType,
            ctypes.c_int, ctypes.c_double,  ctypes.c_int,  ctypes.c_int,  ctypes.c_double,  ctypes.c_int )
except NameError:
    print("Warning: Could not set argtypes of _gpu_treecodeRoutines.  Ignore if not using GPUs.")

try:
    _cpu_treecodeRoutines.BaryTreeInterface.argtypes = ( ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), Kernel, ctypes.c_int, ctypes.POINTER(ctypes.c_double),  Singularity, Approximation, ComputeType,
            ctypes.c_int, ctypes.c_double,  ctypes.c_int,  ctypes.c_int,  ctypes.c_double,  ctypes.c_int )
except NameError:
    print("Could not set argtypes of _cpu_treecodeRoutines.")




def callTreedriver(numTargets, numSources, 
                   targetX, targetY, targetZ, targetValue, 
                   sourceX, sourceY, sourceZ, sourceValue, sourceWeight,
                   kernelName, numberOfKernelParameters, kernelParameters, singularityHandling,
                   approximationName, computeType, order, theta, maxParNode, batchSize, GPUpresent, verbosity, sizeCheck=None):
    '''
    python function which creates pointers to the arrays and calls treedriverWrapper.
    returns the results array.
    '''

    c_double_p = ctypes.POINTER(ctypes.c_double)
    
    targetX_p = targetX.ctypes.data_as(c_double_p)
    targetY_p = targetY.ctypes.data_as(c_double_p)
    targetZ_p = targetZ.ctypes.data_as(c_double_p)
    targetValue_p = targetValue.ctypes.data_as(c_double_p)
    
    sourceX_p = sourceX.ctypes.data_as(c_double_p)
    sourceY_p = sourceY.ctypes.data_as(c_double_p)
    sourceZ_p = sourceZ.ctypes.data_as(c_double_p)  
    sourceValue_p =  sourceValue.ctypes.data_as(c_double_p)
    sourceWeight_p = sourceWeight.ctypes.data_as(c_double_p)

    kernelParameters_p = kernelParameters.ctypes.data_as(c_double_p)
    
    resultArray = np.zeros(numTargets)
    resultArray_p = resultArray.ctypes.data_as(c_double_p)
    
    if not sizeCheck:
        if approximationName == Approximation.LAGRANGE:
            sizeCheck = 1.0
        if approximationName == Approximation.HERMITE:
            sizeCheck = 4.0
    
    if GPUpresent==True:
        _gpu_treecodeRoutines.BaryTreeInterface(ctypes.c_int(numTargets),  ctypes.c_int(numSources),
                                                targetX_p, targetY_p, targetZ_p, targetValue_p,
                                                sourceX_p, sourceY_p, sourceZ_p, sourceValue_p, sourceWeight_p,
                                                resultArray_p, kernelName, ctypes.c_int(numberOfKernelParameters), kernelParameters_p,
                                                singularityHandling, approximationName, computeType,
                                                ctypes.c_int(order), ctypes.c_double(theta), ctypes.c_int(maxParNode), ctypes.c_int(batchSize), ctypes.c_double(sizeCheck), ctypes.c_int(verbosity) )
    elif GPUpresent==False: # No gpu present
        _cpu_treecodeRoutines.BaryTreeInterface(ctypes.c_int(numTargets),  ctypes.c_int(numSources),
                                                targetX_p, targetY_p, targetZ_p, targetValue_p,
                                                sourceX_p, sourceY_p, sourceZ_p, sourceValue_p, sourceWeight_p,
                                                resultArray_p, kernelName, ctypes.c_int(numberOfKernelParameters), kernelParameters_p,
                                                singularityHandling, approximationName, computeType,
                                                ctypes.c_int(order), ctypes.c_double(theta), ctypes.c_int(maxParNode), ctypes.c_int(batchSize), ctypes.c_double(sizeCheck), ctypes.c_int(verbosity) )
    else: 
        print("What should GPUpresent be set to in the wrapper?")
        exit(-1) 
    
    
    return resultArray
