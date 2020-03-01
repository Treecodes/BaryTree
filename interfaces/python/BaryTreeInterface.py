import numpy as np
import ctypes
from mpi4py import MPI
import resource



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
        print("Could not load GPU treecode library.") 
        
    
""" Set argtypes of the wrappers. """
try:
    _gpu_treecodeRoutines.BaryTreeInterface.argtypes = ( ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.POINTER(ctypes.c_double),  ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char),
            ctypes.c_int, ctypes.c_double,  ctypes.c_int,  ctypes.c_int,  ctypes.c_int ) 
except NameError:
    print("Could not set argtypes of _gpu_treecodeRoutines")

try:
    _cpu_treecodeRoutines.BaryTreeInterface.argtypes = ( ctypes.c_int, ctypes.c_int,
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_char), ctypes.c_int, ctypes.POINTER(ctypes.c_double),  ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char),
            ctypes.c_int, ctypes.c_double,  ctypes.c_int,  ctypes.c_int,  ctypes.c_int ) 
except NameError:
    print("Could not set argtypes of _cpu_treecodeRoutines")

print('_treecodeRoutines set.')




def callTreedriver(numTargets, numSources, 
                   targetX, targetY, targetZ, targetValue, 
                   sourceX, sourceY, sourceZ, sourceValue, sourceWeight,
                   kernelName, numberOfKernelParameters, kernelParameters, singularityHandling, approximationName, order, theta, maxParNode, batchSize, GPUversion, verbosity):
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
    
    b_kernelName = kernelName.encode('utf-8')
    b_approximationName = approximationName.encode('utf-8')
    b_singularityHandling = singularityHandling.encode('utf-8')
    
    if GPUversion==True:
        _gpu_treecodeRoutines.BaryTreeInterface(ctypes.c_int(numTargets),  ctypes.c_int(numSources),
                                                     targetX_p, targetY_p, targetZ_p, targetValue_p,
                                                     sourceX_p, sourceY_p, sourceZ_p, sourceValue_p, sourceWeight_p,
                                                     resultArray_p, b_kernelName, ctypes.c_int(numberOfKernelParameters), kernelParameters_p,
                                                     b_singularityHandling, b_approximationName,
                                                     ctypes.c_int(order), ctypes.c_double(theta), ctypes.c_int(maxParNode), ctypes.c_int(batchSize), ctypes.c_int(verbosity) )
    elif GPUversion==False: # No gpu present
        _cpu_treecodeRoutines.BaryTreeInterface(ctypes.c_int(numTargets),  ctypes.c_int(numSources),
                                                     targetX_p, targetY_p, targetZ_p, targetValue_p,
                                                     sourceX_p, sourceY_p, sourceZ_p, sourceValue_p, sourceWeight_p,
                                                     resultArray_p, b_kernelName, ctypes.c_int(numberOfKernelParameters), kernelParameters_p,
                                                     b_singularityHandling, b_approximationName,
                                                     ctypes.c_int(order), ctypes.c_double(theta), ctypes.c_int(maxParNode), ctypes.c_int(batchSize), ctypes.c_int(verbosity) ) 
    else: 
        print("What should GPUversion be set to in the wrapper?")
        exit(-1) 
    
    
    
    return resultArray
