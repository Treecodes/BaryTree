import numpy as np
import sys

if (sys.argv[1]).lower() == 'sources':
    outfile = sys.argv[2]
    N = int(sys.argv[3])

    sources = 2*np.random.rand(N,5)-1.0
    sources[:,-1] = np.ones(N) # set quadrature weights to 1

    sources.tofile(outfile)

elif (sys.argv[1]).lower() == 'targets':
    outfile = sys.argv[2]
    N = int(sys.argv[3])

    targets = 2*np.random.rand(N,4)-1.0
    targets[:,-1] = np.ones(N) # set target weights to 1

    targets.tofile(outfile)
else:
    print("Error! Neither sources nor target output specified. Exiting.")
