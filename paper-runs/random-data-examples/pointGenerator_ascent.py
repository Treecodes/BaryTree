import numpy as np
import sys

n=1
N = int(sys.argv[n]); n+=1


sourcesTXT ='/gpfs/wolf/proj-shared/gen135/BaryTree/randomPoints/S%i.bin' %N



points = 2*np.random.rand(N,5)-1.0
points[:,-1] = np.ones(N) # set quadrature weights to 1

print(points[0:5,:])

# np.savetxt(sourcesTXT, points)
# np.savetxt(targetsTXT, points[:,0:4])
     
points.tofile(sourcesTXT)   