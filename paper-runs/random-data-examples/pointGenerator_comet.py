


import numpy as np
import sys

n=1
N = int(sys.argv[n]); n+=1


<<<<<<< HEAD
sourcesTXT = '/oasis/scratch/comet/njvaughn/temp_project/random/S%ipy.txt' %N
targetsTXT = '/oasis/scratch/comet/njvaughn/temp_project/random/T%ipy.txt' %N
=======
# sourcesTXT = '/scratch/krasny_fluxg/njvaughn/random/S%ipy.txt' %N
# targetsTXT = '/scratch/krasny_fluxg/njvaughn/random/T%ipy.txt' %N

sourcesTXT = '/Users/nathanvaughn/Desktop/randomPoints/S%i.txt' %N
targetsTXT = '/Users/nathanvaughn/Desktop/randomPoints/T%i.txt' %N
>>>>>>> 85c431e862710acb5d23b4d83faa3205134dab38



points = 2*np.random.rand(N,5)-1.0
points[:,-1] = np.ones(N) # set quadrature weights to 1

print(points[0:5,:])

np.savetxt(sourcesTXT, points)
np.savetxt(targetsTXT, points[:,0:4])
   
