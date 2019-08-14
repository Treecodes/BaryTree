import numpy as np
import sys

import bitstring
import struct

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def float_to_bin(num):
    return format(struct.unpack('!I', struct.pack('!f', num))[0], '032b')

# convert (x,y) to d
def xy2d( n, x, y):
    rx=0
    ry=0 
    d=0
    s=n/2
    while s>0:
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(n, x, y, rx, ry);
        s/=2
    return d;

def scat(points, orderedPoints):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    
    # For each set of style and range settings, plot n random points in the box
    # defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
#     for c, m, zlow, zhigh in [('r', 'o', -50, -25), ('b', '^', -30, -5)]:
    xs = points[:,0]
    ys = points[:,1]
    zs = points[:,2]
    ax.scatter(xs, ys, zs, c='r', marker='o')
    
    xs = orderedPoints[:,0]
    ys = orderedPoints[:,1]
    zs = orderedPoints[:,2]
    ax.scatter(xs, ys, zs, c='b', marker='o')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    
#     plt.xlim([-1,1])
#     plt.ylim([-1,1])
#     plt.zlim([-1,1])
    
    plt.show()

def mortonOrdering(points):


    bits = np.zeros(len(points),dtype=np.string_ )
    bits = np.empty(len(points),dtype=np.string_ )
    bits=[]
    bitArrayX = []
    bitArrayY = []
    bitArrayZ = []
    for i in range(len(points)):
        
        bitx = bitstring.BitArray(float=points[i,0], length=64).bin
        bity = bitstring.BitArray(float=points[i,1], length=64).bin
        bitz = bitstring.BitArray(float=points[i,2], length=64).bin

        bitArrayX.append(bitx)
        bitArrayY.append(bity)
        bitArrayZ.append(bitz)
        
        idx = ''
        for j in range(len(str(bitArrayX[i]))):
            idx = idx + bitArrayX[i][j]
            idx = idx + bitArrayY[i][j]
            idx = idx + bitArrayZ[i][j]

        bits.append(idx)
    
    index = np.argsort(bits)

        
    orderedPoints = np.zeros_like(points)
    for i in range(len(points)):
        orderedPoints[i,:] = np.copy(points[index[i],:])
    
#     orderedPoints = points 
    return orderedPoints


def measureMaxDistance(points):
    ravg = 0.0
    rmax = 0.0
    rmax_not_boundary = 0.0
    maxidx = 0
    for i in range(len(points)-1):
        dx = points[i,0]-points[i+1,0]
        dy = points[i,1]-points[i+1,1]
        dz = points[i,2]-points[i+1,2]
        
        r = np.sqrt(dx*dx + dy*dy + dz*dz)
        
#         if r>rmax:
#             rmax = r
#             maxidx = i
#             print(i)
#             print(rmax)
#             print()
        rmax = max(r,rmax)
        ravg += r
        
        if ( (i+1)%(len(points)/64) != 0 ):
           rmax_not_boundary = max(rmax_not_boundary,r) 
    
    ravg /= (len(points)-1)
#     print('Max index: ', maxidx)  
    return rmax, ravg, rmax_not_boundary

def initializePoints(pointsPerOcatant):
    
    
    globalPoints = np.empty((0,5))
    for xlow in [-1,-0.5,0,0.5]:
        for ylow in [-1,-0.5,0,0.5]:
            for zlow in [-1,-0.5,0,0.5]:
                
                xhigh = xlow+0.5
                yhigh = ylow+0.5
                zhigh = zlow+0.5
                
                points = 1/2*np.random.rand(pointsPerOcatant,5)
                points[:,0] += xlow
                points[:,1] += ylow
                points[:,2] += zlow
                points[:,3] *= 2
                points[:,3] -= 1
                points[:,-1] = np.ones(pointsPerOcatant) # set quadrature weights to 1
                
                globalPoints = np.concatenate( (globalPoints,points), axis=0)
                
#     print(globalPoints)
    rmax, ravg, rmax_not_boundary = measureMaxDistance(globalPoints)
    print()
    print('original rmax = ', rmax)
    print('original rmax_not_boundary = ', rmax_not_boundary)
    print('orininal ravg = ', ravg)
    print()
    orderedPoints = mortonOrdering(globalPoints)
    rmax, ravg, rmax_not_boundary = measureMaxDistance(orderedPoints)
#     print(orderedPoints)
    print()
    print('final rmax = ', rmax)
    print('final rmax_not_boundary = ', rmax_not_boundary)
    print('final ravg = ', ravg)
    print()
    return globalPoints, orderedPoints
        
        
if __name__=="__main__":
    n=1
    N = int(sys.argv[n]); n+=1
     
#     sourcesTXT = '/Users/nathanvaughn/Desktop/randomPoints/S%ipy.txt' %N
#     targetsTXT = '/Users/nathanvaughn/Desktop/randomPoints/T%ipy.txt' %N
#       
#     points = 2*np.random.rand(N,5)-1.0
#     points[:,-1] = np.ones(N) # set quadrature weights to 1
#        
# #     print(points[0:5,:])
#        
#        
#     orderedPoints = mortonOrdering(points)
#     rmax, ravg, r= measureMaxDistance(points)
#        
#     print('original rmax = ', rmax)
#     print('orininal ravg = ', ravg)
#      
#     rmax, ravg, r = measureMaxDistance(orderedPoints)
#     print('final rmax = ', rmax)
#     print('final ravg = ', ravg)
#      
#     print(points[300:304,0:3])
#     print(orderedPoints[300:304,0:3])
#      
#     scat(points[40000:40000,:],orderedPoints[40000:80000,:])
# #     
# #     np.savetxt(sourcesTXT, points)
# #     np.savetxt(targetsTXT, points[:,0:4])
# 
# #     print( type( float_to_bin(3.14) ) )
# #     f1 = bitstring.BitArray(float=0.57694, length=64)
# #     print(f1.bin)
# #     f=[1.2,3.2,1.4]
# #     int32bits = np.asarray(f, dtype=np.float64).view(np.int64)
# #     int32bits = np.array([1.2,3.2,1.4]).view(np.int32)
# #     print(int32bits)
# #     print('{:032b}'.format(int32bits))
# #     
#     
    
#     print(xy2d(3,1,2))
#     sourcesTXT = '/Users/nathanvaughn/Desktop/randomPoints/S%i.txt' %(64*N)
#     targetsTXT = '/Users/nathanvaughn/Desktop/randomPoints/T%i.txt' %(64*N)

    sourcesTXT = '/scratch/krasny_fluxg/njvaughn/random/S%i.txt' %(64*N)
    targetsTXT = '/scratch/krasny_fluxg/njvaughn/random/T%i.txt' %(64*N)
    points, orderedPoints = initializePoints(N)
  
    np.savetxt(sourcesTXT, orderedPoints)
    np.savetxt(targetsTXT, orderedPoints[:,0:4])
#     orderedPoints.tofile(sourcesTXT)
#     orderedPoints.tofile(targetsTXT)
#     scat(points[15900:16100,:],orderedPoints[15900:16100,:])
#     scat(points[15995:16005,:],orderedPoints[15995:16005,:])
#     scat(points[16000:17000,:],orderedPoints[16000:17000,:])
#     