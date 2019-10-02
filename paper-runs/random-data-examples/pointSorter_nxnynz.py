import numpy as np
import sys
import time





def sortIntoBoxes(points,nx,ny,nz):
    start = time.time()
    
    Lx = np.max(points[:,0]) - np.min(points[:,0])
    Ly = np.max(points[:,1]) - np.min(points[:,1])
    Lz = np.max(points[:,2]) - np.min(points[:,2])
    
    
    offsets = np.zeros(nx*ny*nz,dtype=np.int)
    sortedPoints = np.zeros_like(points)
    
    
    procID = 0
    for i in range(nx):
        
        xlow = np.min(points[:,0]) + i*Lx/nx
        xhigh = np.min(points[:,0]) + (i+1)*Lx/nx
        
        for j in range(ny):
            
            ylow = np.min(points[:,1]) + j*Ly/ny
            yhigh = np.min(points[:,1]) + (j+1)*Ly/ny
            
            for k in range(nz):
                
                zlow = np.min(points[:,2]) + k*Lz/nz
                zhigh = np.min(points[:,2]) + (k+1)*Lz/nz
                
                temp = np.copy(points)
                
                if i==0:  # Usually points at the boundary get accepted into the high side.  Except on leftmost processor (i=0), in which case boundary points must be included
                    temp = temp[temp[:,0]>=xlow]
                else:
                    temp = temp[temp[:,0]>xlow]
                if j==0:
                    temp = temp[temp[:,1]>=ylow]
                else:
                    temp = temp[temp[:,1]>ylow]
                if k==0:
                    temp = temp[temp[:,2]>=zlow]
                else:
                    temp = temp[temp[:,2]>zlow]

                temp = temp[temp[:,0]<=xhigh]
                temp = temp[temp[:,1]<=yhigh]
                temp = temp[temp[:,2]<=zhigh]
                
                
                
#                 print(temp)
                
#                 print(np.shape(temp))
                
                ## Bounding box is [xlow,xhigh] x [ylow,yhigh] x [zlow,zhigh]
                ## Collect all points that lie in this region
#                 print('ijk',i,j,k)
#                 print("Bounding Box: [%1.2f, %1.2f] x [%1.2f, %1.2f] x [%1.2f, %1.2f]" %(xlow,xhigh,ylow,yhigh,zlow,zhigh))
#                 
                if len(temp)>0:
                    assert np.min(temp[:,0]>=xlow),   'xlow bound not respected'
                    assert np.max(temp[:,0]<=xhigh), 'xhigh bound not respected'
                    assert np.min(temp[:,1]>=ylow),   'ylow bound not respected'
                    assert np.max(temp[:,1]<=yhigh), 'yhigh bound not respected'
                    assert np.min(temp[:,2]>=zlow),   'zlow bound not respected'
                    assert np.max(temp[:,2]<=zhigh),   'zhigh bound not respected'
#                     print(temp)
#                 print(len(temp))
#                 print()
#                 if ( (procID > 0) and (procID+1)<nx*ny*nz):
                if procID+1<nx*ny*nz:
#                 if procID > 0:
#                     print('procID=%i, setting offests[%i] = %i+%i' %(procID,procID+1,offsets[procID],len(temp)))
                    offsets[procID+1] = offsets[procID]+len(temp)
#                 elif procID==0:
#                     offsets[procID]=0
                    
#                 sortedPoints[offsets[procID]:offsets[procID]+len(temp),:] = np.copy(temp)
                offset = offsets[procID]
#                 print(offset)
                np.copyto(sortedPoints[offset:offset+len(temp),:],temp)
#                 print(sortedPoints)
#                 print()
                
                procID+=1
                
#     print(offsets)
#     print(sortedPoints)
    
#     for i in range(len(offsets)):
#         print('Processor %i gets: ' %i)
#         offset=offsets[i]
#         if i<(len(offsets)-1):
#             end=offsets[i+1]
#         else:
#             end=len(sortedPoints)
#         
#         print(sortedPoints[offset:end,:])
#         print()
    
#     sortedPoints.tofile(sourcesTXT)
#     offsets.tofile(targetsTXT)

    print(sortedPoints[:8,:])
    print(sortedPoints[-8:,:])
    assert abs(np.sum(sortedPoints)-np.sum(points))/np.sum(points)<1e-12, 'Checksum failed..., sum sorted = %f, sum original = %f'%(np.sum(sortedPoints),np.sum(points))
                
    end=time.time()
    print("Time to sort %i particles into %i boxes: %f" %(npoints,nx*ny*nz,end-start))
    
    return sortedPoints, offsets
                
                
    
        
        
if __name__=="__main__":
#     npoints=3*10**2
#     points = np.random.rand(npoints,3)
#     print(np.shape(points))
#     sortIntoBoxes(points,2,2,2)
    
    argn=1
    npoints     = int(sys.argv[argn]); argn+=1
    nx          = int(sys.argv[argn]); argn+=1
    ny          = int(sys.argv[argn]); argn+=1
    nz          = int(sys.argv[argn]); argn+=1
    directory   = str(sys.argv[argn]); argn+=1
    
    
    sourcesBin = directory + 'S%i.bin' %(npoints)
    points = np.fromfile(sourcesBin)# read in sources
    points = np.reshape(points, (int(len(points)/5),5))
    
    sortedPoints,offsets = sortIntoBoxes(points,nx,ny,nz)
    
    sortedPoints.tofile(directory+'S%i_%ix_%iy_%iz.bin' %(npoints,nx,ny,nz))
    sortedPoints[:,0:4].tofile(directory+'T%i_%ix_%iy_%iz.bin' %(npoints,nx,ny,nz))
    offsets.tofile(directory+'offsets%i_%ix_%iy_%iz.bin' %(npoints,nx,ny,nz))

    
    
    