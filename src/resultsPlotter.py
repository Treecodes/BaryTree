import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})        # use LaTeX

from matplotlib import rcParams
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern']
rcParams.update({'font.size': 18})




Header = ["Sources", "Targets", "DirectSumComparison","NumSources","NumTargets","Theta","Order","TreeType","MaxParNode", "BatchSize", "Kappa",
          "PotentialType","sflag","pflag","p","time_preproc", "TreeBuildTime", "time2", "time3", "time4", "time5",
          "time6", "time7", "time8", "time9", "time10", "time11", "TreecodeTime", "time13",
          "dpengglob", "tpengglob", "abs_pengerr", "RelativeError", "inferr", "relinferr", "n2err", "reln2err"]


def AversusB(df,A,B,save=False):
    fig, ax = plt.subplots(figsize=(8,6))
    fig.suptitle('%s versus %s' %(A,B))
    df.plot(x=B, y=A, style='o',ax=ax)

    if save == True:
        saveID = A+'Vs'+B
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()

def AversusBcolorbyC(df,A,B,C,save=False):
    fig, ax = plt.subplots(figsize=(8,6))
    fig.suptitle('%s versus %s colored by %s' %(A,B,C))
    grouped = df.groupby(C)
    for name,group in grouped:
##        group.plot(x=B, y=A, style='o', ax=ax, label='%s = %.2f'%(C,name))
        group.plot(x=B, y=A, style='o', ax=ax, label='%s = %s'%(C,name))
    plt.legend(loc = 'best')

    if save == True:
        saveID = A+'Vs'+B+'ColoredBy'+C
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()

def logAversusBcolorbyC(df,A,B,C,save=False):
    fig, ax = plt.subplots(figsize=(8,6))
    fig.suptitle('%s versus %s colored by %s' %(A,B,C))
    grouped = df.groupby(C)
    for name,group in grouped:
        group.plot(x=B, y=A, logy=True, style='o', ax=ax, label='%s = %.2f'%(C,name))
    plt.legend(loc = 'best')

    if save == True:
        saveID = 'log'+A+'Vs'+B+'ColoredBy'+C
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()
    
def logAversusLogBcolorbyC(df,A,B,C,save=False):
    fig, ax = plt.subplots(figsize=(8,6))
#     fig.suptitle('%s versus %s colored by %s' %(A,B,C))
    grouped = df.groupby(C)
    for name,group in grouped:
        if isinstance(name,str):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %s'%(C,name))
        elif isinstance(name,float):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %f'%(C,name))
        elif isinstance(name,int):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %i'%(C,name))
        
    plt.xlabel(B)
    plt.ylabel(A)
#     plt.axhline(y=49.17,color='g', label='Coulomb DS Time')    # GPU Direct Sum for 636000 points
#     plt.axhline(y=61.402522,color='r', label='Yukawa DS Time')    # GPU Direct Sum for 636000 points Yukawa 0.5
    
#     plt.axhline(y=78.54,color='g', label='Coulomb_SS DS Time')    # GPU Direct Sum for 636000 points
#     plt.axhline(y=75.96,color='r', label='Yukawa_SS DS Time')    # GPU Direct Sum for 636000 points
    
#     plt.axhline(y=1600,color='r')     # CPU Direct Sum for 636000 points
#     plt.axhline(y=940,color='r')      # CPU Treecode achieveing 1e-6 relative error
    plt.grid()
#     plt.ylim([20,100])
    
    plt.legend(loc = 'best')

    if save == True:
        saveID = 'log'+A+'VsLog'+B+'ColoredBy'+C
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()
    
def plotScaling(df, df_gpu, save=False):
#     fig, ax = plt.subplots(figsize=(8,6))
    dstimes = np.array([2.101, 25.879, 493.65, 1673.85, 1673.85*(3719492/636608)**2])
    sizes = np.array([21952, 79576, 348488, 636608, 3719492])
    
    quadratic = sizes**2 * (dstimes[0]/sizes[0]**2)*1.3
    nlogn = sizes*np.log(sizes) * (dstimes[0]/( sizes[0] * np.log(sizes[0] ))) *0.7
    nref = sizes * (dstimes[0]/( sizes[0] )) *0.7
    
    df = df.loc[df['NumSources'] != 987125]
    df7 = df.loc[df['Theta']==0.7]
    df7 = df7.loc[df7['Order']==7] 
     
    df9 = df.loc[df['Theta']==0.9]
    df9 = df9.loc[df9['Order']==7]
#     df = df.loc[df['NumSources']==79576]
#     logAversusLogBcolorbyC(df,'TreecodeTime','NumSources', 'Order')
 
    fig, ax = plt.subplots(figsize=(8,6))
#     fig.suptitle('%s versus %s colored by %s' %(A,B,C))
#     df7.plot(x='NumSources', y='TreecodeTime', style='o', ax=ax, loglog=True,label='Treecode Order %i, Theta = %1.1f'%(7,0.7))
#     df9.plot(x='NumSources', y='TreecodeTime', style='o', ax=ax, loglog=True,label='Treecode Order %i, Theta = %1.1f'%(7,0.9))
    df9.plot(x='NumSources', y='TreecodeTime',markerSize=12, style='o', ax=ax, loglog=True,label='Treecode')
#     df_gpu.plot(x='NumSources', y='TreecodeTime',markerSize=12, style='o', ax=ax, loglog=True,label='Treecode - GPU')
    
    plt.plot(sizes,dstimes,'ro',markerSize=12,label='Direct Sum')
    plt.plot(sizes,quadratic,'k--',label=r'$O(N^2)$ Reference')
#     plt.plot(sizes,nlogn,'k-.',label=r'$O(N\log N)$ Reference')
    plt.plot(sizes,nref,'k-.',label=r'$O(N)$ Reference')
    plt.xlabel('Quadrature Points')
    plt.ylabel('Convolution Time (s)')
    plt.legend()
    plt.grid()
    
    if save == True:
        saveID = 'scaling'
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()
 
def gpuSpeedup(df_gpu, df_cpu, theta, order, numPars):
    print('Problem Size: ', numPars)
    df_gpu = df_gpu.loc[df_gpu['Theta']==theta]
    df_gpu = df_gpu.loc[df_gpu['Order']==order]
    df_gpu = df_gpu.loc[df_gpu['NumSources']==numPars]
    
#     print(df_gpu)
    
    df_cpu = df_cpu.loc[df_cpu['Theta']==theta]
    df_cpu = df_cpu.loc[df_cpu['Order']==order]
    df_cpu = df_cpu.loc[df_cpu['NumSources']==numPars]
    
#     print(df_cpu)
    print('CPU time: ', df_cpu['TreecodeTime'])
    print('GPU time: ',df_gpu['TreecodeTime'])
    print('Speedup: ', df_cpu['TreecodeTime']/df_gpu['TreecodeTime'])
    print()
    print()

if __name__=="__main__":
#     resultsDir='/Users/nathanvaughn/Desktop/TreecodeTests/OxygenAtomTests/'
#     resultsFile = 'out636608_noBatchSizeColumn.csv'
#     resultsFile = 'out636608_MaxParNode_32k.csv'
#     resultsFile = 'out636608.csv'
#     resultsFile='out636608_yukawa0p5_versus_Coulomb.csv'
#     resultsFile='out636608_SS_testing.csv'
#     resultsFile='out636608_SS_testing_weighingOutput.csv'




    #### MICDE DATA ####
    resultsFile='tc.csv'
    resultsDir='/Users/nathanvaughn/Documents/synchronizedDataFiles/MICDE_Data_2019/gpu_treecode/'
    df_gpu = pd.read_csv(resultsDir + resultsFile, names=Header)
    resultsFile='tc.csv'
    resultsDir='/Users/nathanvaughn/Documents/synchronizedDataFiles/MICDE_Data_2019/cpu_treecode/'
    df_cpu = pd.read_csv(resultsDir + resultsFile, names=Header)

    
    plotScaling(df_cpu,df_gpu)

    gpuSpeedup(df_gpu, df_cpu, 0.9,5,21952)
    gpuSpeedup(df_gpu, df_cpu, 0.9,5,79576)
    gpuSpeedup(df_gpu, df_cpu, 0.9,5,348488)
    gpuSpeedup(df_gpu, df_cpu, 0.9,5,636608)
    gpuSpeedup(df_gpu, df_cpu, 0.9,5,3719492)
# #     

#     df_cpu = df_cpu.loc[df_cpu['NumSources']==636608]
#     logAversusLogBcolorbyC(df_cpu,'TreecodeTime','RelativeError', 'NumSources')

    