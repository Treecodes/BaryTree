import pandas as pd
import matplotlib.pyplot as plt

Header = ["Sources", "Targets", "DirectSumComparison","NumSources","NumTargets","Theta","Order","TreeType","MaxParNode", "BatchSize", "Kappa",
          "PotentialType","sflag","pflag","p","time_preproc", "TreeBuildTime", "time2", "time3", "time4", "time5",
          "time6", "time7", "time8", "time9", "time10", "time11", "TreecodeTime", "time13",
          "dpengglob", "tpengglob", "abs_pengerr", "RelativeError", "inferr", "relinferr", "n2err", "reln2err"]

# Header = ["Sources", "Targets", "DirectSumComparison","NumSources","NumTargets","Theta","Order","TreeType","MaxParNode", "Kappa",
#           "PotentialType","sflag","pflag","p","time_preproc", "TreeBuildTime", "time2", "time3", "time4", "time5",
#           "time6", "time7", "time8", "time9", "time10", "time11", "TreecodeTime", "time13",
#           "dpengglob", "tpengglob", "abs_pengerr", "RelativeError", "inferr", "relinferr", "n2err", "reln2err"]

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
    fig.suptitle('Log %s versus %s colored by %s' %(A,B,C))
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
    fig.suptitle('%s versus %s colored by %s' %(A,B,C))
    grouped = df.groupby(C)
    for name,group in grouped:
        if isinstance(name,str):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %s'%(C,name))
        elif isinstance(name,float):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %f'%(C,name))
        elif isinstance(name,int):
            group.plot(x=B, y=A, style='o', ax=ax, loglog=True,label='%s = %i'%(C,name))
        
    plt.legend(loc = 'best')
    plt.xlabel(B)
    plt.ylabel(A)
    plt.axhline(y=49.17,color='r')    # GPU Direct Sum for 636000 points
#     plt.axhline(y=1600,color='r')     # CPU Direct Sum for 636000 points
#     plt.axhline(y=940,color='r')      # CPU Treecode achieveing 1e-6 relative error
    plt.grid()
    
    
    if save == True:
        saveID = 'log'+A+'VsLog'+B+'ColoredBy'+C
        plt.savefig(plotsDir+saveID+'.pdf', bbox_inches='tight',format='pdf')
    plt.show()

if __name__=="__main__":
    resultsDir='/Users/nathanvaughn/Desktop/TreecodeTests/OxygenAtomTests/'
#     resultsFile = 'out636608_noBatchSizeColumn.csv'
    resultsFile = 'out636608_MaxParNode_32k.csv'
    df = pd.read_csv(resultsDir + resultsFile, names=Header)
    print(df)
    
#     df = df.loc[df['Theta'].isin([0.6])]
#     df = df.loc[df['Order'].isin([5])]
#     df = df.loc[df['MaxParNode'].isin([32000])]
#     df = df.loc[df['BatchSize'].isin([4000])]
    
#     logAversusLogBcolorbyC(df,'TreecodeTime','RelativeError', 'Order')
    logAversusLogBcolorbyC(df,'TreecodeTime','RelativeError', 'MaxParNode')
#     logAversusLogBcolorbyC(df,'TreecodeTime','RelativeError', 'Theta')
    logAversusLogBcolorbyC(df,'TreecodeTime','RelativeError', 'BatchSize')
#     logAversusLogBcolorbyC(df,'TreeBuildTime','RelativeError', 'Order')
    