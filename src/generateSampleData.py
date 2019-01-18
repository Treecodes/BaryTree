import numpy as np
import os


# Generate Data
numSources = 100000
numTargets = 100000

sourcesTXT = '../examples/S%ipy.txt' %numSources
targetsTXT = '../examples/T%ipy.txt' %numTargets

sourcesTXT = '/Users/nathanvaughn/Documents/testData/H2Sources.txt'
targetsTXT = '/Users/nathanvaughn/Documents/testData/H2Targets.txt'


Sources = np.random.rand(numSources,5)
Targets = np.random.rand(numTargets,4)


# Save as .txt files
np.savetxt(sourcesTXT, Sources)
np.savetxt(targetsTXT, Targets)


# # Call the txt2bin executable
# os.system('../bin/txt2bin.exe %s' %sourcesTXT)
# os.system('../bin/txt2bin.exe %s' %targetsTXT)
# print('Done converting to ..bin')

print('Done.')