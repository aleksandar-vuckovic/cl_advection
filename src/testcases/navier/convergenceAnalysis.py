    # -*- coding: utf-8 -*-

import sys
import os
import numpy as np
from matplotlib import pyplot as plt

#testcases are in the same folder as the script
folder = "./"

DeltaX = []
maxDifferenceAngle = []

theoreticalPlotted = True

#def reference(t):
 #   c1 = 0.1
  #  c2 = -2
   # return np.pi/2 + np.arctan(-1/np.tan(theta0) * np.exp(2*c1*t) + c2 * (np.exp(2*c1*t) -1)/(2*c1))

#def reference(t):
 #   return np.pi/2 + np.arctan(-1/np.sqrt(3) * np.exp(0.2*t) + 10*(np.exp(0.2*t) - 1))
for case in os.listdir(folder):
    if (os.path.isfile(folder + case)):
        continue
    angleFile = open(folder + case + "/contactAngle.csv")
    angleData = np.genfromtxt(folder + case + "/contactAngle.csv", delimiter = ",")
    
    
    #if (not theoreticalPlotted):
     #   plt.plot(curvatureData[:, 0], curvatureData[:, 2])
      #  plt.xlabel("Time in arbitrary units")
        #plt.ylabel("Curvature")
       # theoreticalPlotted = True
    


    DeltaX.append(1/float(case[7:]))
    
    timesteps = []
    theoreticalAngle = []
    actualAngle = []
    
    try:
        length = len(angleData[:, 0])
    except IndexError:
        length = 1
        
    for i in range(length):
        time = angleData[i, 0]
        timesteps.append(time)
	theoreticalAngle.append(angleData[i, 2])
	actualAngle.append(angleData[i, 1])

    theoreticalAngle = np.array(theoreticalAngle)
    actualAngle = np.array(actualAngle)
    
    maxDiffAngle = max(np.abs(theoreticalAngle - actualAngle))
    
    maxDifferenceAngle.append(maxDiffAngle)

print(maxDifferenceAngle)
print(DeltaX)
plt.plot(DeltaX, maxDifferenceAngle, linestyle="None", marker ="x")
plt.xlabel("Cell Width")
plt.ylabel("max difference in angle")
#plt.savefig("maxDiffAngle.pdf")
plt.show()     
