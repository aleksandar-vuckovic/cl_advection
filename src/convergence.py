#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import shutil
import sys
import os
import subprocess

DeltaXArray = np.linspace(0.001, 0.00001, 100)
lenX = 1.0
time = 1.0

try:
    if (sys.argv[1] == "True"):
        submitJobs = True
    else:
        submitJobs = False
except Exception:
    submitJobs = False
    

def DeltaT(DeltaX):
    vNorm = 1
    CFL = 0.2
    return CFL*DeltaX/vNorm

for dx in DeltaXArray:
    directory = "generatedTestcases/" + str(dx)
    
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    InputfileTemplate = open("templates/InputfileTemplate", 'r')
    shutil.copyfile("templates/InputfileTemplate", directory+"/Inputfile")
    shutil.copyfile("templates/jobscriptTemplate", directory+"/jobscript")
    outputfile = open(directory+"/Inputfile", 'w')
    for line in InputfileTemplate:
        line = line.replace("$numX", str(int(lenX/dx)))
        line = line.replace("$numY", str(int(lenX/dx)))
        line = line.replace("$lenX", str(lenX))
        line = line.replace("$lenY", str(lenX))
        line = line.replace("$lenZ", str(dx))
        line = line.replace("$timesteps", str(int(time/DeltaT(dx))))
        line = line.replace("$time", str(time))
        line = line.replace("$writesteps", str(int(time/DeltaT(dx))//10))    
        outputfile.write(line)
        
    if (submitJobs):
        subprocess.run(["sbatch", "jobscript"],  cwd = directory)
        
