# Case Definition
# Three functions read_subcases(), read_case_files() and read_subcase_files()
# have to be implemented.

import copy
import math

def build_subcase(source_active, mesh, setup, cfl):
	subcase = {}
	subcase["solver_data"] = {}
	subcase["jobscript_data"] = {}
	
	## Dictionaries for additional files
	#subcase["subcase_analysis.py"] = {} 
	#subcase["case_analysis.py"] = {} 
	
	subcase["solver_data"]["cfl"] = cfl
	subcase["jobscript_data"]["dx"] = 0.5/mesh
	
	if(source_active==1):
	  subcase["solver_data"]["applySourceTerm"] = "true"
	else:
	  subcase["solver_data"]["applySourceTerm"] = "false"
	  
	  
	geometry1 =  """
centerX=0.5
centerY=-0.15
centerZ=0.0
radius=0.3

expcpX=0.7598076211353315
expcpY=0.0
expcpZ=0.0
expAngle=60.0""" 

	geometry2 =  """
centerX=0.5
centerY=0.15
centerZ=0.0
radius=0.3

expcpX=0.7598076211353315
expcpY=0.0
expcpZ=0.0
expAngle=120.0""" 
	
	subcase["solver_data"]["numY"] = mesh
	subcase["solver_data"]["numX"] = mesh*2
	subcase["solver_data"]["lenZ"] = 1.0/mesh
	

	## define the field
	if(setup=="periodic"):
	
		subcase["solver_data"]["time"]=1.0
	 	subcase["solver_data"]["geometry"] = geometry1
	
		subcase["solver_data"]["field"] = """
v0=-0.2
c1=0.1
c2=-2
tau=0.4
field=timeDependentNavierField"""

	elif(setup=="linear"):
		subcase["solver_data"]["time"]=1.0
	 	subcase["solver_data"]["geometry"] = geometry1
		subcase["solver_data"]["field"] = """
v0=+0.8
c1=-1.0
c2=0
field=navierField"""

	elif(setup=="shear"):
		subcase["solver_data"]["time"]=1.0
	 	subcase["solver_data"]["geometry"] = geometry1
		subcase["solver_data"]["field"] = """
v0=-0.2
field=shearField"""

	
	elif(setup=="largeAngle"):
		subcase["solver_data"]["time"]=1.0
	 	subcase["solver_data"]["geometry"] = geometry2
		subcase["solver_data"]["field"] = """
v0=-0.2
field=shearField"""

	
	return subcase	
	

def build_subcase_name(source_active, mesh,setup,cfl):
	return 'cfl_'+str(cfl)+'/'+setup+"/"+str(source_active)+'/'+str(mesh)

	
def read_subcases():

	# Initialize case_data
	case_data = {}
	
	####################################################################
	
	####### Construct subcases via loops ########

	## Arrays to loop over
	setups=["periodic","shear"] #,"largeAngle","linear"]
        sources= [0,1]
        meshes = [100,200,400]
        cfls = [0.2,0.5]
        ##
        
        for source_active in sources:
        	for mesh in meshes:
        		for setup in setups:
        			for cfl in cfls:
        				subcase_name = build_subcase_name(source_active,mesh,setup,cfl)
	        			case_data[subcase_name] = build_subcase(source_active,mesh,setup,cfl)
        		
       #######################################################################
      

	return case_data

def read_case_files():
	return [] #["case_analysis.py"]

def read_subcase_files():
	return ["plot.gp","readData.py"]
