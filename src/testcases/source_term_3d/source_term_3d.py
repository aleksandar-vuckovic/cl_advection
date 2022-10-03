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
	
	if(source_active==1):
	  subcase["solver_data"]["applySourceTerm"] = "true"
	else:
	  subcase["solver_data"]["applySourceTerm"] = "false"
	  
 
	## define the mesh
	if(mesh==50):
		subcase["solver_data"]["numZ"] = 50
		subcase["solver_data"]["numX"] = 50
		subcase["solver_data"]["numY"] = 15
	elif(mesh==100):
		subcase["solver_data"]["numZ"] = 100
		subcase["solver_data"]["numX"] = 100
		subcase["solver_data"]["numY"] = 30
	elif(mesh==200):
		subcase["solver_data"]["numZ"] = 200
		subcase["solver_data"]["numX"] = 200
		subcase["solver_data"]["numY"] = 60

	## define the field
	if(setup=="strawberry"):
	
		subcase["solver_data"]["time"]=2.0
	
		subcase["solver_data"]["field"] = """
v0=0.3
w0=0.4
c1=0.1
c2=0.1
c3=-0.2
c4=0.3
c5=-0.1
c6=0.1
field=strawberryField"""

	return subcase	
	

def build_subcase_name(source_active, mesh,setup,cfl):
	return 'cfl_'+str(cfl)+'/'+setup+"/"+str(source_active)+'/'+str(mesh)

	
def read_subcases():

	# Initialize case_data
	case_data = {}
	
	####################################################################
	
	####### Construct subcases via loops ########

	## Arrays to loop over
	setups=["strawberry"]
        sources= [0,1]
        meshes = [50,100,200]
        cfls = [0.2]
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
	return ["plot.gp"] #["subcase_analysis.py"]
