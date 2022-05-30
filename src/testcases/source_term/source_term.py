# Case Definition
# Three functions read_subcases(), read_case_files() and read_subcase_files()
# have to be implemented.

import copy

def build_subcase(source_active, mesh):
	subcase = {}
	subcase["solver_data"] = {}
	subcase["jobscript_data"] = {}
	
	## Dictionaries for additional files
	#subcase["subcase_analysis.py"] = {} 
	#subcase["case_analysis.py"] = {} 
	
	if(source_active==1):
	  subcase["solver_data"]["applySourceTerm"] = "true"
	else:
	  subcase["solver_data"]["applySourceTerm"] = "false"
	
	subcase["solver_data"]["numY"] = mesh
	subcase["solver_data"]["numX"] = mesh*10
	
	return subcase
	

def build_subcase_name(source_active, mesh):
	return str(source_active)+'/'+str(mesh)

	
def read_subcases():

	# Initialize case_data
	case_data = {}
	
	####################################################################
	
	####### Construct subcases via loops ########

	## Arrays to loop over
        sources= [0,1]
        meshes = [10,25,50]
        ##
        
        for source_active in sources:
        	for mesh in meshes:
        		subcase_name = build_subcase_name(source_active,mesh)
        		case_data[subcase_name] = build_subcase(source_active,mesh)
        		
       #######################################################################
      

	return case_data

def read_case_files():
	return [] #["case_analysis.py"]

def read_subcase_files():
	return ["plot.gp"] #["subcase_analysis.py"]
