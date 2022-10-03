import pandas as pd
import sys

if(len(sys.argv)>1): # check if command line argument for dx was provided
	dx = sys.argv[1]
	print("dx: "+str(dx))

# Read data
curvData = pd.read_csv('curvature.csv', usecols=[0,1,2], names=['time', 'curv', 'curv_ref'])
max_error_curvature = max(abs(curvData['curv']-curvData['curv_ref']))
print("Maximum error curvature: "+str(max_error_curvature))
# Write Max error to disk
f = open("max_error_curvature.csv","w")
g = open("../max_error_curvature.csv","a")
f.write(str(dx)+","+str(max_error_curvature)+'\n')
g.write(str(dx)+","+str(max_error_curvature)+'\n')
f.close()
g.close()


# Read data
angleData = pd.read_csv('contactAngle.csv', usecols=[0,1,2], names=['time', 'angle', 'angle_ref'])
max_error_angle = max(abs(angleData['angle']-angleData['angle_ref']))
print("Maximum error angle: "+str(max_error_angle))
# Write Max error to disk
f = open("max_error_contactAngle.csv","w")
g = open("../max_error_contactAngle.csv","a")
f.write(str(dx)+","+str(max_error_angle)+'\n')
g.write(str(dx)+","+str(max_error_angle)+'\n')
f.close()
g.close()

# Read data
positionData = pd.read_csv('position.csv', usecols=[0,1,2], names=['time', 'position', 'position_ref'])
max_error_position = max(abs(positionData['position']-positionData['position_ref']))
print("Maximum error position: "+str(max_error_position))
# Write Max error to disk
f = open("max_error_position.csv","w")
g = open("../max_error_position.csv","a")
f.write(str(dx)+","+str(max_error_position)+'\n')
g.write(str(dx)+","+str(max_error_position)+'\n')
f.close()
g.close()

# Read data
gradientData = pd.read_csv('gradientNormAtContactPoint.csv', usecols=[0,1], names=['time', 'gradientNorm'])
max_gradient_deviation = max(abs(gradientData['gradientNorm']-1.0))
print("Maximum gradient deviation: "+str(max_gradient_deviation))
# Write Maximum deviation to disk
f = open("max_gradient_deviation.csv","w")
g = open("../max_gradient_deviation.csv","a")
f.write(str(dx)+","+str(max_gradient_deviation)+'\n')
g.write(str(dx)+","+str(max_gradient_deviation)+'\n')
f.close()
g.close()
