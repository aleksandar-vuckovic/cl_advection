#!/bin/bash

toplevel=$(readlink -f ../../)
echo $toplevel
cd $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_200
$toplevel/src/cl_advection -t 4 -w false -o $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_200
$toplevel/src/cl_advection -t 4 -w false -r 0.5 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_100
$toplevel/src/cl_advection -t 4 -w false -r 0.25 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_50
$toplevel/src/cl_advection -t 4 -w false -r 0.125 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_25
cd $toplevel/src/testcaseScripts; ./shear_2D;
if [ $? -eq 0 ]; then
    shear_2D_success=true
else
    shear_2D_success=false
fi

cd $toplevel/src/testcases/2D/curvature/convergenceAnalysis/shear/shear_200											
$toplevel/src/cl_advection -t 4 -w false -f navierField -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/navier/navier_200					
$toplevel/src/cl_advection -t 4 -w false -f navierField -r 0.5 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/navier/navier_100				
$toplevel/src/cl_advection -t 4 -w false -f navierField -r 0.25 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/navier/navier_50			      
$toplevel/src/cl_advection -t 4 -w false -f navierField -r 0.125 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/navier/navier_25
cd $toplevel/src/testcaseScripts; ./navier_2D;
if [ $? -eq 0 ]; then
    navier_2D_success=true
else
    navier_2D_success=false
fi

#cd $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/quadratic_100										
#$toplevel/src/cl_advection -t 4 -w false -f quadraticField -r 2 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/quadratic/quadratic_200	
#$toplevel/src/cl_advection -t 4 -w false -f quadraticField -r 1 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/quadratic/quadratic_100	
#$toplevel/src/cl_advection -t 4 -w false -f quadraticField -r 0.5 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/quadratic/quadratic_50		
#$toplevel/src/cl_advection -t 4 -w false -f quadraticField -r 0.25 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/quadratic/quadratic_25	

#cd $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/navier_100										
#$toplevel/src/cl_advection -t 4 -w false -f navierField -r 2 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/navier/navier_200		
#$toplevel/src/cl_advection -t 4 -w false -f navierField -r 1 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/navier/navier_100		
#$toplevel/src/cl_advection -t 4 -w false -f navierField -r 0.5 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/navier/navier_50		       
#$toplevel/src/cl_advection -t 4 -w false -f navierField -r 0.25 -o  $toplevel/src/testcases/2D/curvature/convergenceAnalysis/paper_examples/navier/navier_25



cd $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_200
$toplevel/src/cl_advection -t 4 -w false -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_200
$toplevel/src/cl_advection -t 4 -w false -r 0.5 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_100
$toplevel/src/cl_advection -t 4 -w true -r 0.25 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_50
$toplevel/src/cl_advection -t 4 -w false -r 0.125 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_25
cd $toplevel/src/testcaseScripts; ./shear_3D;
if [ $? -eq 0 ]; then
    shear_3D_success=true
else
    shear_3D_success=false
fi

cd $toplevel/src/testcases/3D/curvature/convergenceAnalysis/shear/shear_200
$toplevel/src/cl_advection -t 4 -w false -f strawberryField -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/strawberry/strawberry_200
$toplevel/src/cl_advection -t 4 -w false -f strawberryField -r 0.5 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/strawberry/strawberry_100
$toplevel/src/cl_advection -t 4 -w true -f strawberryField -r 0.25 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/strawberry/strawberry_50
$toplevel/src/cl_advection -t 4 -w false -f strawberryField -r 0.125 -o $toplevel/src/testcases/3D/curvature/convergenceAnalysis/strawberry/strawberry_25
cd $toplevel/src/testcaseScripts; ./strawberry_3D;
if [ $? -eq 0 ]; then
    strawberry_3D_success=true
else
    strawberry_3D_success=false
fi


if $shear_2D_success ; then
    echo "Testcase shear_2D succeeded."
else
    echo "Testcase shear_2D failed."
fi

if $navier_2D_success ; then
    echo "Testcase navier_2D succeeded."
else
    echo "Testcase navier_2D failed."
fi

if $shear_3D_success ; then
    echo "Testcase shear_3D succeeded."
else
    echo "Testcase shear_3D failed."
fi

if $strawberry_3D_success ; then
    echo "Testcase strawberry_3D succeeded."
else
    echo "Testcase strawberry_3D failed."
fi
