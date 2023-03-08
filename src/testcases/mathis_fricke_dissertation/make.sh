#!/bin/bash

cd ../../
make
cd testcases/mathis_fricke_dissertation

cl_advection="../../../../cl_advection"

mkdir results

cd navier_2D
find -name "navier*" | while read line; do cd $line; $cl_advection; cd ..; done
gnuplot convergenceAnalysis.gp
# python (xx) ERRORPLOT
cp "fig_6.2(a).pdf" "../results/fig_6.2(a).pdf"
cp "fig_6.2(b).pdf" "../results/fig_6.2(b).pdf"
#und andere cps für die anderen beiden plots
cd ..

cd shear_3D
find -name "shear*" | while read line; do cd $line; $cl_advection; cd ..; done
# gnuplot (xx) GNUPLOT
#und entsprechende cps
cd ..
