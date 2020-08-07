#! /bin/bash

find -name "navier*" | while read line; do cd $line; ../../../../cl_advection; cd ..; done
gnuplot plot.gp