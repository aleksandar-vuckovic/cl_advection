set datafile separator ","
set term pdf
set key above
set xlabel "time"

testcase = "initEllipsoid"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle_".testcase.".pdf"

plot testcase."_25/contactAngle.csv" using 1:2 title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/contactAngle.csv" using 1:2 title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/contactAngle.csv" using 1:2 title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:2 title "dx = 0.0005" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:3 title "Reference" with line black

set ylabel "Error in Contact Angle / deg"
set output "error_angle_".testcase.".pdf"

plot testcase."_25/contactAngle.csv" using 1:(abs($2-$3)) title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/contactAngle.csv" using 1:(abs($2-$3)) title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/contactAngle.csv" using 1:(abs($2-$3)) title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:(abs($2-$3)) title "dx = 0.0005" with points pointsize 0.5,\

set output "results_position_".testcase.".pdf"
set ylabel "Position"

plot testcase."_25/position.csv" using 1:2 title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/position.csv" using 1:2 title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/position.csv" using 1:2 title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/position.csv" using 1:2 title "dx = 0.0005" with points pointsize 0.5,\
     testcase."_200/position.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Position"
set output "error_position_".testcase.".pdf"

plot testcase."_25/position.csv" using 1:(abs($2-$3)) title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/position.csv" using 1:(abs($2-$3))title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/position.csv" using 1:(abs($2-$3)) title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/position.csv" using 1:(abs($2-$3)) title "dx = 0.0005" with points pointsize 0.5,\

set output "results_curvature_".testcase.".pdf"
set ylabel "Curvature"

plot testcase."_25/curvature.csv" using 1:2 title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/curvature.csv" using 1:2 title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/curvature.csv" using 1:2 title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:2 title "dx = 0.0005" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Curvature"
set output "error_curvature_".testcase.".pdf"

plot testcase."_25/curvature.csv" using 1:(abs($2-$3)) title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/curvature.csv" using 1:(abs($2-$3)) title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/curvature.csv" using 1:(abs($2-$3)) title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:(abs($2-$3)) title "dx = 0.0005" with points pointsize 0.5,\

