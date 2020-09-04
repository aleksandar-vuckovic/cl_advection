set datafile separator ","
set term pdf
set key above
set xlabel "time"

testcase = "strawberry"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle_".testcase.".pdf"

plot testcase."_50/contactAngle.csv" using 1:2 every 4 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/contactAngle.csv" using 1:2 every 8 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:2 every 16 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/contactAngle.csv" using 1:2 every 32 title "dx = 0.0025" with points pointsize 0.5,\
     testcase."_400/contactAngle.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Contact Angle / deg"
set output "error_angle_".testcase.".pdf"

plot testcase."_50/contactAngle.csv" using 1:($2-$3) every 4 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/contactAngle.csv" using 1:($2-$3) every 8 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:($2-$3) every 16 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/contactAngle.csv" using 1:($2-$3) every 32 title "dx = 0.0025" with points pointsize 0.5,\
     

set output "results_position_".testcase.".pdf"
set ylabel "Position"

plot testcase."_400/position.csv" using 1:2 every 32 title "x, dx = 0.0025" with points pointsize 0.5,\
     testcase."_400/position.csv" using 1:4 every 32 title "z, dx = 0.0025" with points pointsize 0.5,\

set ylabel "Error in Position"
set output "error_position_".testcase.".pdf"

plot testcase."_50/position.csv" using 1:($2-$3) every 4 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/position.csv" using 1:($2-$3) every 8 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/position.csv" using 1:($2-$3) every 16 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/position.csv" using 1:($2-$3) every 16 title "dx = 0.0025" with points pointsize 0.5,\

set output "results_curvature_".testcase.".pdf"
set ylabel "Curvature"

plot testcase."_50/curvature.csv" using 1:2 every 1 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/curvature.csv" using 1:2 every 1 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:2 every 1 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/curvature.csv" using 1:2 every 1 title "dx = 0.0025" with points pointsize 0.5,\
     testcase."_400/curvature.csv" using 1:3 every 1 title "2D Reference" with line black

set ylabel "Error in Curvature"
set output "error_curvature_".testcase.".pdf"

plot testcase."_50/curvature.csv" using 1:($2-$3) every 4 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/curvature.csv" using 1:($2-$3) every 8 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:($2-$3) every 16 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/curvature.csv" using 1:($2-$3) every 32 title "dx = 0.0025" with points pointsize 0.5,\

set output "results_curvatureDerivative_".testcase.".pdf"
set ylabel "CurvatureDerivative"

plot testcase."_50/curvatureDerivative.csv" using 1:2 title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/curvatureDerivative.csv" using 1:2 title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/curvatureDerivative.csv" using 1:2 title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/curvatureDerivative.csv" using 1:2 title "dx = 0.0025" with points pointsize 0.5,\
     testcase."_200/curvatureDerivative.csv" using 1:3 title "Reference" with points pointsize 0.5,\


set output "error_curvatureDerivative_".testcase.".pdf"
set ylabel "Absolute Error CurvatureDerivative"
plot testcase."_50/curvatureDerivative.csv" using 1 : (abs($2-$3)) title "dx = 0.02" with points pointsize 0.5,\
     testcase."_100/curvatureDerivative.csv" using 1 : (abs($2-$3)) title "dx = 0.01" with points pointsize 0.5,\
     testcase."_200/curvatureDerivative.csv" using 1 : (abs($2-$3)) title "dx = 0.005" with points pointsize 0.5,\
     testcase."_400/curvatureDerivative.csv" using 1 : (abs($2-$3)) title "dx = 0.0025" with points pointsize 0.5,\