set datafile separator ","
set term pdf
set key above
set xlabel "time"


set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle_navier.pdf"

plot "navier_25/contactAngle.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/contactAngle.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/contactAngle.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/contactAngle.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "navier_25/contactAngle.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Contact Angle / deg"
set output "error_angle_navier.pdf"

plot "navier_25/contactAngle.csv" using 1:($2-$3) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/contactAngle.csv" using 1:($2-$3) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/contactAngle.csv" using 1:($2-$3) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/contactAngle.csv" using 1:($2-$3) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set output "results_position_navier.pdf"
set ylabel "Position"

plot "navier_25/position.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/position.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/position.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/position.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "navier_25/position.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Position"
set output "error_position_navier.pdf"

plot "navier_25/position.csv" using 1:($2-$3) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/position.csv" using 1:($2-$3) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/position.csv" using 1:($2-$3) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/position.csv" using 1:($2-$3) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set output "results_curvature_navier.pdf"
set ylabel "Curvature Divergence"

plot "navier_25/curvature.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/curvature.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/curvature.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/curvature.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "navier_25/curvature.csv" using 1:4 every 1 title "Reference" with line black,\
     "results_fs3d_postprocessing.txt" using 1:2 every 1 title "fs3d_curvature" with points pointsize 0.5,\

set ylabel "Error in Curvature (divergence)"
set output "error_curvature_divergence_navier.pdf"

plot "navier_25/curvature.csv" using 1:($2-$4) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/curvature.csv" using 1:($2-$4) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/curvature.csv" using 1:($2-$4) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/curvature.csv" using 1:($2-$4) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set ylabel "Error in Curvature (height)"
set output "error_curvature_height_navier.pdf"

plot "navier_25/curvature.csv" using 1:($3-$4) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/curvature.csv" using 1:($3-$4) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/curvature.csv" using 1:($3-$4) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/curvature.csv" using 1:($3-$4) every 16 title "dx = 0.0005" with points pointsize 0.5,\


