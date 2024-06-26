set datafile separator ","
set term pdf
set key above
set xlabel "time"


set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle_quadratic.pdf"

plot "quadratic_25/contactAngle.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/contactAngle.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/contactAngle.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/contactAngle.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "quadratic_25/contactAngle.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Contact Angle / deg"
set output "error_angle_quadratic.pdf"

plot "quadratic_25/contactAngle.csv" using 1:($2-$3) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/contactAngle.csv" using 1:($2-$3) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/contactAngle.csv" using 1:($2-$3) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/contactAngle.csv" using 1:($2-$3) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set output "results_position_quadratic.pdf"
set ylabel "Position"

plot "quadratic_25/position.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/position.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/position.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/position.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "quadratic_25/position.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Position"
set output "error_position_quadratic.pdf"

plot "quadratic_25/position.csv" using 1:($2-$3) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/position.csv" using 1:($2-$3) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/position.csv" using 1:($2-$3) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/position.csv" using 1:($2-$3) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set output "results_curvature_divergence_quadratic_.pdf"
set ylabel "Curvature (Divergence)"

plot "quadratic_25/curvature.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/curvature.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/curvature.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/curvature.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "quadratic_25/curvature.csv" using 1:4 every 1 title "Reference" with line black

set ylabel "Error in Curvature (Divergence)"
set output "error_curvature_divergence_quadratic.pdf"

plot "quadratic_25/curvature.csv" using 1:($2-$4) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/curvature.csv" using 1:($2-$4) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/curvature.csv" using 1:($2-$4) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/curvature.csv" using 1:($2-$4) every 16 title "dx = 0.0005" with points pointsize 0.5,\

set output "results_curvature_height_quadratic.pdf"
set ylabel "Curvature (Height)"

plot "quadratic_25/curvature.csv" using 1:3 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/curvature.csv" using 1:3 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/curvature.csv" using 1:3 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/curvature.csv" using 1:3 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "quadratic_25/curvature.csv" using 1:4 every 1 title "Reference" with line black

set ylabel "Error in Curvature (Height)"
set output "error_curvature_height_quadratic.pdf"

plot "quadratic_25/curvature.csv" using 1:($3-$4) every 2 title "dx = 0.004" with points pointsize 0.5,\
     "quadratic_50/curvature.csv" using 1:($3-$4) every 4 title "dx = 0.002" with points pointsize 0.5,\
     "quadratic_100/curvature.csv" using 1:($3-$4) every 8 title "dx = 0.001" with points pointsize 0.5,\
     "quadratic_200/curvature.csv" using 1:($3-$4) every 16 title "dx = 0.0005" with points pointsize 0.5,\


