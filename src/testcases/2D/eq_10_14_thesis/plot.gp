set datafile separator ","
set term pdf
set key above
set xlabel "time"


set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle_eq_10_14_thesis.pdf"

plot "navier_32/contactAngle.csv" using 1:2 every 2 title "numY = 32" with points pointsize 0.5,\
     "navier_64/contactAngle.csv" using 1:2 every 4 title "numY = 64" with points pointsize 0.5,\
     "navier_128/contactAngle.csv" using 1:2 every 8 title "numY = 128" with points pointsize 0.5,\
     "navier_256/contactAngle.csv" using 1:2 every 16 title "numY = 256" with points pointsize 0.5,\
     "navier_512/contactAngle.csv" using 1:2 every 16 title "numY = 512" with points pointsize 0.5,\
     "navier_512/contactAngle.csv" using 1:3 every 1 title "Reference" with line black,\

set ylabel "Error in Contact Angle / deg"
set output "error_angle_eq_10_14_thesis.pdf"

plot "navier_32/contactAngle.csv" using 1:(abs($2-$3)) every 2 title "numY = 32" with points pointsize 0.5,\
     "navier_64/contactAngle.csv" using 1:(abs($2-$3)) every 4 title "numY = 64" with points pointsize 0.5,\
     "navier_128/contactAngle.csv" using 1:(abs($2-$3)) every 8 title "numY = 128" with points pointsize 0.5,\
     "navier_256/contactAngle.csv" using 1:(abs($2-$3)) every 16 title "numY = 256" with points pointsize 0.5,\
     "navier_512/contactAngle.csv" using 1:(abs($2-$3)) every 32 title "numY = 512" with points pointsize 0.5,\

set ylabel "Curvature (Divergence method)"
set output "results_curvature_eq_10_14_thesis.pdf"

plot "navier_32/curvature.csv" using 1:2 every 2 title "numY = 32" with points pointsize 0.5,\
     "navier_64/curvature.csv" using 1:2 every 4 title "numY = 64" with points pointsize 0.5,\
     "navier_128/curvature.csv" using 1:2 every 8 title "numY = 128" with points pointsize 0.5,\
     "navier_256/curvature.csv" using 1:2 every 16 title "numY = 256" with points pointsize 0.5,\
     "navier_512/curvature.csv" using 1:2 every 16 title "numY = 512" with points pointsize 0.5,\
     "navier_512/curvature.csv" using 1:3 every 1 title "Reference" with line black,\

set ylabel "Error in Curvature"
set output "error_curvature_eq_10_14_thesis.pdf"

plot "navier_32/curvature.csv" using 1:(abs($2-$3)) every 2 title "numY = 32" with points pointsize 0.5,\
     "navier_64/curvature.csv" using 1:(abs($2-$3)) every 4 title "numY = 64" with points pointsize 0.5,\
     "navier_128/curvature.csv" using 1:(abs($2-$3)) every 8 title "numY = 128" with points pointsize 0.5,\
     "navier_256/curvature.csv" using 1:(abs($2-$3)) every 16 title "numY = 256" with points pointsize 0.5,\
     "navier_512/curvature.csv" using 1:(abs($2-$3)) every 32 title "numY = 512" with points pointsize 0.5,\
