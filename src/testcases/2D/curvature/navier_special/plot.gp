set datafile separator ","
set term pdf
set output "results_angle.pdf"
set key above
set xlabel "time"
set ylabel "Contact Angle / deg"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

plot "contactAngle.csv" using 1:2 every 10 title "LevelSet" with points pointsize 0.5,\
     "contactAngle.csv" using 1:3 every 10 title "Reference" with line

set output "results_position.pdf"
set ylabel "Position"

plot "position.csv" using 1:2 every 10 title "LevelSet" with points pointsize 0.5,\
     "position.csv" using 1:3 every 10 title "Reference" with line

set output "results_curvature.pdf"
set ylabel "Curvature"

plot "curvature.csv" using 1:2 every 10 title "LevelSet Divergence" with points pointsize 0.5, \
     "curvature.csv" using 1:3 every 10 title "LevelSet Height" with points pointsize 0.5, \
     "curvature.csv" using 1:4 every 10 title "Reference Euler" with line,\
     "kappa_analytic_reference.txt" using 1:2 every 10 title "analytic Reference" with points pointsize 0.5