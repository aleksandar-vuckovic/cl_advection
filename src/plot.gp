set datafile separator ","
set term pdf
set output "results_angle.pdf"
set key above
set xlabel "time / s"
set ylabel "Contact Angle / deg"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

plot "contactAngle.csv" using 1:2 title "LevelSet" with points pointsize 0.5,\
     "contactAngle.csv" using 1:3 title "reference" with points pointsize 0.1,\

set output "results_curvature.pdf"
set ylabel "curvature"

plot "curvature.csv" using 1:2 title "LevelSet" with points pointsize 0.5, \
     "curvature.csv" using 1:3 title "reference" with points pointsize 0.1 \

set output "results_position.pdf"
set ylabel "position"

plot "position.csv" using 1:2 title "LevelSet" with points pointsize 0.5, \
     "position.csv" using 1:3 title "reference" with points pointsize 0.1 \