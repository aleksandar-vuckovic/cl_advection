set datafile separator ","
set term pdf
set key above
set xlabel "time"


set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle / deg"
set output "results_angle.pdf"

plot "navier_25/contactAngle.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/contactAngle.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/contactAngle.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/contactAngle.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "navier_25/contactAngle.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Contact Angle / deg"
set output "error_angle.pdf"

plot "navier_25/contactAngle.csv" using 1:($2-$3) every 2 title "LevelSet 25" with points pointsize 0.5,\
     "navier_50/contactAngle.csv" using 1:($2-$3) every 4 title "LevelSet 50" with points pointsize 0.5,\
     "navier_100/contactAngle.csv" using 1:($2-$3) every 8 title "LevelSet 100" with points pointsize 0.5,\
     "navier_200/contactAngle.csv" using 1:($2-$3) every 16 title "LevelSet 200" with points pointsize 0.5,\

set output "results_position.pdf"
set ylabel "Position"

plot "navier_25/position.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     "navier_50/position.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     "navier_100/position.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     "navier_200/position.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     "navier_25/position.csv" using 1:3 every 1 title "Reference" with line black

set ylabel "Error in Position"
set output "error_position.pdf"

plot "navier_25/position.csv" using 1:($2-$3) every 2 title "LevelSet 25" with points pointsize 0.5,\
     "navier_50/position.csv" using 1:($2-$3) every 4 title "LevelSet 50" with points pointsize 0.5,\
     "navier_100/position.csv" using 1:($2-$3) every 8 title "LevelSet 100" with points pointsize 0.5,\
     "navier_200/position.csv" using 1:($2-$3) every 16 title "LevelSet 200" with points pointsize 0.5,\

