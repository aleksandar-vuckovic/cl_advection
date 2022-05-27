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

plot "source_on/contactAngle.csv" using 1:2 title "Source On" with line,\
     "source_off/contactAngle.csv" using 1:2 title "Source Off" with line,\
     "source_on/contactAngle.csv" using 1:3 title "Reference" with line

set output "results_position.pdf"
set ylabel "Position"

plot "source_on/position.csv" using 1:2 title "Source On" with line,\
     "source_off/position.csv" using 1:2 title "Source Off" with line,\
     "source_on/position.csv" using 1:3 title "Reference" with line
