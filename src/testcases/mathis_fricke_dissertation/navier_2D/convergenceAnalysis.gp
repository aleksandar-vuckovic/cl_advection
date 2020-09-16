set datafile separator ","
set term pdf
set key above
set xlabel "time"

testcase = "navier"

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set ylabel "Contact Angle [deg]"
set output "fig_6.2(a)".".pdf"

plot testcase."_25/contactAngle.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/contactAngle.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/contactAngle.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/contactAngle.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     testcase."_25/contactAngle.csv" using 1:3 every 1 title "Reference" with line black

set output "fig_6.2(b)".".pdf"
set ylabel "Curvature"

plot testcase."_25/curvature.csv" using 1:2 every 2 title "dx = 0.004" with points pointsize 0.5,\
     testcase."_50/curvature.csv" using 1:2 every 4 title "dx = 0.002" with points pointsize 0.5,\
     testcase."_100/curvature.csv" using 1:2 every 8 title "dx = 0.001" with points pointsize 0.5,\
     testcase."_200/curvature.csv" using 1:2 every 16 title "dx = 0.0005" with points pointsize 0.5,\
     testcase."_25/curvature.csv" using 1:3 every 1 title "2D Reference" with line black
