set datafile separator ","
set term pdf
set key above

set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set style line 2 \
    linecolor rgb '#000000' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5 \
    dashtype 4

nplot = 1
#####

f(x) = a*x

#### Plot error ###
set ylabel "Maximum Position Error"
set xlabel "dx"

set output "max_error_position.pdf"

fit f(x) "max_error_position.csv" using 1:2 via a

plot [] [0:] "max_error_position.csv" using 1:2 title "Level Set" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2

#### Plot error ###
set ylabel "Maximum Contact Angle Error (Deg)"
set xlabel "dx"

set output "max_error_contactAngle.pdf"

fit f(x) "max_error_contactAngle.csv" using 1:2 via a

plot [] [0:] "max_error_contactAngle.csv" using 1:2 title "Level Set" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2
      
#####
set output "max_error_curvature.pdf"
set ylabel "Maximum Curvature Error"

fit f(x) "max_error_curvature.csv" using 1:2 via a

plot [] [0:] "max_error_curvature.csv" using 1:2 title "Level Set" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2

#####
set output "max_gradient_deviation.pdf"
set ylabel "max(||grad phi|-1|)"

fit f(x) "max_gradient_deviation.csv" using 1:2 via a

plot [] [0:] "max_gradient_deviation.csv" using 1:2 title "Level Set" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2

