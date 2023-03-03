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
g(x) = b*x

#### Plot error ###
set ylabel "Maximum Position Error"
set xlabel "{/Symbol D}x"

set format x "%2.1t{/Symbol \264}10^{%L}"

set logscale x
set logscale y

set xrange[0.002:0.012]
set xtics(0.0025, 0.005, 0.01)

set output "max_error_position.pdf"

fit f(x) "0/max_error_position.csv" using 1:2 via a
fit g(x) "1/max_error_position.csv" using 1:2 via b

plot  "0/max_error_position.csv" using 1:2 title "Source off" smooth unique with linespoints,\
      "1/max_error_position.csv" using 1:2 title "Source on" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2,\
      g(x) notitle with line linestyle 2

#### Plot error ###
set ylabel "Maximum Contact Angle Error (Deg)"

set output "max_error_contactAngle.pdf"

fit f(x) "0/max_error_contactAngle.csv" using 1:2 via a
fit g(x) "1/max_error_contactAngle.csv" using 1:2 via b

plot "0/max_error_contactAngle.csv" using 1:2 title "Source off" smooth unique with linespoints,\
      "1/max_error_contactAngle.csv" using 1:2 title "Source on" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2,\
      g(x) notitle with line linestyle 2
      
#####
set output "max_error_curvature.pdf"
set ylabel "Maximum Curvature Error"

fit f(x) "0/max_error_curvature.csv" using 1:2 via a
fit g(x) "1/max_error_curvature.csv" using 1:2 via b

plot "0/max_error_curvature.csv" using 1:2 title "Source off" smooth unique with linespoints,\
      "1/max_error_curvature.csv" using 1:2 title "Source on" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2,\
      g(x) notitle with line linestyle 2

#####
set output "max_gradient_deviation.pdf"
set ylabel "max(||grad {/symbol f}|-1|)"

fit f(x) "1/max_gradient_deviation.csv" using 1:2 via a

plot  "0/max_gradient_deviation.csv" using 1:2 title "Source off" smooth unique with linespoints,\
      "1/max_gradient_deviation.csv" using 1:2 title "Source on" smooth unique with linespoints,\
      f(x) title "1st order fit" with line linestyle 2

