set datafile separator ","
set term pdf
set key above



set style line 128 \
    linecolor rgb '#0060ad' \
    linetype 1 linewidth 1.5 \
    pointtype 2 pointsize 0.5

set style line 256 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1.5 \
    pointtype 5 pointsize 0.5

set style line 512 \
    linecolor rgb '#29c524' \
    linetype 1 linewidth 1.5 \
    pointtype 7 pointsize 0.5

set style line 1024 \
    linecolor rgb '#7D72F9' \
    linetype 1 linewidth 1.5 \
    pointtype 3 pointsize 0.5

set style line 2048 \
    linecolor 'orange' \
    linetype 1 linewidth 1.5 \
    pointtype 9 pointsize 0.5

set style line 1 \
    linecolor rgb '#000000' \
    linetype 1 linewidth 1.5 \
    pointtype 4 pointsize 0.5
    
set style line 2 \
    linecolor rgb '#000000' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5 \
    dashtype 4

set colorsequence default # default, podo

set output "results_gradient_norm.pdf"

set xlabel "time"
set ylabel "|grad phi|"

plot for [mesh in "100 200 400"] '0/'.mesh.'/gradientNormAtContactPoint.csv' using ($1):($2) title 'Source-0/Mesh-'.mesh with points pointsize 0.3,\
     for [mesh in "100 200 400"] '1/'.mesh.'/gradientNormAtContactPoint.csv' using ($1):($2) title 'Source-1/Mesh-'.mesh with points pointsize 0.3 


############################

   
set output "results.pdf"

do for [quantity in "position contactAngle curvature"]{

do for [source in "0 1"]{

set xlabel "time"

set ylabel quantity

plot for [mesh in "100 200 400"] source.'/'.mesh.'/'.quantity.'.csv' using ($1):(abs($2)) title source.'/'.mesh with line

}
}
############################

do for [quantity in "position contactAngle curvature"]{

set output 'comparison_'.quantity.'.pdf'

set xlabel "time"

set ylabel quantity

plot '0/400/'.quantity.'.csv' using ($1):($2) every 10 title "Source off" with points pointsize 0.7,\
     '1/400/'.quantity.'.csv' using ($1):($2) every 10 title "Source on" with points pointsize 0.6,\
     '1/400/'.quantity.'.csv' using ($1):($3) title "Reference" with line linestyle 2
}

set output "comparison_gradientNormAtContactPoint.pdf"

set ylabel "|grad {/symbol f}|"

plot '0/400/gradientNormAtContactPoint.csv' using ($1):($2) every 1 title "Source off" with points pointsize 0.7,\
     '1/400/gradientNormAtContactPoint.csv' using ($1):($2) every 1 title "Source on" with points pointsize 0.6

    
    
####################

set logscale y 10
set format y "10^{%L}"
set logscale x

do for [quantity in "position contactAngle curvature"]{

do for [source in "0 1"]{

set output 'error_'.quantity.'_'.source.'.pdf'

set xlabel "time"

set ylabel quantity." error"

plot for [mesh in "100 200 400"] source.'/'.mesh.'/'.quantity.'.csv' using ($1):(abs($2-$3)) title 'Source-'.source.'/Mesh-'.mesh with line

}
}



