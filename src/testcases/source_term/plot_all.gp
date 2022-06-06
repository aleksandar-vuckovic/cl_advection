set datafile separator ","
set term pdf
set key above



set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set output "results_gradient_norm.pdf"

do for [source in "0 1"]{

set xlabel "time"
set ylabel "|grad phi|"

plot for [mesh in "50 100 200"] source.'/'.mesh.'/gradientNormAtContactPoint.csv' using ($1):($2) title source.'/'.mesh

}

    
####################

set logscale y
set logscale x

set output "results.pdf"

do for [source in "0 1"]{

do for [quantity in "position contactAngle curvature"]{

set xlabel "time"
set ylabel quantity

plot for [mesh in "50 100 200"] source.'/'.mesh.'/'.quantity.'.csv' using ($1):(abs($2-$3)) title source.'/'.mesh

}
}



