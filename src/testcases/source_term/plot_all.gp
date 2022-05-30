set datafile separator ","
set term pdf
set output "results.pdf"
set key above
set xlabel "time"


set style line 1 \
    linecolor rgb '#dd181f' \
    linetype 1 linewidth 1 \
    pointtype 2 pointsize 0.5

set logscale y
set logscale x

do for [source in "0 1"]{

do for [quantity in "position contactAngle curvature"]{

set ylabel quantity

plot for [mesh in "10 25 50"] source.'/'.mesh.'/'.quantity.'.csv' using ($1):(abs($2-$3)) title source.'/'.mesh

#file = source.'/'.mesh.'/'.quantity.'.csv'

#plot file using ($1):(abs($2-$3)) title quantity.'/'.source.'/'.mesh


}
}

#list="0/25/contactAngle.csv 0/50/contactAngle.csv 1/25/contactAngle.csv 1/50/contactAngle.csv"
#plot for [file in list] file using ($1):(abs($2-$3))

