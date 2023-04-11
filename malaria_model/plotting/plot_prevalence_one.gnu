#DIRECTORIES = system("ls -dl ../batch/S01*/")

#set term png small
#set output "prevalence_10.png"

set terminal qt size 1300,600 font ",7"

tt = 1

mobility = 1-tt*0.1
each_title = sprintf("Mobility = %.1f", mobility)

#set origin 0.0,(6.0-(tt*0.6))
#set size 0.9, 0.5


set title each_title
set datafile separator ";"
set xrange [730:1460]
set yrange [0:20]
set format y "%.0f%%"


set ylabel "Prevalence"


set xtics ('MDA Start' 730, 'MDA Finish' 994, '+1Y' 1095, '+2Y' 1460)
set xlabel 'Time'
set grid
FILES = system("ls -1 ../batch/S01E0".(tt+1)."/outputs/*_prevalence_daily.csv")
plot for [file in FILES] file every ::730::1460 using 1:(100*$6) with lines lw 1 lc rgb "#99222222" notitle, for [file in FILES] file every ::730::1460 using 1:(100*$3) with lines lw 1 lc rgb "#99662222" notitle

