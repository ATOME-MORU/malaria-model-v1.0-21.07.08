set terminal qt size 400,600 font ",7"
set multiplot layout 5, 2 title " " font ",8" columnsfirst downwards


set datafile separator ";"
set style fill solid 0.25 border -1
#set style boxplot outliers pointtype 7
#set style data boxplot

#set yrange [0.02:0.06]
#set ytics 0.01

set xtics ('MDA Finish + 3 months' 1)

set grid

do for [tt=1:10] {

mobility = 1-tt*0.1
each_title = sprintf("Mobility = %.1f", mobility)

set title each_title offset 0,-1

#set style data boxplot

#plot "../batch/S01E0".(tt+1)."/prevalence_at_1087.out" using (1):6 with boxplot notitle

set jitter overlap first 2
set style data points

#set table $density
plot "../batch/S01E0".(tt+1)."/prevalence_at_1087.out" using (1):6 with lines



}