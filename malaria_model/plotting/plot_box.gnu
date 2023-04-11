set terminal qt size 400,600 font ",7"
set multiplot layout 5, 2 title " " font ",8" columnsfirst downwards


set datafile separator ";"
set style fill solid 0.25 border -1
set style boxplot outliers pointtype 7
set style data boxplot

set yrange [0.02:0.11]
set ytics 0.02

set xtics ('MDA Start' 1, 'MDA Finish + 3 months' 2)

set grid

#do for [tt=1:10] {

#mobility = 1-tt*0.1
#each_title = sprintf("Mobility = %.1f", mobility)

#set title each_title offset 0,-1

#plot "../batch/S01E0".(tt+1)."/prevalence_at_730.out" using (1):6 notitle, "../batch/S01E0".(tt+1)."/prevalence_at_1087.out" using (2):6 notitle

plot for [tt = 1:10] "../batch/S01E0".(tt+1)."/prevalence_at_1087.out" using (tt):6 notitle


#}