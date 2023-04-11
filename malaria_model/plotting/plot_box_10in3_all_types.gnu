set terminal qt size 800,600 font ",7"
set multiplot layout 3, 1 title " All Types " font ",8" rowsfirst downwards


set datafile separator ";"
set style fill solid 0.25 border -1
set style boxplot outliers pointtype 7
set pointsize 0.5
set boxwidth 0.125
set style data boxplot




#set yrange [0.07:0.11]
set ytics 0.01
set ylabel "Prevalence at MDA Start"
set format y "%.2f%%"

set xrange [0:11]
set xtics ('MDA Start' 1, 'MDA Finish + 3 months' 2)
set xlabel " "
set xtics textcolor rgb "#ffffffff"

set grid

do for [tt=1:10] {

mobility = 1-tt*0.1
each_title = sprintf("%.1f", mobility)

set xtics add (each_title tt)

}


plot for [tt = 1:10] "../batch/S01E0".(tt+1)."/prevalence_at_730.out" using (tt):6 notitle






#set yrange [0.02:0.06]
set ytics 0.01
set ylabel "Prevalence at MDA Finish + 3 months"
set format y "%.2f%%"

set xrange [0:11]
set xlabel " "
#set xtics textcolor rgb "#black"

set grid

do for [tt=1:10] {

mobility = 1-tt*0.1
each_title = sprintf("%.1f", mobility)

set xtics add (each_title tt)

}


plot for [tt = 1:10] "../batch/S01E0".(tt+1)."/prevalence_at_1087.out" using (tt):6 notitle





#set yrange [0.03:0.11]
set ytics 0.01
set ylabel "Prevalence at MDA Finish + 1 year"
set format y "%.2f%%"

set xrange [0:11]
set xlabel "Reduced Mobility"
set xtics textcolor rgb "#black"

set grid

do for [tt=1:10] {

mobility = 1-tt*0.1
each_title = sprintf("%.1f", mobility)

set xtics add (each_title tt)

}


plot for [tt = 1:10] "../batch/S01E0".(tt+1)."/prevalence_at_1359.out" using (tt):6 notitle