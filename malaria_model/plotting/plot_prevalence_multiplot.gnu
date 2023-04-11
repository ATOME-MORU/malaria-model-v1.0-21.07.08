#DIRECTORIES = system("ls -dl ../batch/S01*/")

#set term png small
#set output "prevalence_10.png"

set terminal qt size 900,600 font ",7"

set tmargin 5
set bmargin 2

set multiplot layout 5, 2 title " " font ",8" columnsfirst downwards

do for [tt=1:10]{

mobility = 1-tt*0.1
each_title = sprintf("Mobility = %.1f", mobility)

set tics in

#set origin 0.0,(6.0-(tt*0.6))
#set size 0.9, 0.5

if (tt<6) {
	set ytics textcolor rgb "black"
} else {
	set ytics textcolor rgb "#ffffffff"
}

if (tt > 5) {
	set lmargin 5
	set rmargin 8
} else {
	set lmargin 8
	set rmargin 5
}



set title each_title offset 0,-3
set datafile separator ";"
set xrange [730:1460]
set yrange [0:20]
set format y "%.0f%%"

if (tt==3) {
	set ylabel "Prevalence"
} else {
	set ylabel " "
}

#set xlabel 'Time (End of Year #)'

set tmargin 0
set bmargin 2

if (tt==5 || tt==10) {
	#set xtics ('MDA Start' 730, 'MDA Finish' 994, '+1Y' 1095, '+2Y' 1460)
	#set xtics ('S' 730, 'F' 994, '+1M' 1025, '+3M' 1087, '+1Y' 1359, 'Y2' 1460)
	set xtics ('MDA Start' 730, 'Finish' 994, '+3M' 1087, '+1Y' 1359, 'Y2' 1460)

	set xtics textcolor rgb "#black"



	#set label 11 center at graph 0.2,char 1 "- MDA -" font ",7"

} else {
	#set xtics ('MDA Start' 730, 'MDA Finish' 994, '+1Y' 1095, '+2Y' 1460)
	#set xtics ('MDA Start' 730, 'MDA Finish' 994, '+1M' 1025, '+3M' 1087, '+1Y' 1359, 'Y2' 1460)
	set xtics ('MDA Start' 730, 'Finish' 994, '+3M' 1087, '+1Y' 1359, 'Y2' 1460)

	set xtics textcolor rgb "#ffffffff"

	unset xlabel

	#unset label 11
}

if (tt == 6) {
	set label 11 "All Types" at graph 0.85,0.6 font ",7" textcolor rgb "#33222222"
	set label 12 "RA" at graph 0.85,0.3 font ",7" textcolor rgb "#33662222"
} else {
	unset label 11
	unset label 12
}


set grid

FILES = system("ls -1 ../batch/S01E0".(tt+1)."/outputs/*_prevalence_daily.csv")
plot for [file in FILES] file every ::730::1460 using 1:(100*$6) with lines lw 1 lc rgb "#99222222" notitle, for [file in FILES] file every ::730::1460 using 1:(100*$3) with lines lw 1 lc rgb "#99662222" notitle



}



unset multiplot


