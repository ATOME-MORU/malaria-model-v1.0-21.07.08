#!/bin/bash
day_of_interest=1359

for ee in $(seq 1 11); do

	if [ -f ../batch/S01E0${ee}/prevalence_at_$day_of_interest.out ]; then
		rm ../batch/S01E0${ee}/prevalence_at_$day_of_interest.out
	fi

	for i in ../batch/S01E0${ee}/outputs/*_prevalence_daily.csv; do
	    awk -F ";" '{ if($1 == '$(($day_of_interest))') { print } }' $i >> ../batch/S01E0${ee}/prevalence_at_$day_of_interest.out

	done

	echo ${ee}"done"

done