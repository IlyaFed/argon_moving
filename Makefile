all:
	g++ main.cpp -o argon_moving
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr
	./argon_moving 100 0 100 10000 50 60 10
	#gnuplot gnuplot_energy_t.txt
	
big_run:
	g++ main.cpp -o planet_count
	step="0"
	while [ "$step" -lt "100" ]
	do
	echo "$step"
	step=$($step+10)
	done
