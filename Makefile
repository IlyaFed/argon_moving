all:
	g++ main.cpp -o argon_model
	./argon_model 10
	#cp energy_by_time\(100\).txt energy_by_time.txt
	gnuplot gnuplot_energy_t.txt
	
big_run:
	g++ main.cpp -o planet_count
	step="0"
	while [ "$step" -lt "100" ]
	do
	echo "$step"
	step=$($step+10)
	done
