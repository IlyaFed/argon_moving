#!/bin/sh
# PBS -N run_dtime
cd /home/tmipt10/argon_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

	mkdir -p all

# influence of number of particles
mkdir -p all/dtime
# wall = 60
mkdir -p all/dtime/wall_60
step=10
while [ $step -le 300 ]
do
	mkdir -p term_dtime
	mkdir -p image_dtime
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving "$step" 20000 100 20000 50 60 8 10 3>term_dtime/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/dtime/wall_60/dtime_"$step"

	mv -r image_dtime all/dtime/wall_60/dtime_"$step"/image
	mv -r term_dtime all/dtime/wall_60/dtime_"$step"/term
	
	let step=$step+50
done

# wall = 120
mkdir -p all/dtime/wall_120
step=10
while [ $step -le 300 ]
do
	mkdir -p term_dtime
	mkdir -p image_dtime
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving "$step" 20000 100 20000 50 60 8 10 3>term_dtime/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/dtime/wall_120/dtime_"$step"

	mv -r image_dtime all/dtime/wall_120/dtime_"$step"/image
	mv -r term_dtime all/dtime/wall_120/dtime_"$step"/term
	
	let step=$step+50
done