#!/bin/sh
# PBS -N Number
cd /home/tmipt10/argon_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

	
	mkdir -p all

# influence of number of particles
mkdir -p all/N_inf
# wall = 60
mkdir -p all/N_inf/wall_60
step=5
while [ $step -le 20 ]
do
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 20000 100 20000 50 60 8 "$step" 3>term/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/N_inf/wall_60/N_"$step"

	mv  image all/N_inf/wall_60/N_"$step"/
	mv  term all/N_inf/wall_60/N_"$step"/
	
	let step=$step+5
done

# wall = 120
mkdir -p all/N_inf/wall_120
step=5
while [ $step -le 20 ]
do
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 20000 100 20000 50 120 8 "$step" 3>term/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/N_inf/wall_120/N_"$step"

	mv  image all/N_inf/wall_120/N_"$step"/
	mv  term all/N_inf/wall_120/N_"$step"/
	
	let step=$step+5
done
