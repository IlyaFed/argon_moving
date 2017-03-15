#!/bin/sh
# PBS -N r_kr
cd /home/tmipt10/argon_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

	mkdir -p all

# influence of number of particles
mkdir -p all/r_kr
# wall = 60
mkdir -p all/r_kr/wall_60
step=2
while [ $step -le 10 ]
do
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 20000 100 20000 50 60 "$step" 10 3>term/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/r_kr/wall_60/r_kr_"$step"

	mv -r image all/r_kr/wall_60/r_kr_"$step"/
	mv -r term all/r_kr/wall_60/r_kr_"$step"/
	
	let step=$step+1
done

# wall = 120
mkdir -p all/r_kr/wall_120
step=2
while [ $step -le 10 ]
do
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 20000 100 20000 50 60 "$step" 10 3>term/log.txt
	gnuplot state_plot.txt
	gnuplot animate_plot.txt

	mkdir -p all/r_kr/wall_120/r_kr_"$step"

	mv -r image all/r_kr/wall_120/r_kr_"$step"/
	mv -r term all/r_kr/wall_120/r_kr_"$step"/
	
	let step=$step+1
done