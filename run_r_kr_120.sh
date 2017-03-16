#!/bin/sh
# PBS -N Number
cd /home/tmipt10/argon_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

	
	mkdir -p all

# influence of number of particles
mkdir -p all/r_kr
cd all/r_kr

# wall = 120
mkdir wall_120
cd wall_120

for step in 2 6 10
do

	mkdir -p r_kr_"$step"
	cd r_kr_"$step"
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./../../../../argon_moving 100 50000 100 50000 50 120 "$step" 10 >term/log.txt
	gnuplot ./../../../../state_plot.txt
	gnuplot ./../../../../animate_plot.txt

	cd ..
	
	let step=$step+1
done
cd ..