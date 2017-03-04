#!/bin/sh
# PBS -N planet_moving
#cd /home/tmipt10/planet_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

#small dtime only energy
step=10
while [ $step -le 1000 ]
do
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr
	./argon_moving "$step" 0000 100 10000 10 60 10 3>log.txt
	cp energy_by_time\("$step"\).txt energy_by_time.txt
	gnuplot gnuplot_energy_t.txt
	mv energy_by_time.png energy_by_time\("$step"\).png
	let step=$step+50
done

#large dtime only energy
step=1000
while [ $step -le 10000 ]
do
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr
	./argon_moving "$step" 0000 100 10000 10 60 10 3>log.txt
	cp energy_by_time\("$step"\).txt energy_by_time.txt
	gnuplot gnuplot_energy_t.txt
	mv energy_by_time.png energy_by_time\("$step"\).png
	let step=$step+500
done

gnuplot gnuplot_energy_dtime.txt

#./planet_count 10 withplot

#gnuplot gnuplot_animate.txt
