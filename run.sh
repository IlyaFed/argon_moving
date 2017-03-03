#!/bin/sh
# PBS -N planet_moving
cd /home/tmipt10/planet_moving

g++ main.cpp -o planet_count
step=10
while [ $step -le 100000 ]
do
./planet_count "$step"
cp energy_by_time\("$step"\).txt energy_by_time.txt
gnuplot gnuplot_energy_t.txt
mv graph.png energy_by_time\("$step"\).png
let step=$step*10
done

gnuplot gnuplot_energy_dtime.txt

./planet_count 10 withplot

gnuplot gnuplot_animate.txt
