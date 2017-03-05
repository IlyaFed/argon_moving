#!/bin/sh
# PBS -N planet_moving
#cd /home/tmipt10/planet_moving

g++ -std=c++11 main.cpp -o argon_moving 
echo "compilation complete"

step=2
while [ $step -le 10 ]
do
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 20000 100 0000 10 60 "$step" 10 3>log.txt
	gnuplot gnuplot_animate.txt

	mv image/graph.gif image/animate\(100,"$step",60\).png
	let step=$step+1
done

step=20
while [ $step -le 100 ]
do
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 0000 100 20000 10 "$step" 10 10 3>log.txt
	gnuplot gnuplot_animate.txt

	mv image/graph.gif image/animate\(100,10,"$step"\).png

	let step=$step+10
done
# #small dtime only energy
# step=2
# while [ $step -le 10 ]
# do
# 	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
# 	./argon_moving 100 0000 100 20000 10 60 "$step" 10 3>log.txt
# 	cp term/energy_by_time\(100,"$step",60\).txt term/energy_by_time.txt
# 	cp term/g\(r\)_\(100,"$step",60\).txt term/g\(r\).txt
# 	gnuplot gnuplot_energy_t.txt

# 	mv image/energy_by_time.png image/energy_by_time\(100,"$step",60\).png
# 	mv image/g\(r\).png image/g\(r\)_\(100,"$step",60\).png
# 	let step=$step+1
# done

# step=20
# while [ $step -le 100 ]
# do
# 	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
# 	./argon_moving 100 0000 100 20000 10 "$step" 10 10 3>log.txt
# 	cp term/energy_by_time\(100,10,"$step"\).txt term/energy_by_time.txt
# 	cp term/g\(r\)_\(100,10,"$step"\).txt term/g\(r\).txt
# 	gnuplot gnuplot_energy_t.txt

# 	mv image/energy_by_time.png image/energy_by_time\(100,10,"$step"\).png
# 	mv image/g\(r\).png image/g\(r\)_\(100,10,"$step"\).png

# 	let step=$step+10
# done

# # #large dtime only energy
# step=1000
# while [ $step -le 10000 ]
# do
# 	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
# 	./argon_moving "$step" 0000 100 10000 10 60 10 10 3>log.txt
# 	cp energy_by_time\("$step"\).txt energy_by_time.txt
# 	gnuplot gnuplot_energy_t.txt
# 	mv energy_by_time.png energy_by_time\("$step"\).png
# 	let step=$step+500
# done

# gnuplot gnuplot_energy_dtime.txt

#./planet_count 10 withplot

#gnuplot gnuplot_animate.txt
