all:
	g++ -std=c++11 main.cpp -o argon_moving 
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 0000 100 4000 50 60 10 10 3>log.txt
	#gnuplot gnuplot_energy_t.txt
	cp term/energy_by_time\(100,10,60\).txt term/energy_by_time.txt
	cp term/g\(r\)_\(100,10,60\).txt term/g\(r\).txt
	gnuplot gnuplot_energy_t.txt

	mv image/energy_by_time.png image/energy_by_time\(100,10,60\).png
	mv image/g\(r\).png image/g\(r\)_\(100,10,60\).png

big_run:
	bash run.sh
