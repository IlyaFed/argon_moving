all:
	g++ -std=c++11 main.cpp -o argon_moving 
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr
	./argon_moving 100 0000 100 10000 50 60 10 3>log.txt
	#gnuplot gnuplot_energy_t.txt
	
big_run:
	bash run.sh
