all:
	g++ -std=c++11 main.cpp -o argon_moving 
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 10000 100 10000 50 60 8 5 3>term/log.txt
	gnuplot state_plot.txt

big_run:
	bash run.sh

hard:
	g++ -std=c++11 main.cpp -o argon_moving 
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr (N)^{1/3}
	./argon_moving 100 200000 100 0000 50 114 10 10 3>log.txt
	