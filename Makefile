all:
	g++ -std=c++11 -ggdb lib/atom.cpp main.cpp -o argon_moving 
	mkdir -p term
	mkdir -p image
	# ./argon_moving dtime max_photo photo_step max_energy energy_step wall r_kr N
	./argon_moving 1 000000 1 500000 100000 60 8 32 3> term/log.txt
	gnuplot state_plot.txt

clean:
	rm -f argon_moving
	