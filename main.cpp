// ./a.out dtime max_photo photo_step max_energy energy_step wall r_kr
#include <stdlib.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <cmath>
#define MAXNAME 10

using namespace std;


//constant for Lennard_Jones
double mass_Ar = 39.948 / 1000/6.022e23; // moll
double sigma_Ar = 3.405*1e-10; //A
double epsilon_Ar = 119.8*1.23e-23; // *k_Bol
//

double E;
double E_max = 0;
double E_min = 0;

int N = 1000;


class zero_cond{
public:
	//double N_part = 1000;
	double dtime;
	double wall;
	int max_photo, photo_step, max_energy, energy_step, r_kr;

	zero_cond(const int argc, char* argv[]){
		if (argc > 7) {
		dtime = atoi(argv[1])*1.0e-16;
		max_photo = atoi(argv[2]);
		photo_step = atoi(argv[3]);
		max_energy = atoi(argv[4]);
		energy_step = atoi(argv[5]);
		wall = atoi(argv[6])*1.0e-12;
		r_kr = atoi(argv[7]);
		printf("dtime = %5.2f fm\n", dtime*1e15);
	} else {
		fprintf(stderr, "Not enought arguments!\n");
	}
	}
};

class Atom {
public:
	double coord_p[3], coord_0[3]; // dynamic quantity in m
	double accel[3];
	Atom(const double* coord_of, const double* displace){
		for( int i = 0; i < 3; i++) coord_0[i] = coord_of[i]; 
		for( int i = 0; i < 3; i++) coord_p[i] = coord_0[i]-displace[i]; 
	}
	string print(){
		char buf[100];
		sprintf(buf, "%lf \t %lf \t %lf", coord_0[0]*1e10, coord_0[1]*1e10, coord_0[2]*1e10);
		string string_planet = buf;
		return string_planet;
	}
};


class Sys_atom
{	
public:

	Atom **a = new Atom* [N];	
	double sigma;
	double epsilon; 
	double mass;
	int dtime;
	double wall[3];

	double r_kr; //max rad of interraction

	Sys_atom(const double* wall_get, int* number_of_atom, const double sigma_get, const double epsilon_get, const double mass_get, const int dtime_get, const double r_kr_get ){
		for (int l = 0; l < 3; l++) wall[l] = wall_get[l];
		sigma = sigma_get;
		epsilon = epsilon_get;
		mass = mass_get;
		dtime = dtime_get;
		r_kr = r_kr_get;
		//N = (int) aa*b*c;
		//*a = new Atom* [N]; //create massive of atoms
		double x[3];
		double v[3] = {0, 0, 0}; //initial velocity
		double steps[3];
		for (int l = 0; l < 3; l++) steps[l] = wall[l]/number_of_atom[l];


		int k = 0; //counter of atom
		for (x[0] = steps[0]/2; x[0] < wall[0]; x[0] += steps[0]){
			for (x[1] = steps[1]/2; x[1] < wall[1]; x[1] += steps[1]){
				for (x[2] = steps[2]/2; x[2] < wall[2]; x[2] += steps[2]){
					a[k] = new Atom(x,v);
					k++;
				}
			}
		}
		if (k == N){
			fprintf(stderr, "initial complete!\n");
		} else {
			fprintf(stderr, "initial fail!\n");
		}
	}
	~Sys_atom(){
		for (int k = 0; k < N; k++) delete[] a[k];
		delete[] a;
		fprintf(stderr, "Memory clear!\n");
	}

	int displace(const int tag_photo, FILE* foutput, const int time){
		//acceleration collect
		for (int i = 0; i < N; i++)  for (int j = 0; j < 3; j++) a[i]->accel[j] = 0;

		double r;
		double F_r;
		double acceleration;
		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				if (r  < r_kr){
					F_r = 24*epsilon*(2*pow(sigma/r,12) - pow(sigma/r, 6))/r/r/mass; // F/r
					
					for (int l = 0; l < 3; l++) {
						acceleration = F_r*(a[i]->coord_0[l] - a[j]->coord_0[l]);
						a[i]->accel[l] -= acceleration;
						a[j]->accel[l] += acceleration;
					}

				}
			
			}
		}
		//make displacement
		double x[3];
		
		for (int i = 0; i < N; i ++) {
			for (int l = 0; l < 3; l++) {
				x[l] = 2*a[i]->coord_0[l] - a[i]->coord_p[l]+a[i]->accel[l]*dtime*dtime; //collect new coordinate
				
				if (x[l] < 0) {
					x[l] = -x[l];
					a[i]->coord_0[l] = - a[i]->coord_0[l];
				}
				if (x[l] > 1) {
					x[l] = 2*wall[l] - x[l];
					a[i]->coord_0[l] = 2*wall[l] - a[i]->coord_0[l];
				}	
			}
			
			if (tag_photo ){
				fprintf(foutput, "%d %s \t %d \n", i, a[i]->print().c_str(), time);
			}
			for (int l = 0; l < 3; l++) {
				a[i]->coord_p[l] = a[i]->coord_0[l];
				a[i]->coord_0[l] = x[l];
			}
		}
		return 0;
	}

	int energy_collect(FILE* fenergy, const int time)
		//collect Energy quantity
		{
		E = 0;
		double T = 0;
		double r;
		for (int i = 0; i < N; i++){
			T += mass*(pow(a[i]->coord_0[0]-a[i]->coord_p[0],2.0) + pow(a[i]->coord_0[1]-a[i]->coord_p[1],2.0) + pow(a[i]->coord_0[2]-a[i]->coord_p[2],2.0))/pow(dtime,2.0); // *
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				//if (r  < r_kr){
					E += 4*epsilon*(pow(sigma/r,12) - pow(sigma/r, 6)); // F/r	
			//	}
			}
		}
		E += T;
		
		T = 2.0/3/1.23e-23*T;

		fprintf(fenergy, "%lf\t%d %lf\n", E, time, T);

		if ( (E_max == 0) || (E_max < E) ) E_max = E;
		if ( (E_min == 0) || (E_min > E) ) E_min = E;

		return 0;
	}
	
};

int main(int argc, char* argv[]){
	FILE *finput, *foutput, *fenergy;


	zero_cond initial(argc, argv); // collect zero conditions

	if (initial.max_photo && ((foutput = fopen("output_photo.txt", "w") ) == 0)) { 
			fprintf(stderr, "Error! Can't create output_photo file!\n");
			return 1;
	}


	
	// dependense Energy by time

	char *name_file = new char [100];
	sprintf(name_file, "energy_by_time(%d).txt", int(initial.dtime*1e12));

	if ( (fenergy = fopen(name_file, "w")) == 0) { 
		fprintf(stderr, "Error! Can't create energy file!\n");
		return 1;
	}
	fprintf(stdout, "%s create!\n", name_file);
	delete[] name_file;
	fprintf(fenergy, "E\ttime\t T\n");

	//start
	double wall[3] = {initial.wall, initial.wall, initial.wall};

	int number_of_atom[3] = {10,10,10};
	Sys_atom *argon = new Sys_atom(wall, number_of_atom, sigma_Ar, epsilon_Ar, mass_Ar, initial.dtime, initial.r_kr*sigma_Ar);
	
	int step = 0;
	int percentage = 0;
	int max_step = max(initial.max_photo, initial.max_energy);
	bool plot = 0;
	if (initial.max_photo) plot = 1;

	while (step <= max_step){
		if ( initial.max_photo && (!(step % initial.photo_step)) && (step <= initial.max_photo)) argon->displace(1, foutput, step);
		else argon->displace(0, foutput, step);
		if ( initial.max_energy && (!(step % initial.energy_step)) && (step <= initial.max_energy)) argon->energy_collect(fenergy, step);

		step++;
		if (100.0*step/max_step > percentage){
			printf("%d %% \n", percentage);
			percentage += 5;
		}
	}
/*
	FILE *energy_fluct;
	energy_fluct = fopen("energy_fluct.txt", "a+");
	fprintf(energy_fluct, "%d %lf\n", dtime, E_max-E_min);
	fclose(energy_fluct);
*/
	fclose(foutput);
	fclose(fenergy);
	delete argon;
	return 0;
}