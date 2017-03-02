#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <cmath>
#define MAXNAME 10

using namespace std;


//constant for Lennard_Jones
double mass_Ar = 39.948; // moll
double sigma_Ar = 3.405; //A
double epsilon_Ar = 119.8; // *k_Bol
//

double E;
double E_max = 0;
double E_min = 0;

int N = 1000;
class Atom {
public:
	double coord_p[3], coord_0[3]; // dynamic quantity in billion km per second
	double accel[3];
	Atom(const double* coord_of, const double* displace){
		for( int i = 0; i < 3; i++) coord_0[i] = coord_of[i]; //conver to km
		for( int i = 0; i < 3; i++) coord_p[i] = coord_0[i]-displace[i]; //conver to km
	}
	string print(){
		char buf[100];
		sprintf(buf, "%lf \t %lf \t %lf", coord_0[0], coord_0[1], coord_0[2]);
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
	double scale;

	double r_kr; //max rad of interraction
	Sys_atom(double aa, double b, double cc, const double sigma_get, const double epsilon_get, const double mass_get, const int dtime_get, const double scale_get){
		
		scale = scale_get;
		sigma = sigma_get/scale;
		epsilon = epsilon_get*24*1.38*6.022/scale;

		mass = mass_get;
		dtime = dtime_get;
		printf("epsilon = %lf\n", epsilon);
		r_kr = 2*sigma;
		//N = (int) aa*b*c;
		//*a = new Atom* [N]; //create massive of atoms
		int k = 0;
		aa = 1/aa; 
		b = 1/b;
		cc = 1/cc;
		double x[3];
		double v[3] = {0, 0, 0};
		for (x[0] = aa/2; x[0] < 1; x[0] += aa){
			for (x[1] = b/2; x[1] < 1; x[1] += b){
				for (x[2] = cc/2; x[2] < 1; x[2] += cc){
					a[k] = new Atom(x,v);
					k++;
				}
			}
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
		double k;
		double acceleration;
		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				if (r  < r_kr){
					k = epsilon*(2*pow(sigma/r,12) - pow(sigma/r, 6))/r/r/mass/scale/10000000; // F/r
					
					for (int l = 0; l < 3; l++) {
						acceleration = k*(a[i]->coord_0[l] - a[j]->coord_0[l]);
						a[i]->accel[l] += acceleration;
						//fprintf(foutput,"[%d,%d] accel= %lf,\t r = %lf, \t time = %d, \t [a,b] = [%lf, %lf]\n", i,j,  a[i]->accel[l], r, time, a[i]->coord_0[l], a[j]->coord_0[l]);
						a[j]->accel[l] -= acceleration;
					}

				}
			
			}
		}
		//make displacement
		for (int i = 0; i < N; i ++) {
			double x[3];
			for (int l = 0; l < 3; l++) {
				x[l] = 2*a[i]->coord_0[l] - a[i]->coord_p[l]+a[i]->accel[l]*dtime*dtime; //collect new coordinate
				int k = 0;
				if (x[l] < 0) {
					x[l] = -x[l];
					a[i]->coord_0[l] = - a[i]->coord_0[l];
				}
				if (x[l] > 1) {
					x[l] = 2 - x[l];
					a[i]->coord_0[l] = 2 - a[i]->coord_0[l];
				}	
				//if (int(x[l]) % 2) x[l] = 1 - (x[l] - int(x[l]));
				//else x[l] -= int(x[l]); 
				/*
				if (x[l] < 0) x[l] = - x[l];
				while (x[l] > 1){
					k ++;
					if (x[l] > 1) x[l] = 1 - x[l];
					if (x[l] < 0) x[l] = - x[l];
					if (k > 1000) printf("x[l] = %lf", x[l]);
				}*/
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
/*
	int energy_collect(const int tag_energy, FILE* fenergy)
		//collect Energy quantity
		{
			E = 0;
			for (int i = 0; i < nplanet; i++){
				E += 2.985*sun_system[i]->mass*(pow(sun_system[i]->coord_0[0]-sun_system[i]->coord_p[0],2.0)+ pow(sun_system[i]->coord_0[1]-sun_system[i]->coord_p[1],2.0) )/pow(dtime,2.0); // *E30
				for (int j = i+1; j < nplanet; j++) {
					E -= 2377248.03/pow(pow(sun_system[j]->coord_0[0]-sun_system[i]->coord_0[0],2.0)+ pow(sun_system[j]->coord_0[1]-sun_system[i]->coord_0[1],2.0),0.5)*sun_system[j]->mass*sun_system[i]->mass;
				}
			}
			now_energy += step_energy;

			fprintf(fenergy, "%lf\t%d\n", E, time);

			if ( (E_max == 0) || (E_max < E) ) E_max = E;
			if ( (E_min == 0) || (E_min > E) ) E_min = E;
		}
		if (100.0*time/finish_time > percentage){
			printf("%d %% \n", percentage);
			percentage += 5;
		}
		//printf("Make one step %d\n",  time);

		time += dtime;
	} */
	
};
int main(int argc, char* argv[]){
	FILE *finput, *foutput, *fenergy;

	int dtime = 100; //*10^(-15)
	if (argc > 1) {
		dtime = atoi(argv[1]);
		printf("dtime = %d\n", dtime);
	}
/*
	bool plot = 0;
	if (argc > 2) {
		if (!strcmp(argv[2], "withplot")) {
			plot = 1;
			printf("with plot!r\n");
		}
	}*/

	if ((foutput = fopen("output.txt", "w") ) == 0) { 
			fprintf(stderr, "Error! Can't create output_photo file!\n");
			return 0;
	}
/*
	if ( (finput = fopen("zero_conditions.txt", "r")) == 0){
		fprintf(stderr, "Error! Can't read file!\n");
		return 0;
	}
	// data for visualizate planet moving
	if (plot) {
		char *name_file = new char [100];
		sprintf(name_file, "output_photo(%d).txt", dtime);

		if (plot && ((foutput = fopen(name_file, "w")) == 0)) { 
			fprintf(stderr, "Error! Can't create output_photo file!\n");
			return 0;
		}
		delete[] name_file;
		fprintf(foutput, "name\tx\ty\ttime\n");
	}


	int nplanet;
	if (fscanf(finput, "quantity_of_planet= %d\n", &nplanet) == 0) {
		fprintf(stderr, "Error! Can't read quantity of planet!\n");
		return 0;
	}


	int finish_time;
	if (fscanf(finput, "max_photo_time= %d\n", &finish_time) == 0) {
		fprintf(stderr, "Error! Can't read max_photo_time\n");
		return 0;
	}

	int photo_step;
	if (fscanf(finput, "step_photo= %d\n", &photo_step) == 0) {
		fprintf(stderr, "Error! Can't read photo_step!\n");
		return 0;
	}

	int max_energy_time = 0;
	if (fscanf(finput, "max_energy_time= %d\n", &max_energy_time) == 0) {
		fprintf(stderr, "Error! Can't read dtimesubl!\n");
		return 0;
	}

	int step_energy = 0;
	if (fscanf(finput, "step_energy= %d\n", &step_energy) == 0) {
		fprintf(stderr, "Error! Can't read dtimesubl!\n");
		return 0;
	}
*/
	/*
	// dependense Energy by time

	char *name_file = new char [100];
	sprintf(name_file, "energy_by_time(%d).txt", dtime);

	if ( (fenergy = fopen(name_file, "w")) == 0) { 
		fprintf(stderr, "Error! Can't create energy file!\n");
		return 0;
	}
	delete[] name_file;
	fprintf(fenergy, "E\ttime\n");

	char name[MAXNAME];
	double coord[2], mass, speed;
	Planet **sun_system = new Planet* [nplanet]; //create massive of class planet
	for (int i = 0; i < nplanet; i++){
		fscanf(finput, "%s %lf %lf %lf %lf", name, &mass, &coord[0], &coord[1], &speed);
		sun_system[i] = new Planet(name, mass, coord, speed);
	}

	//for (int i = 0; i < nplanet; i++) sun_system[i]->print(stdout);
	int time = 0;
	int now_photo = 0;
	int now_energy = 0;
	int percentage = 0;

	if (plot == 0) {
		finish_time = max_energy_time;
		now_photo = finish_time;
	}
*/
	int a = 10;
	double scale = 10*5.6;
	Sys_atom *argon = new Sys_atom(a, a, a, sigma_Ar, epsilon_Ar, mass_Ar, dtime, scale);
	
	int step = 0;
	int maxstep = 5000;
	int percentage = 0;
	while (step <= maxstep){
		if (step % 100) argon->displace(0, foutput, step);
		else argon->displace(1, foutput, step);
		step++;
		if (100.0*step/maxstep > percentage){
			printf("%d %% \n", percentage);
			percentage += 10;
		}
	}
/*
	FILE *energy_fluct;
	energy_fluct = fopen("energy_fluct.txt", "a+");
	fprintf(energy_fluct, "%d %lf\n", dtime, E_max-E_min);
	fclose(energy_fluct);*/
//	fclose(finput);
	fclose(foutput);
//	fclose(fenergy);
	delete argon;
	return 0;
}