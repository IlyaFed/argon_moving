// ./a.out dtime max_photo photo_step max_energy energy_step wall r_kr
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#define MAXNAME 10
#define GIST_STEPS 100

using namespace std;


//constant for Lennard_Jones
double mass_Ar = 39.948 / 1000/6.022e23; // moll
double sigma_Ar = 3.405*1e-10; //A
double epsilon_Ar = 119.8*1.23e-23; // *k_Bol
//


double E_max = 0;
double E_min = 0;

int N = 1000;


class zero_cond{
public:
	//double N_part = 1000;
	double dtime;
	double wall;
	int max_photo, photo_step, max_energy, energy_step, r_kr, range;

	zero_cond(const int argc, char* argv[]){
		if (argc > 8) {
		dtime = atoi(argv[1])*1.0e-16;
		max_photo = atoi(argv[2]);
		photo_step = atoi(argv[3]);
		max_energy = atoi(argv[4]);
		energy_step = atoi(argv[5]);
		wall = atoi(argv[6])*1.0e-10;
		r_kr = atoi(argv[7]);
		range = atoi(argv[8]);
		cout << "dtime = " << dtime << endl;
	} else {
		cerr << "Not enought arguments!" << endl;
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
		return to_string(coord_0[0]*1e10) + "\t" + to_string(coord_0[1]*1e10) + "\t" + to_string(coord_0[2]*1e10);
	}
};


class Sys_atom
{	
public:
	double E, T, r;
	Atom **a = new Atom* [N];	
	double sigma;
	double epsilon; 
	double mass;
	double dtime;
	double wall[3];
	double a_r;
	double acceleration;
	double x[3];

	double r_kr; //max rad of interraction

	Sys_atom(const double* wall_get, int* number_of_atom, const double sigma_get, const double epsilon_get, const double mass_get, const double dtime_get, const double r_kr_get ){
		for (int l = 0; l < 3; l++) {
			wall[l] = wall_get[l];
			cout << "wall: " << wall[l] << endl;
		}
		sigma = sigma_get;
		epsilon = epsilon_get;
		mass = mass_get;
		dtime = dtime_get;
		r_kr = r_kr_get;
		printf("r_kr = %e\n", r_kr );
		//N = (int) aa*b*c;
	//	*a = new Atom* [N]; //create massive of atoms
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
			cerr << "initial complete!" << endl;
		} else {
			cerr << "initial fail!" << endl;
		}
	}
	~Sys_atom(){
		for (int k = 0; k < N; k++) delete[] a[k];
		delete[] a;
		cerr << "Memory clear!" << endl;
	}

	int displace(const int tag_photo, ofstream &foutput, const int time){
		//acceleration collect
		for (int i = 0; i < N; i++)  for (int l = 0; l < 3; l++) a[i]->accel[l] = 0;

		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				if (r  < r_kr){
					a_r = -24.0*epsilon*(2*pow(sigma/r,12.0) - pow(sigma/r, 6.0))/r/r/mass; // F/r
					
					for (int l = 0; l < 3; l++) {
						acceleration = a_r*(a[i]->coord_0[l] - a[j]->coord_0[l]);
						//fprintf(foutput, "%d, delta, accel = %e , dtime = %e \n",time, a[i]->accel[l], dtime);
						
						a[i]->accel[l] -= acceleration;
						a[j]->accel[l] += acceleration;
					}

				}
			
			}
		}
		//make displacement
		
		for (int i = 0; i < N; i ++) {
			for (int l = 0; l < 3; l++) {
				x[l] = 2*a[i]->coord_0[l] - a[i]->coord_p[l]+a[i]->accel[l]*pow(dtime,2.0); //collect new coordinate
				
				if (x[l] < 0.0) {
					x[l] = -x[l];
					a[i]->coord_0[l] = - a[i]->coord_0[l];
				}
				if (x[l] > wall[l]) {
					x[l] = 2*wall[l] - x[l];
					a[i]->coord_0[l] = 2*wall[l] - a[i]->coord_0[l];
				}	
			}
			
			if (tag_photo ){
				foutput << a[i]->print() << time << endl;
			}
			for (int l = 0; l < 3; l++) {
				a[i]->coord_p[l] = a[i]->coord_0[l];
				a[i]->coord_0[l] = x[l];
			}
		}
		return 0;
	}

	int energy_collect( ofstream &fenergy, const int time, const bool fluct_flag)
		//collect Energy quantity
		{
		E = 0;
		T = 0.0;
		for (int i = 0; i < N; i++){
			T += mass*(pow(a[i]->coord_0[0]-a[i]->coord_p[0],2.0) + pow(a[i]->coord_0[1]-a[i]->coord_p[1],2.0) + pow(a[i]->coord_0[2]-a[i]->coord_p[2],2.0))/pow(dtime,2.0); // *
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				//if (r  < r_kr){
					E -= 4*epsilon*(pow(sigma/r,12) - pow(sigma/r, 6)); // F/r	
			//	}
			}
		}
		E += T;
		
		T = 2.0/3/1.23e-23*T/N;
	//	fprintf(stdout, "%e \t %d %lf\n", E, time, T);
		fenergy << E << "\t" << T << "\t" << time << endl;
		if (fluct_flag) {
			if ( (E_max == 0) || (E_max < E) ) E_max = E;
			if ( (E_min == 0) || (E_min > E) ) E_min = E;
		}
		return 0;
	}

	int g_r(double* mass){
		//acceleration collect
		double max_len = wall[0];// pow( pow(wall[0],2.0) + pow(wall[0],2.0) + pow(wall[0],2.0),0.5);

		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				r = pow( pow(a[j]->coord_0[0]-a[i]->coord_0[0],2.0) + pow(a[j]->coord_0[1]-a[i]->coord_0[1],2.0) + pow(a[j]->coord_0[2]-a[i]->coord_0[2],2.0),0.5); // distance between atoms
				mass[int(r/max_len*GIST_STEPS)] += pow(max_len,3)/N/4/3.14/pow(r,2)*1e10;
			}
			
		}
		//for (int i = 0; i < GIST_STEPS; i++) mass[i] = mass[i]/pow(i,2);

		return 0;
		}
		

	
};

int main(int argc, char* argv[]){
	ofstream foutput, fenergy;

	//create g(r)
	double g_r[GIST_STEPS];
	for (int i = 0; i < GIST_STEPS; i++) g_r[i] = 0;
	//

	zero_cond initial(argc, argv); // collect zero conditions

	if (initial.max_photo) foutput.open("term/output_photo.txt");
	// dependense Energy by time

	string name_file = "term/energy_by_time("+to_string(int(initial.dtime*1e16))+","+to_string(int(initial.r_kr))+","+to_string(int(initial.wall*1e10))+").txt" ;

	if (initial.max_energy) fenergy.open(name_file.c_str());

	//start
	double wall[3] = {initial.wall, initial.wall, initial.wall};
	N = pow(initial.range,3);
	int number_of_atom[3] = {initial.range, initial.range, initial.range};

	Sys_atom *argon = new Sys_atom(wall, number_of_atom, sigma_Ar, epsilon_Ar, mass_Ar, initial.dtime, initial.r_kr*sigma_Ar);
	
	
	int step = 0;
	int percentage = 0;
	int max_step = max(initial.max_photo, initial.max_energy);
	bool plot = 0;
	if (initial.max_photo) plot = 1;

	while (step <= max_step){
		if ( initial.max_photo && (!(step % initial.photo_step)) && (step <= initial.max_photo)) argon->displace(1, foutput, step);
		else argon->displace(0, foutput, step);
		if ( initial.max_energy && (!(step % initial.energy_step)) && (step <= initial.max_energy)) if (max_step-step > 10) argon->energy_collect(fenergy, step, 0);
			else argon->energy_collect(fenergy, step, 1);

		if (step == max_step) argon->g_r(g_r);

		if (100.0*step/max_step > percentage){
			printf("%d %% \n", percentage);
			percentage += 5;
		}


		step++;
	}
/*
	FILE *energy_fluct;
	energy_fluct = fopen("energy_fluct.txt", "a+");
	fprintf(energy_fluct, "%d %lf\n", dtime, E_max-E_min);
	fclose(energy_fluct);
*/

	foutput.close();
	fenergy.close();
	
	//write fluctuation of energy
	ofstream ffluct; //fluct of energy
	if (initial.max_energy) ffluct.open("term/energy_fluct.txt", ios::app);
	ffluct << initial.dtime << "\t" << E_max-E_min << endl;
	ffluct.close();

	//write g(r)
	ofstream g_rfile;
	name_file = "term/g(r)_("+to_string(int(initial.dtime*1e16))+","+to_string(int(initial.r_kr))+","+to_string(int(initial.wall*1e10))+").txt";
	g_rfile.open(name_file.c_str());
	for (int i = 0; i < GIST_STEPS; i++) g_rfile <<  i << "\t" << g_r[i] << endl;
	g_rfile.close();

	delete argon;
	cerr << "Everything ok!\n Finish!" << endl;
	return 0;
}