// ./a.out dtime max_photo photo_step max_energy energy_step wall r_kr
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>


#define MAXNAME 10
#define GIST_STEPS 100
#define V_STEPS 100
#define MSD_STEP 400
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
		r_kr = min(r_kr, int(wall/sigma_Ar/2));
		range = atoi(argv[8]);
		cout << "dtime = " << dtime << endl;
	} else {
		cerr << "Not enought arguments!" << endl;
	}
	}
};


class Atom {
public:
	double coord_p[3], coord_0[3], coord_initial[3]; // dynamic quantity in m
	double accel[3], velocity_initial[3];
	int flag;
	Atom(const double* coord_of, const double* displace){
		for( int l = 0; l < 3; l++) {
			coord_0[l] = coord_of[l]; 
			coord_p[l] = coord_0[l]-displace[l]; 
			flag = 1;
			set_initial();
		}
	}

	int set_initial(){
		for( int l = 0; l < 3; l++) {
			coord_initial[l] = coord_0[l]; 
			velocity_initial[l] = (coord_0[l]-coord_p[l]); 
		}
		return 0;
	}
	string print(const double* wall){
		return to_string(coord_0[0]/wall[0]) + " " + to_string(coord_0[1]/wall[1]) + " " + to_string(coord_0[2]/wall[2]);
	}
};


class Sys_atom
{	
public:
	double E, T, P, r;
	Atom **a = new Atom* [N];	
	double sigma;
	double epsilon; 
	double mass;
	double dtime;
	double wall[3];
	double a_r;
	int k = 0;
	double acceleration;
	double x[3];
	double x_rand[3];
	double VC = 0; //velocity_correlation
	double MSD = 0;
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
		//N = (int) aa*b*c;
	//	*a = new Atom* [N]; //create massive of atoms
		double v[3] = {0, 0, 0}; //initial velocity
		double steps[3];
		for (int l = 0; l < 3; l++) steps[l] = wall[l]/number_of_atom[l];


		k = 0; //counter of atom
		for (x[0] = steps[0]/2; x[0] < wall[0]/2; x[0] += steps[0]){
			for (x[1] = steps[1]/2; x[1] < wall[1]; x[1] += steps[1]){
				for (x[2] = steps[2]/2; x[2] < wall[2]; x[2] += steps[2]){

					//too mach for save sum impuls = 0
					for (int l = 0; l < 3; l++) v[l] = (1.0*rand()/RAND_MAX - 0.5)*steps[l]/2000;
					a[k] = new Atom(x,v);
					if (k == 0) a[k]->flag = 3;
					k++;
					for (int l = 0; l < 3; l++) {
						v[l] = -v[l];
						x_rand[l] = x[l]; 
					}
					x_rand[0] = wall[0] - x[0];
					a[k] = new Atom(x_rand,v);
					k++;
					//cout << k <<endl;
				}
			}
		}
		if (k == N){
			cout << "initial complete!" << endl;
		} else {
			cout << "initial fail!" << endl;
		}
	}
	~Sys_atom(){
		for (int k = 0; k < N; k++) delete[] a[k];
		delete[] a;
		cerr << "Memory clear!" << endl;
	}

	double module(const double a, const double b) {
		if  (abs(a-b) < wall[0]/2.0) {
			return a-b;	                                       // TODO make it for all wall
		} else { 
			if (a < b) return a-b+wall[0];
			else return a-b-wall[0];
		} 
	}

	int displace(const int tag_photo, ofstream &foutput, const int time){
		//acceleration collect
		for (int i = 0; i < N; i++)  {
			for (int l = 0; l < 3; l++) a[i]->accel[l] = 0;
			a[i]->flag = 1;
		}
		a[0]->flag = 3;
		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_0[l], a[j]->coord_0[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0),0.5); // distance between atoms
				if (r  < r_kr){
					a_r = 24.0*epsilon*(2*pow(sigma/r,12.0) - pow(sigma/r, 6.0))/r/r/mass; // F/r/m
					if (i == 0) a[j]->flag = 2;
					for (int l = 0; l < 3; l++) {
						acceleration =  a_r*x[l];						
						a[i]->accel[l] += acceleration;
						a[j]->accel[l] -= acceleration;
					}

				}
			
			}
		}
		//make displacement
		k = 0;
		for (int i = 0; i < N; i ++) {
			for (int l = 0; l < 3; l++) {
				x[l] = 2*a[i]->coord_0[l] - a[i]->coord_p[l]+a[i]->accel[l]*pow(dtime,2.0); //collect new coordinate
				
				if (x[l] < 0.0) {
					x[l] = wall[l] + x[l];
					a[i]->coord_0[l] = wall[l] + a[i]->coord_0[l];
				}
				if (x[l] > wall[l]) {
					x[l] = x[l] - wall[l];
					a[i]->coord_0[l] = wall[l] - a[i]->coord_0[l];
				}	
			}
			k++;

			if (tag_photo ){
				foutput << k << " " << a[i]->flag  << " " << a[i]->print(wall) << endl;
			}
			for (int l = 0; l < 3; l++) {
				a[i]->coord_p[l] = a[i]->coord_0[l];
				a[i]->coord_0[l] = x[l];
			}
			
			//if (time == MSD_STEP) for (int i = 0; i < N; i++) a[i]->set_initial();
			MSD = 0;
			//if (time > MSD_STEP) {
				for (int i = 0; i < N; i++) for (int l = 0; l < 3; l ++) VC += (a[i]->coord_0[l]-a[i]->coord_p[l])*a[i]->velocity_initial[l]/dtime;
				for (int i = 0; i < N; i++) for (int l = 0; l < 3; l ++) MSD+= pow(a[i]->coord_0[l]-a[i]->coord_initial[l],2.0);
				MSD = MSD/N/6/time/dtime;
			//}
			
		}
		return 0;
	}

	int energy_collect( ofstream &fenergy, const int time, const bool fluct_flag)
		//collect Energy quantity
		{
		P = 0.0;
		T = 0.0;
		for (int i = 0; i < N; i++){
			T += mass*(pow(a[i]->coord_0[0]-a[i]->coord_p[0],2.0) + pow(a[i]->coord_0[1]-a[i]->coord_p[1],2.0) + pow(a[i]->coord_0[2]-a[i]->coord_p[2],2.0))/pow(dtime,2.0)/2; // *
			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_0[l], a[j]->coord_0[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0),0.5); // distance between atoms                             TODO not true
				//if (r  < r_kr){
					P += 4*epsilon*(pow(sigma/r,12) - pow(sigma/r, 6)); // F/r	
				//}
			}
		}
		P = P/N/1.23e-23;
		T= T/N/1.23e-23;
		E = (P + T);
		
		T = 2.0/3*T;
	
		fenergy << time << "\t" << T << "\t" << E << "\t" << MSD << "\t" << VC/N/3 << "\t" << P << endl;
		if (fluct_flag) {
			if ( (E_max == 0) || (E_max < E) ) E_max = E;
			if ( (E_min == 0) || (E_min > E) ) E_min = E;
		}
		return 0;
	}

	int state(double* g, double* v, double* T_get, double* P, double* MSD, double* VC_get){
		//acceleration collect
		double max_len = wall[0]/2;// pow( pow(wall[0],2.0) + pow(wall[0],2.0) + pow(wall[0],2.0),0.5);
		double velocity;
		double max_velocity = pow(3*1.23e-23*T/mass,0.5)*2;
		*T_get = T;
		*P = 0;
		for (int i = 0; i < N; i++) {
			velocity = pow(pow(a[i]->coord_0[0]-a[i]->coord_p[0],2.0) + pow(a[i]->coord_0[1]-a[i]->coord_p[1],2.0) + pow(a[i]->coord_0[2]-a[i]->coord_p[2],2.0), 0.5)/dtime;
			v[int(min(velocity, max_velocity)/max_velocity*V_STEPS)] += 1;

			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_0[l], a[j]->coord_0[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0),0.5); // distance between atoms                             TODO not true
				if (r < max_len) g[int(r/max_len*GIST_STEPS)] += 1;
				*P += 24.0*epsilon*(2*pow(sigma/r,12.0) - pow(sigma/r, 6.0)); // F/r
			}
		}


		*P = (N*1.23e-23*T - 1.0/3*(*P)/N)/pow(wall[0],3);


		for (int i = 0; i < GIST_STEPS; i++) g[i] = 16.0*g[i]/(4/3*3.14*(pow(i+1,3)-pow(i,3)))*pow(GIST_STEPS,3)/N/N;
		for (int i = 0; i < V_STEPS; i++) v[i] = v[i]/N;

		//MSD
		*MSD = 0;
		for (int i = 0; i < N; i++) *MSD += pow(a[i]->coord_0[0]-a[i]->coord_initial[0],2.0) + pow(a[i]->coord_0[1]-a[i]->coord_initial[1],2.0) + pow(a[i]->coord_0[2]-a[i]->coord_initial[2],2.0);
		*MSD = *MSD/N/3;
		*VC_get = VC/N/3;
		return 0;
	}

	int velocity_reverse() {
		double temp;
		for (int i = 0; i < N; i++){
			for (int l = 0; l < 3; l++){
				temp = a[i]->coord_p[l];
				a[i]->coord_p[l] = a[i]->coord_0[l];
				a[i]->coord_0[l] = temp;
			}
		}
		return 0;
	}
		

	
};

int main(int argc, char* argv[]){
	ofstream foutput, fenergy;

	//create g(r)
	double g_r[GIST_STEPS];
	for (int i = 0; i < GIST_STEPS; i++) g_r[i] = 0;
	//
	//create n(v)
	double v[V_STEPS];
	for (int i = 0; i < V_STEPS; i++) v[i] = 0;
	//
	double P,T, MSD, VC;

	zero_cond initial(argc, argv); // collect zero conditions
	
	//string name_file = "term/output_photo("+to_string(int(initial.dtime*1e16))+","+to_string(int(initial.r_kr))+","+to_string(int(initial.wall*1e10))+").txt" ;

	if (initial.max_photo) foutput.open("term/output_photo.txt");


	// dependense every_thing by time

	//name_file = "term/date_file("+to_string(int(initial.dtime*1e16))+","+to_string(int(initial.r_kr))+","+to_string(int(initial.wall*1e10))+").txt" ;

	if (initial.max_energy) fenergy.open("term/date_file.txt");
	fenergy << "#  time, fc \t E, J \t T, K \t MSD (Einstein) \t MSD (GREBO-CUBO) \t P, J" << endl; 
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
		if ( initial.max_photo && (!(step % initial.photo_step)) && (step <= initial.max_photo)){
			foutput << "ITEM: TIMESTEP\n" << step << "\nITEM: NUMBER OF ATOMS\n" << N << "\nITEM: BOX BOUNDS pp pp pp\n" << 0e0  << " " << wall[0]*1e10 << "\n" << 0e0  << " "  << wall[1]*1e10 << "\n" << 0e0  << " "  << wall[2]*1e10 << "\nITEM: ATOMS id type xs ys zs" << endl;
			argon->displace(1, foutput, step);
		} 
		else argon->displace(0, foutput, step);
		if ( initial.max_energy && (!((step+1) % initial.energy_step)) && (step <= initial.max_energy)) {if (max_step-step > 10) argon->energy_collect(fenergy, step, 0);}
		else argon->energy_collect(fenergy, step, 1);

		if (step == max_step) argon->state(g_r, v, &T, &P, &MSD, &VC);
		//if (step == max_step/2) argon->velocity_reverse();

		if (100.0*step/max_step >= percentage){
			cout << percentage << endl;
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
	
	/*
	//write fluctuation of energy
	ofstream ffluct; //fluct of energy
	if (initial.max_energy) ffluct.open("term/energy_fluct.txt", ios::app);
	ffluct << initial.dtime << "\t" << E_max-E_min << endl;
	ffluct.close();
*/
	//write state
	fenergy << "\n\n# P = " << P << "\t T = " << T << endl;
	fenergy << "# PV/NkT = " << P*pow(wall[0],3)/(N*1.23e-23*T) << endl;
	fenergy << "# MSD = " << MSD << "\t VC = " << VC << endl;
	//write g(r)
	fenergy << "# r \t g(r)" << endl;
	for (int i = 0; i < GIST_STEPS; i++) fenergy <<  i*wall[0]/2/GIST_STEPS << "\t" << g_r[i] << endl;
    //write n(v)
	fenergy << "\n\n# n \t v" << endl;
	double max_velocity = pow(3*1.23e-23*T/mass_Ar,0.5)*2;
	for (int i = 0; i < V_STEPS; i++) fenergy <<  i*max_velocity/V_STEPS << "\t" << v[i] << endl;
	fenergy.close();


	delete argon;
	cerr << "Everything ok!\n Finish!" << endl;
	return 0;
}