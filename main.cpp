// ./a.out dtime max_photo photo_step max_energy energy_step wall r_kr
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <assert.h>
#include "lib/atom.h"

#define MAXNAME 10
#define GIST_STEPS 50
#define V_STEPS 20
#define MSD_STEP 400
using namespace std;

//constant for Lennard_Jones
double mass_Ar = 39.948 / 1000/6.022e23; // moll
double sigma_Ar = 3.405*1e-10; //A
double epsilon_Ar = 119.8*1.23e-23; // *k_Bol
//




class zero_cond{
public:
	//double N_part = 1000;
	double dtime;
	double wall, r_kr;
	int max_photo, photo_step, max_energy, energy_step, range;

	zero_cond(const int argc, char* argv[]){
		if (argc > 8) {
		dtime = atoi(argv[1])*1.0e-16;
		max_photo = atoi(argv[2]);
		photo_step = atoi(argv[3]);
		max_energy = atoi(argv[4]);
		energy_step = atoi(argv[5]);
		wall = atoi(argv[6])*1.0e-10;
		r_kr = atoi(argv[7])*sigma_Ar;
		r_kr = min(r_kr, wall/2);
		cout << "r_kr_new" << r_kr << endl;
		range = atoi(argv[8]);
		cout << "dtime = " << dtime << endl;
	} else {
		cerr << "Not enought arguments!" << endl;
	}
	}
};




int main(int argc, char* argv[]){
	ofstream foutput, fenergy;
	int N;
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
	N = initial.range;

	Sys_atom *argon = new Sys_atom(wall, sigma_Ar, epsilon_Ar, mass_Ar, initial.dtime, initial.r_kr, N);
	
	


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

		
		if (100.0*step/max_step >= percentage){
			cout << percentage << endl;
			percentage += 5;
		}


		step++;
	}

	cout << "Moving complete, finish energy file!" << endl;
	argon->state(g_r, v, &T, &P, &MSD, &VC);
/*
	FILE *energy_fluct;
	energy_fluct = fopen("energy_fluct.txt", "a+");
	fprintf(energy_fluct, "%d %lf\n", dtime, E_max-E_min);
	fclose(energy_fluct);

	ofstream lj_file;
	lj_file.open("term/lj.txt");
	argon->print_lj(lj_file);
	lj_file.close();
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
	cerr << "# P = " << P << "\t T = " << T << endl;
	cerr << "# PV/NkT = " << P*pow(wall[0],3)/(N*1.23e-23*T) << endl;
	cerr << "# MSD = " << MSD << "\t VC = " << VC << endl;
	//write g(r)
	fenergy << "\n\n# r \t g(r)" << endl;
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