#ifndef __ATOM_CPP__
#define __ATOM_CPP__

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <assert.h>
#include "atom.h"

#define MAXNAME 10
#define GIST_STEPS 50
#define V_STEPS 20
#define MSD_STEP 400
using namespace std;


//implication for class Atom

Atom::Atom(const double* coord_of, const double* displace, const double mass_of){
			for( int l = 0; l < 3; l++) {
				coord_0[l] = coord_of[l]; 
				coord_p[l] = coord_0[l]-displace[l];
				r_coord_0[l] = coord_0[l]; 
				r_coord_p[l] = coord_p[l];
				flag = 1;
				set_initial();
			}
			mass = mass_of;
	}

int Atom::set_initial(){
			for( int l = 0; l < 3; l++) {
				coord_initial[l] = coord_0[l]; 
			}
			return 0;
	}
string Atom::print(const double* wall){
			return to_string(coord_0[0]/wall[0]) + " " + to_string(coord_0[1]/wall[1]) + " " + to_string(coord_0[2]/wall[2]);
	}

//implication for class Sys_atom

Sys_atom::Sys_atom(const double* wall_get, const double sigma_get, const double epsilon_get, const double mass_get, const double dtime_get, const double r_kr_get , const int N_get){
		for (int l = 0; l < 3; l++) {
			wall[l] = wall_get[l];
			cout << "wall: " << wall[l] << endl;
		}
		N = N_get;
		a = new Atom* [N];
		sigma = sigma_get;
		epsilon = epsilon_get;
		mass = mass_get;
		dtime = dtime_get;
		r_kr = r_kr_get;
		cout << "r_kr = " << r_kr << endl;
		//N = (int) aa*b*c;
	//	*a = new Atom* [N]; //create massive of atoms
		double v[3] = {0, 0, 0}; //initial velocity
		double steps[3];
		double number_of_atom = pow(N, 1.0/3);
		if (number_of_atom > int(number_of_atom)) number_of_atom = int(number_of_atom) + 1;
		else number_of_atom = int(number_of_atom);

		for (int l = 0; l < 3; l++) steps[l] = wall[l]/number_of_atom;
		
		k = 0; //counter of atom
		for (x[0] = steps[0]/2; x[0] < wall[0]/2; x[0] += steps[0]){
			for (x[1] = steps[1]/2; x[1] < wall[1]; x[1] += steps[1]){
				for (x[2] = steps[2]/2; x[2] < wall[2]; x[2] += steps[2]){

					//too mach for save sum impuls = 0
					for (int l = 0; l < 3; l++) v[l] = (1.0*rand()/RAND_MAX - 0.5)*steps[l]/200000;
					a[k] = new Atom(x,v, mass_get);
					if (k == 0) a[k]->flag = 3;
					k++;

					if (k >= N) break;
					for (int l = 0; l < 3; l++) {
						v[l] = -v[l];
						x_rand[l] = x[l]; 
					}
					x_rand[0] = wall[0] - x[0];
					a[k] = new Atom(x_rand,v, mass_get);
					k++;

					if (k >= N) break;
					//cout << k <<endl;
				}

					if (k >= N) break;
			}

					if (k >= N) break;
		}
		if (k == N){
			cout << "initial complete!" << endl;
		} else {
			cout << "initial fail! k = "  << k << endl;
			assert(0);
		}
	}
Sys_atom::~Sys_atom(){
		for (int k = 0; k < N; k++) delete[] a[k];
		delete[] a;
		cerr << "Memory clear!" << endl;
	}

double Sys_atom::module(const double a, const double b) {
		if  (abs(a-b) < wall[0]/2.0) {
			return a-b;	                                       // TODO make it for all wall
		} else { 
			if (a < b) return a-b+wall[0];
			else return a-b-wall[0];
		} 
	}

int Sys_atom::displace(const int tag_photo, ofstream &foutput, const int time){
		//acceleration collect
		for (int i = 0; i < N; i++)  {
			for (int l = 0; l < 3; l++) a[i]->f[l] = 0;
			//a[i]->flag = 1;
		}
		//a[0]->flag = 3;
		for (int i = 0; i < N-1; i++) {
			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_0[l], a[j]->coord_0[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0) ,0.5); // distance between atoms
				if (r  < r_kr){
					f = lj(r); // F/r
					//if (i == 0) a[j]->flag = 2;
					for (int l = 0; l < 3; l++) {
						force =  f*x[l]/r;						
						a[i]->f[l] += force;
						a[j]->f[l] -= force;
					}

				}
			
			}
		}
		//make displacement
		k = 0;
			
		for (int i = 0; i < N; i ++) {
			for (int l = 0; l < 3; l++) {
				x[l] = 2*a[i]->coord_0[l] - a[i]->coord_p[l]+a[i]->f[l]*pow(dtime,2.0)/a[i]->mass; //collect new coordinate

				a[i]->velocity[l]= (a[i]->r_coord_0[l]+x[l]-a[i]->coord_0[l] - a[i]->r_coord_p[l])/2;
				
				//set real coordinate (boundary cond)
				a[i]->r_coord_p[l] = a[i]->r_coord_0[l];
				a[i]->r_coord_0[l] = a[i]->r_coord_0[l]+x[l]-a[i]->coord_0[l];
				
				if (time == 1) {
					a[i]->velocity_initial[l] = a[i]->velocity[l];
					a[i]->VC_coeff += a[i]->velocity_initial[l]*a[i]->velocity_initial[l];
					if (l == 2) for(int p = 0; p < 3; p++) a[i]->velocity_initial[p] = a[i]->velocity_initial[p]/a[i]->VC_coeff;
				}
				
				// boundary condition
				if (x[l] < 0.0) {
					x[l] = wall[l] + x[l];
					a[i]->coord_0[l] = wall[l] + a[i]->coord_0[l];
				}
				if (x[l] > wall[l]) {
					x[l] = x[l] - wall[l];
					a[i]->coord_0[l] = a[i]->coord_0[l] - wall[l];
				}	
			}
			k++;

			if (tag_photo ){
				foutput << k << " " << a[i]->flag  << " " << a[i]->print(wall) << endl;
			}

			//update displacement
			for (int l = 0; l < 3; l++) {
				a[i]->coord_p[l] = a[i]->coord_0[l];
				a[i]->coord_0[l] = x[l];
			}
			
		}

		MSD_collect(time);

		return 0;
}

int Sys_atom::MSD_collect(const int time){
		MSD = 0;
		VC = 0;
		for (int i = 0; i < N; i++) for (int l = 0; l < 3; l ++) {
			VC += a[i]->velocity[l]*a[i]->velocity_initial[l];
			MSD+= pow(a[i]->r_coord_0[l]-a[i]->coord_initial[l],2.0);
		}

		MSD = MSD/N; //
		//D_msd = MSD/6/time;
		VC = VC/N;
		//D_gr += VC/3*dtime;

		return 0;
}

int Sys_atom::energy_collect( ofstream &fenergy, const int time, const bool fluct_flag)
		//collect Energy quantity
		{
		P = 0.0;
		T = 0.0;
		for (int i = 0; i < N; i++){
			T += a[i]->mass*(pow(a[i]->velocity[0],2.0) + pow(a[i]->velocity[1],2.0) + pow(a[i]->velocity[2],2.0))/pow(dtime,2.0)/2; // *
			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_p[l], a[j]->coord_p[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0),0.5); // distance between atoms                             TODO not true
				if (r  < r_kr){
					P += 4*epsilon*(pow(sigma/r,12) - pow(sigma/r, 6)); // 
				}
			}
		}
		P = P/N/1.23e-23;
		T = T/N/1.23e-23;
		E = (P + T);
		
		T = 2.0/3*T;
		
		fenergy << time << "\t" << T << "\t" << E << "\t" << MSD << "\t" << VC << "\t" << P << "\t" << "\t" << D_gr << "\t" << D_msd << endl;
		
		return 0;
}

int Sys_atom::state(double* g, double* v, double* T_get, double* P, double* MSD, double* VC_get){
		//acceleration collect
		double max_len = wall[0]/2;// pow( pow(wall[0],2.0) + pow(wall[0],2.0) + pow(wall[0],2.0),0.5);
		double velocity;
		double T_velocity = 3*1.23e-23*T/mass;
		*T_get = T;
		*P = 0;
		for (int i = 0; i < N; i++) {
			velocity = (pow(a[i]->velocity[0],2.0) + pow(a[i]->velocity[1],2.0) + pow(a[i]->velocity[2],2.0))/dtime/dtime/T_velocity;
			
			v[int(min(velocity, 4.0)/4*V_STEPS)] += 1;
			for (int j = i+1; j < N; j++) {
				for (int l = 0; l < 3; l++) x[l] = module(a[i]->coord_p[l], a[j]->coord_p[l]);
				r = pow( pow(x[0],2.0) + pow(x[1],2.0) + pow(x[2],2.0),0.5); // distance between atoms                             TODO not true
				if (r < max_len) g[int(r/max_len*GIST_STEPS)] += 1;
				*P += 24.0*epsilon*(2*pow(sigma/r,12.0) - pow(sigma/r, 6.0)); // F/r
			}
		}


		*P = (N*1.23e-23*T - 1.0/3*(*P)/N)/pow(wall[0],3);


		for (int i = 0; i < GIST_STEPS; i++) g[i] = 16.0*g[i]/(4/3*3.14*(pow(i+1,3)-pow(i,3)))*pow(GIST_STEPS,3)/N/N;
		for (int i = 0; i < V_STEPS; i++) v[i] = v[i]/N;

		//MSD
		*MSD = D_msd;
		*VC_get = D_gr;
		return 0;
}

int Sys_atom::velocity_reverse() {
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
		
double Sys_atom::lj(double rast){
		if (rast <= 0.8*r_kr) return 24.0*epsilon*(2*pow(sigma/rast,12.0) - pow(sigma/rast, 6.0))/rast;
		if (rast < r_kr) return (r_kr - rast)/0.2*24.0*epsilon*(2*pow(sigma/r_kr/0.8,12.0) - pow(sigma/r_kr/0.8, 6.0))/r_kr/0.8/r_kr;
		if (rast >= r_kr) return 0;
		return 5;
}

int Sys_atom::print_lj(ofstream &lj_file){
	for (double i = r_kr*0.5; i < r_kr*1.1; i+= r_kr*0.6/1000) lj_file << i/r_kr << "\t" << lj(i) << endl;
	return 0;
}



#endif // __ATOM_CPP__