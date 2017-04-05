#ifndef __ATOM__
#define __ATOM__

#include <iostream>
#include <fstream>
#include <string>

using namespace std;


class Atom {
public:
	double coord_p[3], coord_0[3], coord_initial[3]; // dynamic quantity in m
	double r_coord_p[3], r_coord_0[3]; //real coord

	double VC_coeff = 0;
	
	double f[3], velocity_initial[3], mass, velocity[3];
	int flag;
	Atom(const double* coord_of, const double* displace, const double mass_of);

	int set_initial();
	
	string print(const double* wall);

};


class Sys_atom
{	
public:
	double E, T, P, r;
	int N;
	Atom **a; // = new Atom* [N];	
	double sigma;
	double epsilon; 
	double mass;
	double dtime;
	double wall[3];
	double f;
	int k = 0;
	double force;
	double x[3];
	double x_rand[3];
	double VC = 0; //velocity_correlation
	double MSD = 0;
	double D_msd = 0;
	double D_gr = 0;
	double VC_coeff;
	double r_kr; //max rad of interraction

	Sys_atom(const double* wall_get, const double sigma_get, const double epsilon_get, const double mass_get, const double dtime_get, const double r_kr_get , const int N_get);
	~Sys_atom();

	double module(const double a, const double b);

	int displace(const int tag_photo, ofstream &foutput, const int time);

	int MSD_collect(const int time);

	int energy_collect( ofstream &fenergy, const int time, const bool fluct_flag);

	int state(double* g, double* v, double* T_get, double* P, double* MSD, double* VC_get);

	int velocity_reverse();
		
	double lj(double rast);

	int print_lj(ofstream &lj_file);
};

#endif // __ATOM__