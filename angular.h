#ifndef ANGULAR_H
#define ANGULAR_H

#include<vector>

using namespace std;


double mindis(double dx,double dy,double dz,double a,double b,double c);

int Rand_INT(int min,int max);

double Rand_DOUBLE();

double min(double a,double b);

double angle_btwn_3points(vector<double> const &r,int i,int j1,int j2,double a,double b,double c);

void angular_distribution_TiOTi(double Lx, double Ly, double Lz, double **tioti, int a_size, vector<double> const &r, double q1, double q2);

void angular_distribution_OTiO(double Lx, double Ly, double Lz, double **otio, int a_size, vector<double> const &r, double q1, double q2);

#endif
