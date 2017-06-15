#include<iostream>
#include<cmath>
#include<time.h>
#include<fstream>
#include<stdlib.h>
#include<vector>
#include"ewald.h"
#include<mpi.h>
#include "angular.h"

using namespace std;

const double Pi = 3.141593;
const double Unit_Coulomb = 1389.55;
const double Unit_eV = 96.49930;
const double R        = 0.0083144598; // Gas constant
const double Temp     = 2000;         //Temperature of the system
const double Pr       = 0.6023 * ;      //Pressure of the system
const int equ_iter    = 2500000;
const int iter        = 5000000;
const int inter       = 200;
const double kappa    = 0.20;
const double v_c      = 5;

/*
// generates a random  integer between [min,max)
int Rand_INT(int min,int max) //[min,max)
{
  return rand()%(max-min)+min;
}

// generates the random number between 0 and 1
double Rand_DOUBLE() //(0,1)
{
  return (double) rand() / (double) RAND_MAX;
}

// finds the minimum (a,b)
double min(double a,double b)
{
  if(a>b)
  {
    return b;
  }
  else
  {
    return a;
  }
} */


// Maximum of 3 elements
template <class T>
T max3(T a, T b, T c)
{
 return (a < c && b < c) ? c : (a < b) ? b : a;
}

// Short-range interaction - interatomic potentials were taken from Swamy and JD Gale paper
double ShortRange(vector<double> const &r, double Lx, double Ly, double Lz, double q1, double q2, int size, int rank)
{
  double u = 0, rij  = 0;  
  // MPI variables
  int from = 0, to = 0;
  double Ai = 0, Aj = 0, Bi = 0, Bj = 0, Ci = 0, Cj = 0;
  from = rank * r.size() / size ;
  to = (rank + 1) * r.size() / size ;
  if(rank + 1 == size)
  {
    to = to - 4;
  }

  for(unsigned int i = from;i < to ;i += 4)
    {
      for(unsigned int j = i+4;j < r.size();j += 4)
        {
          rij = sqrt( pow((r[i+0]-r[j+0]) - Lx * round((r[i+0]-r[j+0])/Lx),2) 
                    + pow((r[i+1]-r[j+1]) - Ly * round((r[i+1]-r[j+1])/Ly),2) 
                    + pow((r[i+2]-r[j+2]) - Lz * round((r[i+2]-r[j+2])/Lz),2) );

          // Potential truncation
          if(rij < 10)
            {
              if(r[i+3] == r[j+3] && r[i+3] == q1)//Ti - Ti
                {
                  Ai = 1.1823; Aj = 1.1823; Bi = 0.077; Bj = 0.077; Ci = 22.5; Cj = 22.5;
                  u +=   4.184 * (Bi+Bj) * exp((Ai+Aj-rij)/(Bi+Bj)) - (Ci*Cj/pow(rij,6));
                }
              if(r[i+3] == r[j+3] && r[i+3] == q2)//O - O
                {
                  Ai = 1.6339; Aj = 1.6339; Bi = 0.117; Bj = 0.117; Ci = 54.0; Cj = 54.0;
                  u +=   4.184 * (Bi+Bj) * exp((Ai+Aj-rij)/(Bi+Bj)) - (Ci*Cj/pow(rij,6));
                }
              if(r[i+3] != r[j+3]) // Ti - O
                { 
                  Ai = 1.1823; Aj = 1.6339; Bi = 0.077; Bj = 0.117; Ci = 22.5; Cj = 54.0;
                  u +=   4.184 * (Bi+Bj) * exp((Ai+Aj-rij)/(Bi+Bj)) - (Ci*Cj/pow(rij,6));
                }
            }
          /*if(rij < 10)
            {
              if(r[i+3] == r[j+3] && r[i+3] == q1)//Ti - Ti
                {
                  u +=   31120.2 * Unit_eV * exp(-(rij/0.154)) - 5.25 * Unit_eV / pow(rij,6);
                }
              if(r[i+3] == r[j+3] && r[i+3] == q2)//O - O
                {
                  u +=   11782.76 * Unit_eV * exp(-(rij/0.234)) - 30.22 * Unit_eV / pow(rij,6);
                }
              if(r[i+3] != r[j+3]) // Ti - O
                { 
                  u +=   16957.53 * Unit_eV * exp(-(rij/0.194)) - 12.59 * Unit_eV / pow(rij,6);
                }
            }*/
        }
    }
  return u;
}



// Main implementation
int main(int argc, char** argv)
{

  int size, rank, root = 0; 
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Variables
  double L = 0, V = 0, V_new = 0;
  unsigned int nx = 0, ny = 0, nz = 0, mean = 0 ;
  double r0 = 0;
  unsigned int N = 0, N_ = 0;
  unsigned int per = 0;
  double q1 = 0 , q2 = 0;
  vector<double> r;
  vector<double> r_new;

  // Properties variables
  // 1. Angular distribution variables
  int a_size = 360;
  double** otio = new double*[2];
  double** tioti = new double*[2];
  for(int i = 0;i < 2; i++)
    {
      otio[i] = new double[a_size];
      tioti[i] = new double[a_size];
    }    
  int j = 0;
  for(double i = 1;i < 180.5;i = i + 0.5)
    {
      otio[0][j] = i;
      otio[1][j] = 0;
      tioti[0][j] = i;
      tioti[1][j] = 0;
      j = j + 1;
    }

  //2. Radial Distribution Function(RDF)
  unsigned int RDF_size = 700;
  double RDF_h = 0.1;
  double* rad = new double[RDF_size];
  double* RDF_all = new double[RDF_size];
  double* RDF_Ti = new double[RDF_size];
  double* RDF_O  = new double[RDF_size];
  double* RDF_TiO =  new double[RDF_size];
  for(unsigned int i = 1;i < RDF_size ;i++)
    {
      rad[i]          = i * 0.01;
      RDF_all[i]    = 0;
      RDF_Ti[i]     = 0;
      RDF_O [i]     = 0;
      RDF_TiO[i]    = 0;
    }
  ofstream outfile1;
  ofstream outfile2;
  ofstream outfile3;

  // Processor 0 - input/output processing
  if(rank == root) 
    {
      double Lx = 0, Ly = 0, Lz = 0;
      ifstream infile(argv[1]);
      infile >> nx >> ny >> nz ;
      infile >> r0 ;
      infile >> N ;
      infile >> per ;
      per *= nx * ny * nz;
      infile >> q1 >> q2 ;
      infile >> Lx >> Ly >> Lz;
      N_ = N * 4 * nx * ny * nz ;
      r.resize(N*4*nx*ny*nz);
      r_new.resize(N*4*nx*ny*nz);

      for(unsigned int i = 0; i < (N*4); i += 4)
        {
          infile >> r[i+0] >> r[i+1] >> r[i+2] >> r[i+3];
          r[i+0] *= Lx;  r[i+1] *= Ly; r[i+2] *= Lz; 
        }
      infile.close();
      infile.clear();

      unsigned int i = 0;
      for(unsigned int x = 0; x < nx; ++x)
        {
          for(unsigned int y = 0; y < ny; ++y)
            {
              for(unsigned int z = 0; z < nz; ++z)
                {
                  for(unsigned int j = 0;j < (N*4); j += 4)
                    {
                      r[i+0] = Lx * x + r[j+0];
                      r[i+1] = Ly * y + r[j+1];
                      r[i+2] = Lz * z + r[j+2];
                      r[i+3] = r[j+3];
                      //cout << i << " " <<  r[i+0] << " " << r[i+1] << " " << r[i+2] <<" " << r[i+3] << endl;
                      i += 4;
                    }
                }
            }  
        }
      r_new = r;
      Lx *= nx;
      Ly *= ny;
      Lz *= nz; cout << Lx << " " << Ly << " " << Lz << endl;
      L = max3(Lx,Ly,Lz);
      V = pow(L,3);
      V_new = V;
      outfile1.open(argv[2]);
      outfile2.open(argv[3]);
      outfile3.open(argv[4]);    }

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&N_, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  r.resize(N_);
  r_new.resize(N_);
  MPI_Bcast(&q1, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&q2, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&V, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&L, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  MPI_Bcast(&r[0], r.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);

  // Ref equation 3.9 in (https://www.researchgate.net/publication/311138439_A_study_on_the_structure_and_properties_of_silica_glass_and_silica_nanoparticles_via_Monte_Carlo_simulations)
  double U_local = 0, U_global = 0, u = 0, u_new = 0;
  U_local = ( RealandReciprocalSpace(r, L, L, L, kappa, 1, size, rank) 
             + kappa * PointEnergy(r, size, rank) / sqrt(Pi) ) * Unit_Coulomb + ShortRange(r, L, L, L, q1, q2, size, rank);
  MPI_Reduce(&U_local, &u, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
  u += (Unit_Coulomb *  2 * Pi * pow(Dipole(r),2) / (3 * V)) ;

  srand(time(NULL));
  for(int k = 0;k < iter; k++)
    {
      if(rank == root)
        {
          //1. randomly select an atom
          unsigned int P = Rand_INT(0, r.size()/4) * 4;

          //2. move the selected atom randomly
          r_new[P+0] += 0.1*(2*Rand_DOUBLE()-1); 
          r_new[P+1] += 0.1*(2*Rand_DOUBLE()-1);
          r_new[P+2] += 0.1*(2*Rand_DOUBLE()-1);
          
          P = Rand_INT(0, r.size()/4) * 4;
          r_new[P+0] += 0.1*(2*Rand_DOUBLE()-1); 
          r_new[P+1] += 0.1*(2*Rand_DOUBLE()-1);
          r_new[P+2] += 0.1*(2*Rand_DOUBLE()-1);
          
          //3. Volume change
          V_new += v_c * (2*Rand_DOUBLE()-1);
          L = pow(V_new,0.333333);
          double ratio = pow(V_new,0.3333333)/pow(V,0.3333333);
          for(unsigned int i = 0 ; i < r.size();i = i + 4)
            {
              r_new[i+0] = r_new[i+0]*ratio;
              r_new[i+1] = r_new[i+1]*ratio;
              r_new[i+2] = r_new[i+2]*ratio;
            }
          U_global = 0;
        }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Bcast(&L, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&r_new[0], r_new.size(), MPI_DOUBLE, root, MPI_COMM_WORLD);
      MPI_Bcast(&V_new, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);

      //4. Calculate the potential energy of the configurations
      U_global = 0;
      U_local = ( RealandReciprocalSpace(r_new, L, L, L, kappa, 1, size, rank) 
             + kappa * PointEnergy(r_new, size, rank) / sqrt(Pi) ) * Unit_Coulomb + ShortRange(r_new, L, L, L, q1, q2, size, rank);

      MPI_Reduce(&U_local, &U_global, 1, MPI_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
     if (rank == root)
       {
         u_new = U_global + (Unit_Coulomb *  2.0 * Pi * pow(Dipole(r_new),2) / (3.0 * V_new));
        
        //5. Minimum Metropolis condition
         double weight = exp((r.size()/4.0)*log(V_new/V)- (1/(R*Temp))*((u_new-u)+Pr*(V_new-V)));
         if(min(1,weight) >= Rand_DOUBLE())
           {
             r = r_new;
             V = V_new;
             u = u_new;
             L = pow(V,0.33333);
           }
         else
           {
             r_new = r;
             V_new = V;
             L = pow(V,0.33333);
           }
         
         // Averages of the properties
         if(k%1000 == 1)
           {
            outfile1 << k <<" "<< ((r.size()/4) * 44.3 /V) << "  " <<  u/(96.49930*(r.size()/4)/3) <<   "\n";
            cout << k <<" "<< ((r.size()/4) * 44.3 /V) << "  " <<  u/(96.49930*(r.size()/4)/3) <<  "\n";
           }      
         if(k > equ_iter && k%inter == 1)
           {
             mean += 1;
             
             // 1. Angular distribution function
             angular_distribution_OTiO(L, L, L, otio, a_size, r, q1, q2);
             angular_distribution_TiOTi(L, L, L, tioti, a_size, r, q1, q2);

             // 2. Radial Distribution Function
             for(unsigned int i = 1;i < RDF_size-1;++i)
               {
                 unsigned int N_All = 0,N_Ti = 0,N_O = 0,N_TiO = 0;
                 double rij = 0;
                 for(unsigned int j = 0;j < r.size();j = j + 4)
                   {
                     for(unsigned int k_ = 0;k_ < r.size() ;k_ = k_ + 4)
                       {
                         rij = sqrt( pow((r[k_+0]-r[j+0]) - L * round((r[k_+0]-r[j+0])/L),2) 
                                   + pow((r[k_+1]-r[j+1]) - L * round((r[k_+1]-r[j+1])/L),2) 
                                   + pow((r[k_+2]-r[j+2]) - L * round((r[k_+2]-r[j+2])/L),2) );

                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && j != k_ )
                           {
                             N_All += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] == r[j+3] && r[k_+3] == q2 && j != k_)
                           {
                             N_O += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] == r[j+3] && r[k_+3] == q1 && j != k_ )
                           {
                             N_Ti += 1;
                           }
                         if(rij <= rad[i] + (RDF_h/2) && rij > rad[i] - (RDF_h/2) && r[k_+3] != r[j+3] && r[k_+3] == q1 && j != k_)
                           {
                             N_TiO += 1;
                           }
                       }
                   }

                 RDF_all[i] += N_All  /( ((r.size()/4) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_Ti[i]  += N_Ti  / ( (((r.size()/4)*1/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_O [i]  += N_O   / ( (((r.size()/4)*2/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
                 RDF_TiO[i] += N_TiO / ( (((r.size()/4)*1/3) / V) * 4 * Pi * rad[i] * rad[i] *  RDF_h);
               }
           }
      }
    }

  if(rank == root)
    {

      ofstream position("output/data/position.txt");
      for(unsigned int i = 0; i < r.size();i += 4)
        {
          position << r[i+0] << "  " << r[i+1] << "  " << r[i+2] << "  " << r[i+3] << "\n";    
        }
      position.close();

      ofstream length("output/data/length.txt");
      length << L << "   " << L << "   " << L << "\n Volume = " << V << "\n";    
      length.close();

      for(unsigned int i = 0; i < RDF_size; ++i)
        {
          outfile2 << rad[i] << " " << RDF_all[i] / (mean * r.size()/4) << " " << RDF_Ti[i] / (mean * (r.size()/4)/3) 
                   << " " << RDF_O [i] / (mean * (r.size()/4) * 2/3) << " " << RDF_TiO[i] / (mean * (r.size()/4) * 2/3) << "\n";
        }

      // Angular distribution
      int total1 = 0, total2 = 0;
      for(unsigned int i = 0;i < a_size; ++i)
        {
          total1 += tioti[1][i];
          total2 += otio[1][i];
        } 
      for(unsigned int i = 0;i < a_size; ++i)
        {
          tioti[1][i] = tioti[1][i]/total1;
          otio[1][i]  =  otio[1][i]/total2;
          outfile3 << tioti[0][i] << "  " << tioti[1][i] << " " << otio[1][i] << "\n";
        } 
    }

  MPI_Finalize();  
  return 0;
}
