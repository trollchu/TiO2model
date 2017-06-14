#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include "angular.h"
#include "ewald.h"
#include<vector>


double mindis(double dx,double dy,double dz,double a,double b,double c)
{
  return norm(dx - a * round(dx/a),dy - b * round(dy/b),dz - c * round(dz/c));
}

int Rand_INT(int min,int max) //[min,max)
{
  return rand()%(max-min)+min;
}

double Rand_DOUBLE() //(0,1)
{
  return (double) rand() / (double) RAND_MAX;
}

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
}


double angle_btwn_3points(vector<double> const &r,int i,int j1,int j2,double a,double b,double c)
{
  double x1,x2;
  double y1,y2;
  double z1,z2;

  x1 = r[j1+0] - r[i];
  x2 = r[i+0] - r[j2];
  y1 = r[j1+1] - r[i+1];
  y2 = r[i+1] - r[j2+1];
  z1 = r[j1+2] - r[i+2];
  z2 = r[i+2] - r[j2+2]; 

  x1 = x1 - a*round(x1/a);
  x2 = x2 - a*round(x2/a);
  y1 = y1 - b*round(y1/b);
  y2 = y2 - b*round(y2/b);
  z1 = z1 - c*round(z1/c);
  z2 = z2 - c*round(z2/c);

  double top =  (x1*x2+y1*y2+z1*z2) - 0.0001 ;
  double bot =  norm(x1,y1,z1) * norm(x2,y2,z2);
  if( top == bot )
    {
      return 180; 
    }
  else
    {
      return  180 - (acos(top/bot) * 57.296);
    }
}


void angular_distribution_TiOTi(double Lx, double Ly, double Lz, double **tioti, int a_size, vector<double> const &r, double q1, double q2)
{
     const int end  = r.size();
     int size_ = 6;
     int index[size_];

     for(int i = 0 ; i < end;i = i + 4)
      {
        if(r[i+3] == q2)
        {
          int incre = 0;
          for(int j = 0 ; j < end;j = j + 4)
          {
            if(i == j){}
            else if(r[j+3] == q1)
            {
              double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],Lx,Ly,Lz);
              if(rij < 2.5)
              {
                index[incre] = j;
                incre += 1;
              }
            }
          }

          for(int ii = 0; ii <  incre-1; ii++)
          {
            for(int jj = ii+1;jj < incre;jj++)
            {
              double angle  = angle_btwn_3points(r,i,index[ii],index[jj],Lx,Ly,Lz);
              for(int kk = 0;kk < a_size;kk++)
              {
                if(angle <= tioti[0][kk] + 0.25 && angle > tioti[0][kk] - 0.25)
                {
                  tioti[1][kk] += 1;
                }
              }
            }
          } 

        }
      }
}

void angular_distribution_OTiO( double Lx, double Ly, double Lz, double **otio, int a_size, vector<double> const &r, double q1, double q2)
{

     const int end  = r.size();
     int size_ = 15;
     int index[size_];

     for(int i = 0 ; i < end;i = i + 4)
      {
        if(r[i+3] == q1)
        {
          int incre = 0;
          for(int j = 0 ; j < end;j = j + 4)
          {
            if(i == j){}
            else if(r[j+3] == q2)
            {
              double rij = mindis(r[i]-r[j],r[i+1]-r[j+1],r[i+2]-r[j+2],Lx,Ly,Lz);
              if(rij < 2.5)
              {
                index[incre] = j;
                incre += 1;
              }
            }
          }

          for(int ii = 0; ii <  incre-1; ii++)
          {
            for(int jj = ii+1;jj < incre;jj++)
            {
              double angle1  = angle_btwn_3points(r,i,index[ii],index[jj],Lx,Ly,Lz);
              for(int kk = 0;kk < a_size;kk++)
              {
                if(angle1 < otio[0][kk] + 0.25 && angle1 > otio[0][kk] - 0.25)
                  {
                    otio[1][kk] += 1;
                  }
              }
            }
          } 

        }
      }
}
