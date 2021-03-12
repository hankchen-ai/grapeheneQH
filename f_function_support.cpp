//f_function_support.cpp


/* Include standard libraries */
#include <iostream>
#include <sys/stat.h>
#include <set>
#include <cmath>


/* Include user-define libraries */
#include "f_function_support.h"


real_num f_a1(real_num x, real_num y)
{
  
  real_num sum = 0.0_r;
  int nmax = 500;
  
  for(int n=0; n<=nmax; n++)
    {
      real_num Y = 2.0_r * n * y;
      
      for(int sigma=-1; sigma<=1; sigma=sigma+2)
	{
	  real_num D = sqrt(Y) + sigma * x;
	  
	  sum += 1.0_r / sqrt(1.0_r + D * D) - 1.0_r / sqrt(1.0_r + Y);
	  
	}
      
    }
  
  return sum;
  
}


real_num f_a2(real_num y)
{
  real_num stepratio = 1.5_r;
  real_num ds = 0.001_r / stepratio;
  real_num epsilon = 0.00000001_r;
  real_num sum = 0.0_r;
  int nmax = (int) (10000 * stepratio);
  
  for(int n=0; n<=nmax; n++)
    {
      
      real_num s = n * ds + epsilon;
      real_num sy = s * y;
      
      sum += ds / pow(s,1.5_r) * (1.0_r - sy * exp(-s) / tanh(sy));
      
    }
  
  return sum;
  
}



real_num f_a3(real_num x)
{
  real_num ds = 0.002_r;
  real_num epsilon = 0.00000001_r;
  real_num sum = 0.0_r;
  int nmax = 5000;
  
  for(int n=0; n<=nmax; n++)
    {
      
      real_num s = n * ds + epsilon;
      real_num sx = s * x * x;
      
      sum += ds / sqrt(s) * exp(-s) * (1.0_r - exp(-sx));
      
    }
  
  return sum;
  
}



real_num f_f(real_num x, real_num y)
{

  real_num sum = 0.0_r;
  int nmax = 500;
  
  for(int n=0; n<=nmax; n++)
    {
      
      for(int sigma=-1; sigma<=1; sigma=sigma+2)
	{
	  
	  real_num D = sqrt(2.0_r * n * y) + sigma * x;
	  
	  sum += sigma * D / sqrt(1.0_r + D * D);
	  
	}

    }

  sum -= x / sqrt(1.0_r + x * x);
  
  return sum;
  
}



real_num f_1(real_num y)
{
  
  real_num ds = 0.001_r;
  real_num epsilon = 1.0e-8_r;
  real_num sum = 0.0_r;
  int nmax = 10000;
  
  for(int n=0; n<=nmax; n++)
    {
      
      real_num s = n * ds + epsilon;
      real_num sy = s * y;
      
      sum += ds / pow(s,1.5_r) * (1.0_r - sy * exp(-s) / tanh(sy));
      
    }
  
  return sum;
  
}
