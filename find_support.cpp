//find_support.cpp
#include <iostream>
#include <math.h>
#include <cmath>
#include <cstdio>
#include <ctime>


#include "find_support.h"


void find_complicated(real_num Bperp, real_num Bbar,  real_num delta_a, 
		      real_num g, real_num &tt, real_num &c, real_num N,
		      int mem)
{
  
  const real_num pi = atan(1) * 4.0_r;
  real_num step = -0.0005;
  real_num dc =  0.000001_r;
  int stepsize = 100000;
  real_num temp = 0.0_r;
  real_num oldtemp;
  real_num Btot = sqrt(Bperp*Bperp + Bbar*Bbar);
  real_num epsilon;
  if(mem == 0)
    {
      epsilon = 0.000001_r;
    }
  else
    {
      epsilon = N + step;
    }
  
  for(int i=0; i<=stepsize; i++)
    {
      oldtemp = temp;
      c = i * dc + epsilon;
      real_num yN = Bperp / c / c;
      real_num x = g * Btot / c  + tt;
      
      temp = delta_a - c * yN * f_a1(x,yN) + c / sqrt(pi) * f_a2(yN) 
	- c * yN / sqrt(pi) * f_a3(x);


      if(temp * oldtemp < 0)
	{
	  return;	 
	}
      
    }

  return;
}



void find_simple(real_num B, real_num delta, real_num &c)
{
  const real_num pi = atan(1) * 4.0_r;
  real_num dc = 0.00005_r;
  real_num temp = 0.0_r;
  int stepsize = 1000;
  real_num epsilon = 0.0000001_r;
  
  for(int i=0; i<=stepsize; i++)
    {
      real_num oldtemp = temp;
      c = i * dc + epsilon;
      real_num y = B / c / c;
      
      temp = delta + c / sqrt(pi) * f_1(y);

      if(temp * oldtemp < 0)
	{
	  return;
	}

    }

  return;
}

void newton_simple(real_num B, real_num delta, real_num& c)
{
  const real_num pi = atan(1) * 4.0_r;
  const real_num epsilon = 1.0e-6_r;
  const real_num h = 0.0002_r;
  real_num temp = c;
  real_num oldtemp = 0.0_r;
  int count = 0;
  real_num y;
  real_num y_f;
  real_num y_b;
  real_num f;
  real_num f_f;
  real_num f_b;
  real_num f_N;

  while(fabs(oldtemp - temp) > epsilon)
    {
      oldtemp = temp;
      y = B / oldtemp / oldtemp;
      y_f = B / (oldtemp+h) / (oldtemp+h);
      y_b = B / (oldtemp-h) / (oldtemp-h);

      f =  delta + (oldtemp) / sqrt(pi) * f_1(y);
      f_f =  (oldtemp+h) / sqrt(pi) * f_1(y_f);
      f_b =  (oldtemp-h) / sqrt(pi) * f_1(y_b);

      f_N = (f_f - f_b) / h / 2.0_r;

      if(fabs(f_N) < 2.0e-25_r)
	{
	  std::cout << "Vanishing Jacobian found for m = 0 via Newton's method at N = " << oldtemp << ".\n";
	  std::cout << "\t \t \t f_N = " << f_N << ".\n";
	  exit(1);
	}
      
      temp = oldtemp - f/f_N;
      

      count++;
      if(count > 100)
	{
	  if(count%100 == 0)
	    {
	      std::cout << "\t \t \t m = 0 Newton's method stepped " << count << "times. \n";
	      std::cout << "\t \t \t Difference between steps: " << fabs(oldtemp - temp) << ".\n";
	    }	      
	}

      
    }
  
  c = temp;
  std::cout << "\t \t m = 0 root found via Newton's method after " << count + 1 << " number of steps.\n";
  return;

}


void iterate(real_num Bperp, real_num Bbar,  real_num delta_a, 
	     real_num delta_m, real_num g, real_num &N, real_num Nc,
	     real_num &theta, int iteration,int mem)
{
  
  const real_num tolerance = 1.0e-7_r;
  real_num theta_old = tan(theta);
  real_num nt = N;
  real_num Btot = sqrt(Bperp*Bperp+Bbar*Bbar);
  int count = 0;
  
  //for(int count = 0; count < iteration; count++)
  //{
 loop:
      real_num no = nt;
      
      real_num yN = Bperp / no / no;
      real_num x = g * Btot / no + theta_old;
      real_num thetatemp = no * yN / delta_m * f_f(x,yN); 
      
      find_complicated(Bperp, Bbar, delta_a, g, theta_old, nt, Nc,mem);

      theta_old = thetatemp;

      thetatemp = nt * yN / delta_m * f_f(x,yN);
      
      std::cout << "\t \t " << count +1 << " number of iteration(s) through theta performed.\n";

      if(fabs(nt - no) >= tolerance)
      	{
      	  std::cout << "\t \t Differnece between N values for successive steps through theta: " << fabs(nt - no) << ".\n";
      	  count++;
      	  goto loop;
      	}
      
      N = nt;
      theta = atan(thetatemp);
      
      //}
  
  return;
}


void boyden_comp(real_num Bperp, real_num Bbar,  real_num delta_a, 
		 real_num delta_m, real_num g, real_num &N, real_num &theta)
{
  const real_num pi = atan(1) * 4.0_r;
  const real_num epsilon = 5.0e-7_r;
  real_num t = 1.0_r;
  real_num N_temp = N;
  real_num theta_temp = theta;
  real_num Btot = sqrt(Bperp*Bperp + Bbar*Bbar);
  real_num invJac [2][2] = {{1,0},{0,1}};
  int count = 0;
  real_num O_old[2];
  real_num f_old[2];
  real_num p[2];
  real_num N_new;
  real_num theta_new;
  real_num O[2];
  real_num f_new[2];
  real_num old_norm;
  real_num new_norm;
  real_num difference_norm;
  real_num y[2];
  real_num pBy = 0.0_r;
  real_num By[2];
  real_num pB[2];

 loop:
  O_old[0] = Bperp / N_temp / N_temp;
  O_old[1] = g * Btot / N_temp + theta_temp;
  f_old[0] = delta_a - N_temp * O_old[0] * f_a1(O_old[1],O_old[0]) + N_temp / sqrt(pi) * f_a2(O_old[0]) - N_temp * O_old[0] / sqrt(pi) * f_a3(O_old[1]);
  f_old[1] = N_temp * O_old[0] / delta_m * f_f(O_old[1],O_old[0]);
  old_norm = sqrt(f_old[0]*f_old[0] + f_old[1]*f_old[1]);
  
  p[2] = {};
  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
	{
	  p[i] += -(invJac[i][j]*f_old[j]);
	}
    }
  
 restart:
  N_new = N_temp + t*p[0];
  theta_new = theta_temp + t*p[1];
  O[0] = Bperp / N_new / N_new;
  O[1] = g * Btot / N_new + theta_new;
  f_new[0] = delta_a - N_new * f_a1(O[1],O[0]) + N_new / sqrt(pi) * f_a2(O[0]) 
    - N_new * O[0] / sqrt(pi) * f_a3(O[1]);
  f_new[1] = N_new * O[0] / delta_m * f_f(O[1],O[0]);
  new_norm = sqrt(f_new[0]*f_new[0] + f_new[1]*f_new[1]);
  
  if(old_norm <= new_norm)
    {
      t -= 0.05_r;
      goto restart;
    }
  
  difference_norm = sqrt((f_old[0]-f_new[0])*(f_old[0]-f_new[0]) + (f_old[1]-f_new[1])*(f_old[1] - f_new[1]));
  if(difference_norm > epsilon)
    {
      y[2] = {};
      pBy = 0.0_r;
      By[2] = {};
      pB[2] = {};
      for(int n=0;n<2;n++)
	{
	  y[n] = f_new[n] - f_old[n];
	  for(int m=0;m<2;m++)
	    {
	      pBy += p[m]*invJac[m][n]*y[n];
	      By[m] += invJac[m][n]*y[n];
	      pB[m] += p[n]*invJac[n][m];
	    }
	}
      if(fabs(pBy) < 2.0e-25_r)
	{
	  std::cout << "Division by zero encountered for nonzero m via Boyden's method at (N,m) = (" << N_new << "," << theta_new << ").\n";
	  exit(1);	     
	}
      for(int i=0;i<2;i++)
	{
	  for(int j=0;j<2;j++)
	    {
	      invJac[i][j] -= ((By[i] - p[i])*pB[j])/pBy;
	    }
	}
      count++;
      // if(N_new < 0 || theta_new < 0)
      // 	{
      // 	  N_temp = N_new;
      // 	  theta_temp = theta_new;
      // 	  goto loop;
      // 	}

      N_temp = N_new;
      theta = theta_new;
      std::cout << "\t \t Performed  " << count +1 << " steps through Boyden's method.\n";
      std::cout << "\t \t \t N = " << N_new << ", theta = " << theta_new << ".\n";
      std::cout << "\t \t \t Difference of roots between steps: (" << fabs(N_new - N_temp) << "," << fabs(theta_new - theta_temp) << ").\n";
      std::cout << "\t \t \t Norm of difference between steps = " << difference_norm << ".\n";
      goto loop;     
    }
  else
    {
      N = N_new;
      theta = theta_new;
      std::cout << "\t \t Nonzero m roots found via Boyden's method after " << count + 1 << " number of steps.\n";
      std::cout << "\t \t \t N = " << N_new << ", theta = " << theta_new << ".\n";
      return;
    }
}
