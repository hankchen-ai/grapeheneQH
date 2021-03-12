//gap.cpp


#include <iostream>
#include <sys/stat.h>
#include <set>
#include <cmath>



#include "gap.h"

#include "readBulkData.h"
#include "printBulkData.h"

//#include <vector>


int main (int argc, char** argv)
{

    // typedef std::vector< std::vector< std::vector<double> > > 
    // 	vec_3D_dbl;

    // vec_3D_dbl data_tensor;
    // std::vector<int> i_vec, j_vec;

    // readBulkData("Ferromagnetic.dat", data_tensor, i_vec, j_vec);
    // printBulkData(data_tensor, i_vec, j_vec);

    // exit(1);
    
  int nmax = 220;
  int step = 2;
  
  intVector1D n_array;
  for(int array_index = 1;array_index <= nmax;array_index++)
    {
      if(array_index%step==0)
	{
	  n_array.push_back(array_index);
	}
    }

  intVector1D i_array = {7,45};
  intVector1D j_array = {1};

  int n_size = n_array.size();
  int i_size = i_array.size();
  int j_size = j_array.size();

  real_num N_temp;
  realVector3D Gamma(n_size), ferr(n_size), antif(n_size);

  const real_num g = 0.014_r;
  const real_num d_delta_a = 0.005_r;
  const real_num d_delta_m = 1.0_r;

  //Set number of iterations through theta:
  int iteration = 1;

  //Toggle N history:
  int mem = 0;
  
  real_num Bbar;
  realVector1D Bperp = {0.0001_r};
  //{0.000014_r,0.00003_r,0.00005_r,0.00007_r,0.00014_r,0.0003_r,0.0005_r,0.0007_r};
  real_num dB = 0.00002_r;
  
  // int p = 0;
  int B_size = Bperp.size();
  
  for(int p=0;p<B_size;p++)
    {
      for(int n_index = 0; n_index < n_size; n_index++)
	{

	  Gamma[n_index].resize(i_size);
	  antif[n_index].resize(i_size);
	  ferr[n_index].resize(i_size);
      
	  int n = n_array[n_index];
	  Bbar = n * dB;

	  std::cout << "n = " << n << ". \n";
      
	  for(int i_index=0; i_index<i_size; i_index++)
	    {
	  	      
	      Gamma[n_index][i_index].resize(j_size);
	      ferr[n_index][i_index].resize(j_size);
	      antif[n_index][i_index].resize(j_size);
	  
	      int i = i_array[i_index];
	      real_num delta_a = i * d_delta_a;
	  

	      find_simple(Bperp[p], delta_a, N_temp);

	      std::cout << " i = " << i << ".\n";

	      real_num N = N_temp;
	  
	      for(int j_index=0;j_index<j_size; j_index++)
		{	      
		  int  j = j_array[j_index];
		  real_num delta_m = j * d_delta_m;
		  real_num theta = 0.0_r;

		  std::cout << "  j =  "<< j << ". \n";

		  iterate(Bperp[p],Bbar,delta_a, delta_m, g, N_temp, N, theta, iteration,mem);

		  real_num m = N_temp * tan(theta);

		  std::cout << "\t \t \t N =  " << N_temp << ", m = " << m << ".\n";

		  real_num gB = g * sqrt(Bbar * Bbar + Bperp[p] * Bperp[p]);
	      
		  Gamma[n_index][i_index][j_index] = sqrt(N_temp*N_temp + (m + gB)*(m + gB));
		  ferr[n_index][i_index][j_index] = m;
		  antif[n_index][i_index][j_index] = N_temp;

		}
	    }
      
	}
  
      std::string magg = std::to_string(Bperp[p] * 10000).substr(0,3);
      
      test_writing_function(i_array, j_array, Gamma,  ferr,  antif,"Gamma" + magg + ".dat","Ferro" + magg + ".dat","AntiFerro" + magg + ".dat");
    }
  return 0;
  
} // End of main function
