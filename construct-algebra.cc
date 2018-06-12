#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <random>
#include <complex>
#include <iostream>
#include <fstream>
#include <chrono>

//from getopt example. Do I actually need all of these?
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "runtime-check.h"
#include "util.cc"
#include "chebyshev.cc"
#include "globals.h"

// Do we do runtime checks? Defined via makefile
//#define CHECK false
// Do we do slow checks (the ones that require creating a whole dense
// operator)?
// #define SLOW_CHECK true 
			
using namespace itensor;

#ifndef TESTING //test has its own main function
int main(int argc, char **argv)
{

  int L = 8; // System size
  int N = 8; // Chebyshev order. This is very small: be nice to use
	      // more like N = 100 or even N = 1000, but that's not
	      // great for testing purposes (note that this si
  
  int Maxm = 32; // Bond dimension cutoff. This is also pretty small,
		 // though we may be able to get away with it
  uint s = 0;
  double hz = 1.0;
  double cutoff = 1e-14;
  double Q = 1;
  bool profligate = false; //use memory-profligate strategy?
  std::string filename = "/tmp/conductivity";
  std::string model = "rfheis";

  int isweep_flag; //do we multiply with sweeps/fit?
  int nsweeps; //number of sweeps

  int option_index = 0;
  while (1) {
    static struct option long_options[] =
      {
	{"sweep", no_argument, &isweep_flag, 1},
	{"system-size",     required_argument,        0, 'L'},
	{"bond-dimension",  required_argument,        0, 'M'},
	{"chebyshev-order", required_argument,        0, 'N'},
	{"model",   required_argument, 0, 'm'},
	{"output",  required_argument, 0, 'o'},
	{"nsweeps",  required_argument, 0, 'w'},
	{0, 0, 0, 0}
      };

      char c = getopt_long (argc, argv, "L:M:N:h:f:e:Q:m:s:w:pd", long_options, &option_index);

      if (c == -1) break;
      switch (c) {
        case 0: if (long_options[option_index].flag != 0) break;
        case 'L':
	  L = std::stoi(optarg);
	  std::cout << "L = " << L << "\n";
	  break;
	  
        case 'M':
	  Maxm = std::stoi(optarg);
	  std::cout << "Maxm = " << Maxm << "\n";
          break;

        case 'N':
	  N = std::stoi(optarg);
	  std::cout << "N = " << N << "\n";
          break;

	case 'h':
	  hz = std::stod(optarg);
	  std::cout << "hz = " << hz << "\n";
	  break;
	  
	case 'e':
	  cutoff = pow(10, -std::stod(optarg));
	  std::cout << "cutoff = " << cutoff << "\n";
	  break;
	  
        case 's':
	  s = std::stoi(optarg);
	  std::cout << "s = " << s << "\n";
          break;
	  
        case 'w':
	  nsweeps = std::stoi(optarg);
	  std::cout << "nsweeps = " << nsweeps << "\n";
          break;
	  
        case 'f':
	  filename = optarg;
	  std::cout <<  filename << "\n";
          break;

	case 'm':
	  model = optarg;
	  std::cout << model << "\n";
	  break;

	  //todo disorder width
	  
	case 'Q':
	  Q = std::stod(optarg);
	  std::cout << "scale" << Q << "\n";
	  break;
	  
        case '?':
          /* getopt_long already printed an error message. */
          break;

        default:
          abort ();
        }
    }

  //==============================================================
  //Sanity check (are we doing runtime checks)

  bool sweep_flag = isweep_flag;
#ifdef CHECK
  std::cout << "runtime checks inoperative" << "\n";
  exit(1);

  if (8 < L ||  32 < N)
    {
      std::cerr << "runtime checks are a bad idea with large L, N; L, N ";
      std::cerr << L << " " << N << "\n";
    }
    
#endif
  
  //==============================================================
  //construct random field Heisenberg Hamiltonian, current operator
  
  std::complex<double> im(0,1); //1i
  
  //set up random number generation
  //cf https://isocpp.org/files/papers/n3551.pdf
  std::default_random_engine e{s};
  
  SiteSet sites = SpinHalf(L);
  writeToFile(filename+".sites",sites);
  
  std::vector<double> hzs;
  IQMPO H;
  double opnorm_bound;
  //should probably be case/switch
  if     (model == "rfheis") { std::tie(H, hzs, opnorm_bound) = rfheis(sites, hz, e); }
  else if(model == "xx")     { std::tie(H, hzs, opnorm_bound) = XX(sites, hz, e); }
  else if(model == "rpara")  { std::tie(H, hzs, opnorm_bound) = rpara(sites, hz, e); }
  else if(model == "2NJW")   { std::tie(H, hzs, opnorm_bound) = rf_2NJW(sites, hz, e); }
  else {std::cerr << "unknown model " << model ; exit(-1);}
  
  H *= Q/opnorm_bound;

  //write the disorder to file
  OPENE(filename + ".dis", disout_file);
  for(int b = 0; b < L; ++b) { disout_file << hzs[b] << " "; }
  disout_file << "\n";

  //==============================================================
  // compute mu
  t0 = std::chrono::high_resolution_clock::now();
  auto cheb = listandwrite_dangler(H, filename, N, Maxm, cutoff, sweep_flag, nsweeps, 32);
}

#endif //ifndef TESTING
