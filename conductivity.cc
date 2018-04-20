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
  bool dangler = false;
  std::string filename = "/tmp/conductivity.txt";
  std::string model = "rfheis";


  //==============================================================
  //Parse Command Line Arguments
  while (1)
    {
      static struct option long_options[] =
        {
          /* These options don’t set a flag.
             We distinguish them by their indices. */
          {"system-size",     required_argument,        0, 'L'},
          {"bond-dimension",  required_argument,        0, 'M'},
          {"chebyshev-order", required_argument,        0, 'N'},
	  {"output-file",     required_argument,        0, 'f'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      char c = getopt_long (argc, argv, "L:M:N:h:f:e:Q:m:pd", long_options, &option_index);

      /* Detect the end of the options. */
      if (c == -1)
        break;

      switch (c)
        {
        case 0:
          /* If this option set a flag, do nothing else now. */
          if (long_options[option_index].flag != 0)
            break;

	  //alternate strategies
	case 'p':
	  profligate = true;
	  break;
	  
	case 'd':
	  dangler = true;
	  break;
	  
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
  //set some variables, mostly files/filenames


  std::string realmu_filename = filename + ".re";
  std::string imagmu_filename = filename + ".im";
  std::string disout_filename = filename + ".dis";
  std::string timing_filename = filename + ".tim";
  std::string chebbd_filename = filename + ".chM"; //record bond dimensions of Chebyshev polynomials
  
  // macro OPENE declares second arg and initializes with filehandle to
  // firstarg. Defined in util.cc

  OPENE(filename + ".re",  realmu_file);
  OPENE(filename + ".im",  imagmu_file);
  OPENE(filename + ".dis", disout_file);
  OPENE(filename + ".tim", timing_file);
  OPENE(filename + ".chs", chsing_file);
  OPENE(filename + ".chM", chebbd_file);
  OPENE(filename + ".chtrre", chtrre_file);
  OPENE(filename + ".chtrim", chtrim_file);

  //If we need more than 10 digits of precision (say), we're in deep
  //trouble anyway
  realmu_file.precision(15);
  imagmu_file.precision(15);
  disout_file.precision(15);
  timing_file.precision(15); //this really doesn't need this much precision

  //==============================================================
  //Sanity check (are we doing runtime checks)
  
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

  std::vector<double> hzs;
  IQMPO H;
  double opnorm_bound;
  //should probably be case/switch
  if     (model == "rfheis") { std::tie(H, hzs, opnorm_bound) = rfheis(sites, hz, e); }
  else if(model == "xx")     { std::tie(H, hzs, opnorm_bound) = XX(sites, hz, e); }
  else if(model == "rpara")  { std::tie(H, hzs, opnorm_bound) = rpara(sites, hz, e); }
  else {std::cerr << "unknown model " << model ; exit(-1);}
  
  H *= Q/opnorm_bound;

  //write the disorder to file
  for(int b = 0; b < L; ++b) { disout_file << hzs[b] << " "; }
  disout_file << "\n";

  //construct current mpo
  auto j_ampo = AutoMPO(sites);
  for(int b = 1; b < L; ++b)
    {
      auto coeff = im*0.5;
      j_ampo += coeff,"S+",b,"S-",b+1;
      j_ampo += -coeff,"S-",b,"S+",b+1;
    }

  //==============================================================
  // compute mu
  auto j = IQMPO(j_ampo);

  auto t0 = std::chrono::high_resolution_clock::now();
  if (dangler){
    auto mu = dangler_all_double_mu(H, j, realmu_file, imagmu_file, chebbd_file, chsing_file, chtrre_file, chtrim_file, N, Maxm, cutoff, 1);
  } else {
    auto mu = all_double_mu(H, j, realmu_file, imagmu_file, chebbd_file, chsing_file, chtrre_file, chtrim_file, N, Maxm, cutoff, 1);
  }
  auto t1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> computation_time = t1 - t0;
  timing_file << computation_time.count() << "\n";
}

#endif //ifndef TEST
