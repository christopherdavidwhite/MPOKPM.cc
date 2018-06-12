//this should probably be in julia with the C++ ffi

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
#include "chebyshev.cc"
#include "globals.h"
		
using namespace itensor;

#ifndef TESTING //test has its own main function
int main(int argc, char **argv)
{
  std::string output_filename  = "/tmp/conductivity";
  std::string dangler_filename = "/tmp/conductivity.chMPA";
  std::string siteset_filename = "/tmp/conductivity.sites";
  std::string model = "rfheis";
  
  int option_index = 0;
  while (1) {
      static struct option long_options[] =
        {
          {"sites",   required_argument, 0, 't'},
	  {"dangler", required_argument, 0, 'd'},
	  {"model",   required_argument, 0, 'm'},
	  {"output",  required_argument, 0, 'o'},
        };
      
      char c = getopt_long (argc, argv, "t:d:m:o:", long_options, &option_index);
      
      if (c == -1) break;
      switch (c) {
        case 0:   if (long_options[option_index].flag != 0) break;
        case 'o': output_filename  = optarg; break;
	case 'm': model            = optarg; break;
	case 'd': dangler_filename = optarg; break;
	case 't': siteset_filename = optarg; break;
        case '?': break;
        default: abort ();
      }
  }

  SpinHalf sites;
  readFromFile(siteset_filename,sites);
  MPOt<IQTensor> Tn(sites);
  readFromFile(dangler_filename, Tn);
  
  int L = Tn.N();
  std::cout << "L " << L << "\n";

  //this will depend on model
  auto j_ampo = AutoMPO(Tn.sites());
  for(int b = 1; b < L; ++b) {
      auto coeff = im*0.5;
      j_ampo +=  coeff,"S+",b,"S-",b+1;
      j_ampo += -coeff,"S-",b,"S+",b+1;
  }
  auto j = IQMPO(j_ampo);

  IQTensor mu = double_mu(Tn, Tn, j);

  OPENE(output_filename + ".re",  realmu_file);
  OPENE(output_filename + ".im",  imagmu_file);
  write_doubleKPM(mu, realmu_file, imagmu_file );
}

#endif //ifndef TEST
