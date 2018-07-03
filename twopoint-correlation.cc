//this should probably be in julia with the C++ ffi

#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <random>
#include <complex>
#include <iostream>
#include <fstream>
#include <chrono>

//for making directories
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

//from getopt example. Do I actually need all of these?
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>

#include "runtime-check.h"
#include "globals.h"
#include "util.h"
		
using namespace itensor;

int main(int argc, char **argv)
{
  std::string output_filename  = "/tmp/conductivity";
  std::string dangler_filename = "/tmp/conductivity.chMPA";
  std::string siteset_filename = "/tmp/conductivity.sites";
  std::string model = "rfheis";
  std::vector<std::string> quantities;
  
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
  if (optind < argc) {
    while (optind < argc) { quantities.push_back(argv[optind++]); } }
  else {
    Error("Specify correlation to compute"); }
  
  SpinHalf sites;
  readFromFile(siteset_filename,sites);
  MPOt<IQTensor> Tn(sites);
  readFromFile(dangler_filename, Tn);
  
  int L = Tn.N();
  
  for(int i1 = 0; i1 < quantities.size(); i1++) {
    std::string q1 = quantities[i1];
    for(int i2 = 0; i2 < quantities.size(); i2++) {
      std::string q2 = quantities[i2];
      std::string dir = output_filename + q1 + q2;
      /* Unix only. Let's be honest: this is a stupid hack and the
	 right thing to do is to use HDF5.*/
      int returnval = mkdir(dir.c_str(), 0755);
      if(0 != returnval) {
	std::cerr << "Could not make directory: errno " << errno << "\n";
	return(returnval);
      }
      
      //std::vector<std::vector<Tensor>> mu;
      auto mu = twopoint_correlation(Tn, q1, q2);
      assert(mu.size() == L);
      assert(mu[0].size() == L); //just a spot-check...
      
      //ftm assumes one-site observables
      for(int j1 = 0; j1 < L; j1++) {
	for(int j2 = 0; j2 < L; j2++) {
	  auto f = dir + "/" + std::to_string(j1) + "-" + std::to_string(j2);
	  OPENE(f + ".re", realmu_file);
	  OPENE(f + ".im", imagmu_file);
	  write_doubleKPM(mu[j1][j2], realmu_file, imagmu_file);
	}
      }
    }
  }
}