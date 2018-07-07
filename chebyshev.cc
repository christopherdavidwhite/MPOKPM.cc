#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <chrono>
#include <string>

#include "util.h"
#include "algebra.cc"
#include "globals.h"

#ifndef MPOKPM_CHEBYSHEV
#define MPOKPM_CHEBYSHEV
  
using namespace itensor;


template <class Tensor>
MPOt<Tensor>
advance_dangler_chebyshevs(MPOt<Tensor>& cheb, MPOt<Tensor> iter,
			   int n0, int Maxm, double cutoff, bool fit, int sweeps)
{
  //n is the number of Chebyshevs we've already computed.
  //note that this n is notionally a 1-index: "T1" is identity
  int n = n0;
  MPOt<Tensor> store;

  auto i = findtype(cheb.A(1), Dangler);
  auto j = findtype(iter.A(1), Dangler);

  IQIndex id = i.dag();
  IQIndex jd = j.dag();
  IQIndex k  = IQIndex("iq", Index("i",n+1, Dangler, 0), QN(0));

  //bc implementing Cheb. iterations to tensor on
  auto iterbc = Tensor(k, id,jd);
  //keep things multiplied by identity
  for(int l = 1; l <= n; l++) { iterbc.set(k(l),id(l),jd(1), 1.0); }
  iterbc.set(k(n+1), id(n),   jd(2), 2.0);  //2*H*Tn
  iterbc.set(k(n+1), id(n-1), jd(1), -1.0); // - Tn-1

  if (fit)
    fitmultMPAlgebra(cheb, iter, iterbc, store,
		     {"Maxm", Maxm, "Cutoff", cutoff}, sweeps);
  else
    nmultMPAlgebra(cheb, iter, iterbc, store,
		       {"Maxm", Maxm, "Cutoff", cutoff});
  cheb = store; //copies! Bad! Can go do rvalue magic...

  return cheb;
}

template <class Tensor>
MPOt<Tensor>
listandwrite_dangler(MPOt<Tensor> const& H,
		     std::string const& filename,
		     int N,
		     int Maxm, double cutoff, bool fit, int sweeps,
		     int prog_per)
{
  // macro OPENE declares second arg and initializes with filehandle to
  // firstarg. Defined in util.cc
  OPENE(filename + ".tim", timing_file);
  OPENE(filename + ".chM", chebbd_file);

  //simulation start time
  std::chrono::time_point<std::chrono::high_resolution_clock> t0;
  t0 = std::chrono::high_resolution_clock::now();
  std::cout << "using dangler\n";

  int N_sofar = 2;
  bool writep = true;
  
  auto I = eye(H.sites());
  MPOt<Tensor> cheb = oplus(I,H);
  MPOt<Tensor> iter = oplus(I,H);

  //normalization
  int L = H.N();
  cheb = cheb*pow(2.0, -L/2);
  cheb.orthogonalize();
  iter.orthogonalize();

  if(0 < N){
    for(int i = 0; i < N; i++){
      cheb = advance_dangler_chebyshevs(cheb, iter, N_sofar, Maxm, cutoff, fit, sweeps);
      N_sofar++;
      
      //if we've hit the ceiling, quit: results will be junk.
      if(Maxm <= maxM(cheb)) { writep = false; break;}
      
      chebbd_file << N_sofar << " " << maxM(cheb) << "\n" << std::flush;
      
      //performance information
      auto t1 = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> computation_time = t1 - t0;
      timing_file << N_sofar << " " << computation_time.count() << "\n" << std::flush;
      if ((0 == N_sofar % prog_per) && writep) { writeToFile(filename + ".chMPA", cheb); }
    }
  }
  
  if (writep) { writeToFile(filename + ".chMPA", cheb); }
  return cheb;
}

#endif //#ifndef MPOKPM_CHEBYSHEV
