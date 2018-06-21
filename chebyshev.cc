#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <chrono>
#include <string>

#include "util.cc"
#include "algebra.cc"
#include "globals.h"

#ifndef MPOKPM_CHEBYSHEV
#define MPOKPM_CHEBYSHEV
  
using namespace itensor;

// compute overlap of Tn with op.

template <class Tensor>
Tensor
single_mu(MPOt<Tensor> const& Tn,
	  MPOt<Tensor> const& opp )
{

  //just a transpose
  //first copy
  MPOt<Tensor> op = opp;
  //then change primelevels
  op.mapprime(1,2);
  op.mapprime(0,1);
  op.mapprime(2,0);

  int N = op.N();
  assert(N == Tn.N());
    
  auto L = Tn.A(N) * op.A(N);
  for(int i=N-1; i >= 1; --i) { L = L * Tn.A(i) * op.A(i); }
  return L;
}

/* slight modification of overlap(psi, H,K, phi) from itensor mpo.cc.
overlap is almost but not quite what I want: it doesn't get the prime
level quite right, hence messes up the contractions.

nominally computes a single coefficient mu = Tr[j Tn j Tm], but when
applied to an algebra of Chebyshevs gives all the mu */

template <class Tensor>
Tensor
double_mu(MPOt<Tensor> const& Tn,
	  MPOt<Tensor> const& Tm,
	  MPOt<Tensor> const& j)
	  
{
  if(Tm.N() != Tn.N() || j.N() != Tn.N()) Error("Mismatched N in single_mu");
  auto N = Tn.N();
  auto Tndag = Tn;
  
  // I'm pretty sure this copies at least the indices (deduce this
  // from the fact that Miles uses it in overlap, which one wouldn't
  // expect to mess with the prime level of the input)
  auto jp    = j;  
  auto jpdag = j; //TODO this never actually gets daggered
  
  jp.mapprime(1,2);
  jp.mapprime(0,1);

  Tndag.mapprime(1,3);
  Tndag.mapprime(0,2);

  jpdag.mapprime(0,3);
  jpdag.mapprime(1,0); 

  //TODO start from right (so I'm not dragging a NxN matrix all the
  //way across in the algebra case)
  
  //scales as m^2 k^2 d per Miles in mpo.cc. What in my case?
  auto L = jpdag.A(N) * Tndag.A(N) * jp.A(N) * Tm.A(N);
  for(int i = N-1; i >= 1; --i)
    {
      //scales as m^3 k^2 d + m^2 k^3 d^2. What in my case?
      L = L * jpdag.A(i) * Tndag.A(i) * jp.A(i) * Tm.A(i);
    }
  //PrintData(L);
  return L;
}


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
  
  if (writep) { writeToFile(filename + ".chMPA", cheb); }
  return cheb;
}

#endif //#ifndef MPOKPM_CHEBYSHEV
