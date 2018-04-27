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

/*
slight modification of overlap(psi, H,K, phi) from itensor mpo.cc.
overlap is almost but not quite what I want: it doesn't get the
prime level quite right, hence messes up the contractions.

nominally computes a single coefficient mu = Tr[j Tn j Tm], but when
applied to an algebra of Chebyshevs gives all the mu
*/
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
next_chebyshev(MPOt<Tensor> const& Tn,
	       MPOt<Tensor> const& Tnm1,
	       MPOt<Tensor> const& H,
	       int Maxm, double cutoff)
{
  MPOt<Tensor> C;
 
  // Multiply Tn by H and store the result in C. (Fortran idiom.) This
  // is probably the slowest part.

  int L = Tn.N();
  nmultMPO(Tn, H, C, {"Cutoff",cutoff}); 
  C *= 2;
  C.position(L,{"Maxm", Maxm, "Cutoff", 1e-14});
  C.position(1,{"Maxm", Maxm, "Cutoff", 1e-14});

  //Probably second slowest step.
  C.plusEq( (-1)* Tnm1, {"Cutoff",1e-30}) ;

  C.position(L, {"Maxm", Maxm, "Cutoff", 1e-14});
  C.position(1, {"Maxm", Maxm, "Cutoff", 1e-14});
  return C;
}
template <class Tensor>
void
advance_chebyshevs(MPOt<Tensor>& Tn,
		   MPOt<Tensor>& Tnm1,
		   MPOt<Tensor> const& H,
		   int Maxm)
{
  std::cerr << "Called advance_chebyshevs. I think this is buggy!";
  exit(1);
  MPOt<Tensor> C;

  // We don't want to truncate during nmultMPO---that would truncate
  // at each bond as we zip up from left to right, which is very
  // wrong. We do still want nix the numerical noise,
  // though. Frobenius norm is submultiplicative, so we take a cutoff
  // like 1e-14 * (upper bound on Frobenius norm of H*Tn).
  double cutoff = 1e-14;//*norm(H)*norm(Tn); 
  
  // Multiply Tn by H and store the result in C. (Fortran idiom.) This
  // is probably the slowest part.
  
  nmultMPO(Tn, H, C, {"Cutoff",cutoff}); 
  C *= 2;
  C.orthogonalize({"Maxm", Maxm, "Cutoff", 1e-14});

  //Probably second slowest step.
  C.plusEq( (-1)* Tnm1, {"Cutoff",cutoff}) ;
  C.orthogonalize({"Maxm", Maxm, "Cutoff", 1e-14});
  
  /* This next bit is weird but cool to me, the C++ newbie.

     First: take an rvalue reference to Tn. This is like a Rust
     borrow: we decide that we're going to think of the value of Tn as
     a special hidden thing with only one reference. This invalidates
     all other ways of accessing/setting Tn, so we can mess with it
     with impunity.

     We then hand this reference off to Tn. (What happens to the old
     Tn? Does it get deallocated?)

     
     Why not an lvalue? If we were to do this via lvalue reference
     like 

         Tnm1 = static_cast<MPOt<Tensor> &> C 

     then when we set B = C we'll be setting A as well. That would be
     bad and wrong.
  */
  Tnm1 = static_cast<MPOt<Tensor> &&>(Tn);
  Tn   = static_cast<MPOt<Tensor> &&>(C);

  /* So here's my question: When C goes out of scope, does the memory
     get deallocated? I think not, because it's no longer associated
     with C: that's the point of doing static_cast<MPOt<Tensor> &&> C.
     But am I wrong?
  */
  
}

template <class Tensor>
MPOt<Tensor>
advance_dangler_chebyshevs(MPOt<Tensor> cheb, int n0, MPOt<Tensor> H,
			   int N, int Maxm, double cutoff,
			   std::ofstream& chebbd_file,
			   int prog_per)
{
  auto I = eye(H.sites());

  MPOt<Tensor> iter = oplus(I,H);

  cheb.orthogonalize();
  iter.orthogonalize();

  //n is the number of Chebyshevs we've already computed.
  //note that this n is notionally a 1-index: "T1" is identity
  for(int n = n0; n < n0 + N; n++)
    {
      if (n % prog_per == 0) { std::cout << n << " " << maxM(cheb) << "\n"; }
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
      
      nmultMPAlgebra(cheb, iter, iterbc, store, {"Maxm", Maxm, "Cutoff", cutoff}); 
      cheb = store; //copies! Bad! Can go do rvalue magic...

      chebbd_file << n << " " << maxM(cheb) << "\n" << std::flush;
    }
  return cheb;
}

template <class Tensor>
std::vector<MPOt<Tensor>>
list_chebyshevs(MPOt<Tensor> H, int N, int Maxm, double cutoff, int prog_per)
{
  std::vector<MPOt<Tensor>> lst;
  lst.reserve(N);

  auto I = eye(H.sites());

  int L = H.N();
  MPOt<Tensor> Tn = I*pow(2.0,-L/2);
  MPOt<Tensor> Tnm1 = I;

  for(int n = 0; n < N; n++)
    {
      lst.push_back(Tn);
      if (n % prog_per == 0) { std::cout << n << " " << maxM(Tn) << "\n"; }
    
      if(0 == n) {
	Tn   = H*pow(2.0,-L/2);
	Tnm1 = I*pow(2.0,-L/2);
	Tn.orthogonalize();
	Tnm1.orthogonalize();
      } else { // get the Chebyshevs moved up for the next iteration
	MPOt<Tensor> hold = next_chebyshev(Tn, Tnm1, H, Maxm, cutoff);
	Tnm1 = Tn; //these copy: not great
	Tn = hold; 
      }
    }
  return lst;
}

template <class Tensor>
void
check_chebyshevs(MPOt<Tensor> H,
		 MPOt<Tensor> H2,
		 MPOt<Tensor> Tn,
		 MPOt<Tensor> Tnm1,
		 MPOt<Tensor> Tnm2,
		 MPOt<Tensor> Un,
		 MPOt<Tensor> Unm1,
		 MPOt<Tensor> Unm2)
{
  // note that nmultMPO imposes a default cutoff 1e-14, but (AFAICT)
  // nothing else.
  
  // TODO: bandwidth of Tn
  MPOt<Tensor> HTnm1,HUnm1; //really creative naming here
  std::cout << "before H mults\n";
  std::cout << H.N() << "\n";
  std::cout << Tnm1.N() << "\n";
  std::cout << Unm1.N() << "\n";
  nmultMPO(H, Tnm1, HTnm1); //Does this truncate by d
  nmultMPO(H, Unm1, HUnm1); //Does this truncate by d

  std::cout << "H mults finished\n";
  // this is is just a reimplimentation of the recursion relation. Not
  // such an interesting check, probably
  std::cout << "try sum\n";
  std::cout << sum(Tnm2, 2*HTnm1) << "\n";
  std::cout << "try overlap\n";
  std::cout << overlap(Tn, Tnm1) << "\n";
  check(magdiff(Tn, sum(Tnm2, 2*HTnm1)), 0, "recursion relation: T");
  check(magdiff(Un, sum(Unm2, 2*HUnm1)), 0, "recursion relation: U");

  // Pell equation.
  // This'll take more memory than (strictly speaking) it might require.
  MPOt<Tensor> Tn2, Unm12, H2Unm12;
  nmultMPO(Tn,   Tn,   Tn2);
  nmultMPO(Unm1, Unm1, Unm12);
  std::cout << "T,U squares finished";
  check(magdiff(sum(Tn2, Unm12) , sum((-1)*H2Unm12, eye(Tn2.sites()))), 0, "Pell equation");
}


//TODO: when I get runtime checks working, this needs some
template <class Tensor>
std::vector<std::complex<double>>
all_single_mu(MPOt<Tensor> const&H,
	      MPOt<Tensor> const&j,
	      std::ofstream& realmu_file,
	      std::ofstream& imagmu_file,
	      std::ofstream& chebbd_file,
	      std::ofstream& chsing_file,
	      int N, int Maxm, double cutoff, int prog_per)
{
  std::vector<std::complex<double>> mu(N,0);
  IQMPO I = eye(H.sites());

  int L = H.N();
  MPOt<Tensor> Tn = I*pow(2.0,-L/2);
  MPOt<Tensor> Tnm1 = I;

  for(int n = 0; n < N; n++)
    {
      mu[n] = single_mu(Tn, j).cplx();
      realmu_file << real(mu[n]) << " ";
      imagmu_file << imag(mu[n]) << " ";
      chebbd_file << n << " " << maxM(Tn) << "\n" << std::flush;

      auto spectrum = entanglement_spectrum(Tn, Tn.N()/2);
      for(double s : spectrum.eigs()) { if (s > 1e-14){chsing_file  << s << " " << std::flush;}}
      chsing_file << "\n";

      if (n % prog_per == 0) { std::cout << n << " " << maxM(Tn) << "\n"; }
    
      if(0 == n) {
	Tn   = H*pow(2.0,-L/2);
	Tnm1 = I*pow(2.0,-L/2);
	Tn.orthogonalize();
	Tnm1.orthogonalize();
      }
      // else get the Chebyshevs moved up for the next iteration
      else
	{
	  MPOt<Tensor> hold = next_chebyshev(Tn, Tnm1, H, Maxm, cutoff);
	  Tnm1 = Tn; //these copy: not great
	  Tn = hold; 
	}
      
    }
  return mu;
}

// Miles probably has a matrix class that's better than this
// vector<vector<>> construction, but I'm not going to worry about it:
// this is the simplest thing.
template <class Tensor>
std::vector<std::vector<std::complex<double>> >
all_double_mu(MPOt<Tensor> const& H,
	      MPOt<Tensor> const& j,
	      std::ofstream& realmu_file,
	      std::ofstream& imagmu_file,
	      std::ofstream& chebbd_file,
	      std::ofstream& chsing_file,
	      std::ofstream& chtrre_file,
	      std::ofstream& chtrim_file,
	      int N, int Maxm, double cutoff, int prog_per)
{
  //mu should be real for physical reasons. Important not to take real
  //part until we really mean to, though: imag could be clue to
  //instability.
  std::vector<std::vector<std::complex<double>>> mu;
  mu.reserve(N);

  int L = H.N();

  //initialize mu to be 0
  for (int i = 0; i < N; i++) {
    mu.push_back(std::vector<std::complex<double>>(N,0));
  }


  //TODO this is eye
  auto I_ampo = AutoMPO(H.sites());
  I_ampo += 1.0, "Id", 1;
  auto I = IQMPO(I_ampo);
  
  MPOt<Tensor> Tn   = H ;
  MPOt<Tensor> Tnm1 = I;
  MPOt<Tensor> Tm, Tmm1, C;

#ifdef CHECK
  MPOt<Tensor> Tnm2, Un, Unm1, Unm2;
  MPOt<Tensor> Tmm2, Um, Umm1, Umm2;
  MPOt<Tensor> H2;
  nmultMPO(H,H,H2);
#endif
    
  Tn = I*pow(2.0,-L/2);
  Tmm1 = I;
  for (int n = 0; n < N; n++) {
    // we have the chebyshev polynomial Tn. Now walk through the Tms
    // computing tr(Tn j Tm j). We'll advance Tn at the end after
    // we've gone through the Tms.
    
    // these two lines copy. That is as it should be!
    Tm   = Tn;
    Tmm1 = Tnm1;
#ifdef CHECK
    Tnm2 = Tnm1; // Copies. Store these guys so we can use them to
    Unm2 = Unm1; // check that we got the recursion relation right
#endif

    for (int m = 0; m < n; m++)
      {
	realmu_file << 0.0 << " ";
	imagmu_file << 0.0 << " ";
      }
    for (int m = n; m < N; m++) {
      mu[n][m] = mu[m][n] = double_mu(Tn,Tm,j).cplx();
      realmu_file << real(mu[m][n]) << " ";
      imagmu_file << imag(mu[m][n]) << " ";

      // Once again, this 0 == m case code copies. That's fine!
      //
      // probably be clever and *not* copy, but that would require
      // cleverness and care in a whole bunch of places and seems
      // awfully easy to mess up in some subtle way.
      if(0 == m) {
	  Tm   = H*pow(2.0,-L/2);
	  Tmm1 = I*pow(2.0,-L/2);
	  Tm.orthogonalize();
	  Tmm1.orthogonalize();
	}
      // else get the Chebyshevs moved up for the next iteration
      else
	{
	  MPOt<Tensor> hold = next_chebyshev(Tm, Tmm1, H, Maxm, cutoff);
	  Tmm1 = Tm; //these copy: not great
	  Tm = hold;
	}

#ifdef CHECK
      if(0 == m) { Um   = 2*H; Umm1 = I; }
      else
	{
	  
	  MPOt<Tensor> hold = next_chebyshev(Um, Umm1, H, Maxm, cutoff);
	  Umm1 = Tn; //these copy: not great
	  Um = hold; 
	  check_chebyshevs(H, H2, Tn, Tnm1, Tnm2, Um, Umm1, Umm2);
	} 
#endif
    }

    realmu_file << "\n";
    imagmu_file << "\n";
      
    if(0 == n) {
      Tn   = H*pow(2.0,-L/2);
      Tnm1 = I*pow(2.0,-L/2);
      Tn.orthogonalize();
      Tnm1.orthogonalize();
    }

    else       {
      MPOt<Tensor> hold = next_chebyshev(Tn, Tnm1, H, Maxm, cutoff);
      Tnm1 = Tn; //these copy: not great
      Tn = hold; 
    }
#ifdef CHECK
    if(0 == n) { Un   = 2*H; Unm1 = I; }
    else
      {
	advance_chebyshevs(Um, Umm1, H, Maxm);
	check_chebyshevs(H, H2, Tn, Tnm1, Tnm2, Un, Unm1, Unm2);
      } 
#endif
    
    chebbd_file << n << " " << maxM(Tn) << "\n" << std::flush;
    std::complex<double> chtr = single_mu(Tn, I).cplx();
    chtrre_file << real(chtr) << " ";
    chtrim_file << real(chtr) << " " ;
    if (n % prog_per == 0) { std::cout << n << " " << maxM(Tn) << "\n"; }
  }
  return mu;
}


template <class Tensor>
std::vector<std::vector<std::complex<double>> >
memoryprofligate_all_double_mu(MPOt<Tensor> const& H,
			       MPOt<Tensor> const& j,
			       std::ofstream& realmu_file,
			       std::ofstream& imagmu_file,
			       std::ofstream& chebbd_file,
			       std::ofstream& chsing_file,
			       std::ofstream& chtrre_file,
			       std::ofstream& chtrim_file,
			       int N, int Maxm, double cutoff, int prog_per)
{
  //mu should be real for physical reasons. Important not to take real
  //part until we really mean to, though: imag could be clue to
  //instability.
  std::vector<std::vector<std::complex<double>>> mu;
  mu.reserve(N);

  //initialize mu to be 0
  for (int i = 0; i < N; i++) {
    mu.push_back(std::vector<std::complex<double>>(N,0));
  }

  //TODO this is eye
  auto I_ampo = AutoMPO(H.sites());
  I_ampo += 1.0, "Id", 1;
  auto I = IQMPO(I_ampo);
  
  std::vector<MPOt<Tensor>> Tn = list_chebyshevs(H, N, Maxm, cutoff, prog_per);
    
  for (int n = 0; n < N; n++) {
    for (int m = 0; m < n; m++)
      {
	realmu_file << 0.0 << " ";
	imagmu_file << 0.0 << " ";
      }
    chebbd_file << n << " " << maxM(Tn[n]) << "\n" << std::flush;
    for (int m = n; m < N; m++) {
      mu[n][m] = mu[m][n] = double_mu(Tn[n],Tn[m],j).cplx();
      realmu_file << real(mu[m][n]) << " ";
      imagmu_file << imag(mu[m][n]) << " ";
    }
    
    realmu_file << "\n";
    imagmu_file << "\n";
    
    chebbd_file << n << " " << maxM(Tn[n]) << "\n" << std::flush;
    std::complex<double> chtr = single_mu(Tn[n], I).cplx();
    chtrre_file << real(chtr);
    chtrim_file << real(chtr);
    if (n % prog_per == 0) { std::cout << n << "\n"; }
  }
  return mu;
}


template <class Tensor>
std::vector<std::vector<std::complex<double>> >
dangler_all_double_mu(MPOt<Tensor> const& H,
		      MPOt<Tensor> const& j,
		      std::string const& filename,
		      int N, int Maxm, double cutoff, int prog_per)
{

  
  std::vector<std::vector<std::complex<double>>> vecmu;
  
  // macro OPENE declares second arg and initializes with filehandle to
  // firstarg. Defined in util.cc
  OPENE(filename + ".tim", timing_file);
  OPENE(filename + ".chM", chebbd_file);

  std::cout << "using dangler\n";


  int N_save = 16;
  int N_step = floor((N-2)/N_save);
  int N_rem = (N-2) % N_save;
  
  int N_sofar = 2;
  
  auto I = eye(H.sites());
  MPOt<Tensor> cheb = oplus(I,H);

  //normalization
  int L = H.N();
  cheb = cheb*pow(2.0, -L/2);
  
  for(int i = 0; i < N_step; i++){
    
    cheb = advance_dangler_chebyshevs(cheb, N_sofar, H, N_save, Maxm, cutoff, chebbd_file, prog_per);
    N_sofar += N_save;
    
    IQTensor chtr = single_mu(cheb, I);
    IQTensor mu = double_mu(cheb, cheb, j);
    
    //performance information
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> computation_time = t1 - t0;
    timing_file << N_sofar << " " << computation_time.count() << "\n";
    
    OPENE(filename + "." + std::to_string(N_sofar) + ".re",  realmu_file);
    OPENE(filename + "." + std::to_string(N_sofar) + ".im",  imagmu_file);
    OPENE(filename + "." + std::to_string(N_sofar) + ".chtrre", chtrre_file);
    OPENE(filename + "." + std::to_string(N_sofar) + ".chtrim", chtrim_file);
    
    write_singleKPM(chtr, N_sofar, chtrre_file, chtrim_file);
    vecmu = write_doubleKPM(mu, N_sofar, realmu_file, imagmu_file);
    std::cout << "\n";
  }


  if (N_rem != 0) {
    
    cheb = advance_dangler_chebyshevs(cheb, N_sofar, H, N_rem, Maxm, cutoff, chebbd_file, prog_per);
    N_sofar += N_rem;
    IQTensor chtr = single_mu(cheb, I);
    IQTensor mu = double_mu(cheb, cheb, j);
    
    
    //performance information
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> computation_time = t1 - t0;
    timing_file << N_sofar << " " << computation_time.count() << "\n";

    
    OPENE(filename + ".re",  realmu_file);
    OPENE(filename + ".im",  imagmu_file);
    OPENE(filename + ".chtrre", chtrre_file);
    OPENE(filename + ".chtrim", chtrim_file);
    
    write_singleKPM(chtr, N_sofar, chtrre_file, chtrim_file);
    vecmu = write_doubleKPM(mu, N_sofar, realmu_file, imagmu_file);
  }
  
  return vecmu;
}

#endif //#ifndef MPOKPM_CHEBYSHEV
