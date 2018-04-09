#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "util.cc"

#ifndef MPOKPM_CHEBYSHEV
#define MPOKPM_CHEBYSHEV

using namespace itensor;

// compute overlap of Tn with op.

template <class Tensor>
std::complex<double>
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
    
  auto L = Tn.A(1) * op.A(1);
  for(int i=2; i <= N; ++i) { L = L * Tn.A(i) * op.A(i); }
  auto z = L.cplx();
  return z;
}
  
// slight modification of overlap(psi, H,K, phi) from itensor mpo.cc.
// overlap is almost but not quite what I want: it doesn't get the
// prime level quite right, hence messes up the contractions.
template <class Tensor>
std::complex<double>
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
  auto jpdag = j;
  
  jp.mapprime(1,2);
  jp.mapprime(0,1);

  Tndag.mapprime(1,3);
  Tndag.mapprime(0,2);

  jpdag.mapprime(0,3);
  jpdag.mapprime(1,0); 

  //scales as m^2 k^2 d per Miles in mpo.cc. What in my case?
  auto L = jpdag.A(1) * Tndag.A(1) * jp.A(1) * Tm.A(1);
  for(int i = 2; i < N; ++i)
    {
      //scales as m^3 k^2 d + m^2 k^3 d^2. What in my case?
      L = L * jpdag.A(i) * Tndag.A(i) * jp.A(i) * Tm.A(i);
    }
  //scales as m^2 k^2 d
  L = L* jpdag.A(N) * Tndag.A(N) * jp.A(N) * Tm.A(N);
  //PrintData(L);
  auto z = L.cplx();
  return z;
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
      mu[n] = single_mu(Tn, j);
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
  
  MPOt<Tensor> Tn   = H ;
  MPOt<Tensor> Tnm1 = I;
  MPOt<Tensor> Tm, Tmm1, C;

#ifdef CHECK
  MPOt<Tensor> Tnm2, Un, Unm1, Unm2;
  MPOt<Tensor> Tmm2, Um, Umm1, Umm2;
  MPOt<Tensor> H2;
  nmultMPO(H,H,H2);
#endif
    
  Tn = I;
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
      mu[n][m] = mu[m][n] = double_mu(Tn,Tm,j);
      realmu_file << real(mu[m][n]) << " ";
      imagmu_file << imag(mu[m][n]) << " ";

      // Once again, this 0 == m case code copies. That's fine!
      //
      // probably be clever and *not* copy, but that would require
      // cleverness and care in a whole bunch of places and seems
      // awfully easy to mess up in some subtle way.
      if(0 == m) { Tm   = H; Tmm1 = I; }
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
      
    if(0 == n) { Tn   = H; Tnm1 = I; }
    else       {
      MPOt<Tensor> hold = next_chebyshev(Tn, Tnm1, H, Maxm, cutoff);
      Tnm1 = Tn; //these copy: not great
      Tn = hold; 
    }
    chebbd_file << n << " " << maxM(Tn) << "\n" << std::flush;
    
#ifdef CHECK
    if(0 == n) { Un   = 2*H; Unm1 = I; }
    else
      {
	advance_chebyshevs(Um, Umm1, H, Maxm);
	check_chebyshevs(H, H2, Tn, Tnm1, Tnm2, Un, Unm1, Unm2);
      } 
#endif
    
    if (n % prog_per == 0) { std::cout << n << "\n"; }
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
      mu[n][m] = mu[m][n] = double_mu(Tn[n],Tn[m],j);
      realmu_file << real(mu[m][n]) << " ";
      imagmu_file << imag(mu[m][n]) << " ";
    }
    
    realmu_file << "\n";
    imagmu_file << "\n";
    
    if (n % prog_per == 0) { std::cout << n << "\n"; }
  }
  return mu;
}

#endif //#ifndef MPOKPM_CHEBYSHEV
