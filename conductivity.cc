#include "itensor/all.h"
#include <random>
#include <complex>
#include <iostream>
#include "runtime-check.h"
// Do we do runtime checks? Defined via makefile
//#define CHECK false
// Do we do slow checks (the ones that require creating a whole dense
// operator)?
// #define SLOW_CHECK true 
			

using namespace itensor;
MPOt<IQTensor>
eye(SiteSet sites)
{
  auto I_ampo = AutoMPO(sites);
  I_ampo += 1.0, "Id", 1;
  MPOt<IQTensor> I = IQMPO(I_ampo);
  return I;
}

// slight modification of overlap(psi, H,K, phi) from itensor mpo.cc.
// overlap is almost but not quite what I want: it doesn't get the
// prime level quite right, hence messes up the contractions.
template <class Tensor>
std::complex<double>
single_mu(MPOt<Tensor> const& Tn,
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
  for(int i = 2; i < N; i++)
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
Tensor
full(MPOt<Tensor> B)
{
  auto L = B.A(1);
  for(int i = 2; i < B.N(); i++) { L = L*B.A(i); }
  return L;
}


template <class Tensor>
double
operatornorm(MPOt<Tensor> const& B)
{
  //TODO
}

template <class Tensor>
void
advance_chebyshevs(MPOt<Tensor>& Tn,
		   MPOt<Tensor>& Tnm1,
		   MPOt<Tensor> const& H,
		   int Maxm)
{
  MPOt<Tensor> C;
  // Multiply Tn by H and store the result in C. (Fortran idiom.) This
  // is probably the slowest part.
  nmultMPO(Tn, H, C, {"Maxm",Maxm}); 
  C *= 2;

  //Probably second slowest step.
  C.plusEq( (-1)* Tnm1, {"Maxm",Maxm}) ;
  
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

// Miles probably has a matrix class that's better than this
// vector<vector<>> construction, but I'm not going to worry about it:
// this is the simplest thing.
template <class Tensor>
std::vector<std::vector<std::complex<double>> >
all_mu(MPOt<Tensor> const& H,
       MPOt<Tensor> const& j,
       int N, int Maxm, int prog_per)
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
    
  auto I_ampo = AutoMPO(H.sites());
  I_ampo += 1.0, "Id", 1;
  auto I = IQMPO(I_ampo);
  
  MPOt<Tensor> Tn   = H ;
  MPOt<Tensor> Tnm1 = I;
  MPOt<Tensor> Tm, Tmm1, C;

#if CHECK
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
#if CHECK
    Tnm2 = Tnm1; // Copies. Store these guys so we can use them to
    Unm2 = Unm1; // check that we got the recursion relation right
#endif
    for (int m = n; m < N; m++) {
      mu[n][m] = mu[m][n] = single_mu(Tn,Tm,j);

      // Once again, this 0 == m case code copies. That's fine!
      // Those're really small MPOs, and we only do it once! We could
      // probably be clever and *not* copy, but that would require
      // cleverness and care in a whole bunch of places and seems
      // awfully easy to mess up in some subtle way.
      if(0 == m) { Tm   = H; Tmm1 = I; }
      // else get the Chebyshevs moved up for the next iteration
      else       { advance_chebyshevs(Tm, Tmm1, H, Maxm); }
#if CHECK
      if(0 == m) { Um   = 2*H; Umm1 = I; }
      else
	{
	  advance_chebyshevs(Um, Umm1, H, Maxm);
	  check_chebyshevs(H, H2, Tn, Tnm1, Tnm2, Um, Umm1, Umm2);
	} 
#endif
    }
    
    if(0 == n) { Tn   = H; Tnm1 = I; }
    else       { advance_chebyshevs(Tn, Tnm1, H, Maxm); }
#if CHECK
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


    
int main()
{

  int L = 8; // System size
  int N = 8; // Chebyshev order. This is very small: be nice to use
	      // more like N = 100 or even N = 1000, but that's not
	      // great for testing purposes (note that this si
  
  int Maxm = 32; // Bond dimension cutoff. This is also pretty small,
		 // though we may be able to get away with it
    
  int s = 0;

#if CHECK
  std::cout << "runtime checks inoperative" << "\n";
  exit(1);

  if (8 < L ||  32 < N)
    {
      std::cerr << "runtime checks are a bad idea with large L, N; L, N ";
      std::cerr << L << " " << N << "\n";
    }
    
#endif

  
  std::complex<float> im(0,1); //1i
  
  //set up random number generation
  //cf https://isocpp.org/files/papers/n3551.pdf
  
  std::default_random_engine e{s};
  std::uniform_real_distribution<double> d(0.0,1.0);
  
  SiteSet sites = SpinHalf(L);

  //construct our random-field Heisenberg
  //TODO: check factor of two
  auto H_ampo = AutoMPO(sites);
  for(int b = 1; b < L; ++b)
    {
      H_ampo += 0.5,"S+",b,"S-",b+1;
      H_ampo += 0.5,"S-",b,"S+",b+1;
      H_ampo +=     "Sz",b,"Sz",b+1;
    }
  for(int b = 1; b <= L; ++b)
    {
      double hz = d(e);
      H_ampo += hz, "Sz",b;
    }
  auto H = IQMPO(H_ampo);

  //construct current mpo
  auto j_ampo = AutoMPO(sites);
  for(int b = 1; b < L; ++b)
    {
      auto coeff = 0.5;
      j_ampo += coeff,"S+",b,"S-",b+1;
      j_ampo += -coeff,"S-",b,"S+",b+1;
    }
  auto j = IQMPO(j_ampo);
  auto mu = all_mu(H, j, N, Maxm, 1);
  /*
  for (int n = 0; n < L; n++) {
    for (int m = 0; m < L; m++) {
      std::cout << n << " " << m << " " << mu[n][m] << "\n";
    }
  }
  */
}

