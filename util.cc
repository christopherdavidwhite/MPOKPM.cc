// Miscellaneous utilities.

#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"
#include "globals.h"
#include "util.h"

#ifndef MPOKPM_UTIL
#define MPOKPM_UTIL
using namespace itensor;

MPOt<IQTensor>
eye(SiteSet sites)
{
  auto I_ampo = AutoMPO(sites);
  I_ampo += 1.0, "Id", 1;
  MPOt<IQTensor> I = IQMPO(I_ampo);
  return I;
}


template <class Tensor>
Tensor
full(MPOt<Tensor> B)
{
  auto L = B.A(1);
  for(int i = 2; i < B.N(); i++) { L = L*B.A(i); }
  return L;
}

template IQTensor full(MPOt<IQTensor> B);

std::tuple<IQMPO,std::vector<double>, double>
XX(SiteSet sites, double hz, std::default_random_engine e)
{
  //cf https://isocpp.org/files/papers/n3551.pdf
  auto L = sites.N();
  std::vector<double> fields; //need for return...

  auto H_ampo = AutoMPO(sites);
  //++b or b++?
  for(int b = 1; b < L; ++b)
    {
      fields.push_back(0.0);
      H_ampo += 0.5,"S+",b,"S-",b+1;
      H_ampo += 0.5,"S-",b,"S+",b+1;
    }

  auto H = IQMPO(H_ampo);
  
  //operator norm is bounded by sum of operator norms of terms
  double opnorm_bound = (L-1);
  return std::make_tuple(H, fields, opnorm_bound);
}

//construct our random-field Heisenberg Hamiltonian
//returns a tuple of IQMPO, fields, and upper bound on operator norm
std::tuple<IQMPO,std::vector<double>, double>
rfheis(SiteSet sites, double hz, std::default_random_engine e)
{
  //cf https://isocpp.org/files/papers/n3551.pdf
  std::uniform_real_distribution<double> d(0.0,1.0);
  std::vector<double> fields;
  auto L = sites.N();

  auto H_ampo = AutoMPO(sites);
  //++b or b++?
  for(int b = 1; b < L; ++b)
    {
      H_ampo += 0.5,"S+",b,"S-",b+1;
      H_ampo += 0.5,"S-",b,"S+",b+1;
      H_ampo += 1.0,"Sz",b,"Sz",b+1;
    }

  double sumhzj = 0;
  for(int b = 1; b <= L; ++b)
    {
      double hzj = hz * (2*d(e) - 1);
      sumhzj += std::abs(hzj);
      fields.push_back(hzj);
      H_ampo += hzj, "Sz",b;
    }
  auto H = IQMPO(H_ampo);

  //operator norm is bounded by sum of operator norms of terms
  double opnorm_bound = 3*(L-1)*0.25 + sumhzj*0.5;
  return std::make_tuple(H, fields, opnorm_bound);
}

std::tuple<IQMPO,std::vector<double>, double>
rf_2NJW(SiteSet sites, double W, std::default_random_engine e)
{
  std::normal_distribution<double> d;
  std::vector<double> fields;
  auto L = sites.N();

  auto H_ampo = AutoMPO(sites);
  //++b or b++?
  for(int b = 1; b < L; ++b)
    {
      H_ampo += 0.5,"S+",b,"S-",b+1;
      H_ampo += 0.5,"S-",b,"S+",b+1;
      H_ampo += 1.0,"Sz",b,"Sz",b+1;
    }
  for(int b = 1; b < L-1; ++b)
    {
      H_ampo += 1.0,"S+",b,"Sz",b+1,"S-",b+2;
      H_ampo += 1.0,"S-",b,"Sz",b+1,"S+",b+2;
    }
  
  double sumhzj = 0;
  for(int b = 1; b <= L; ++b)
    {
      double hzj = W*d(e);
      sumhzj += std::abs(hzj);
      fields.push_back(hzj);
      H_ampo += hzj, "Sz",b;
    }
  auto H = IQMPO(H_ampo);

  //operator norm is bounded by sum of operator norms of terms
  double opnorm_bound = 3*0.25*(L-1) + 2*0.25*(L-2) + sumhzj*0.5;
  return std::make_tuple(H, fields, opnorm_bound);
}


std::tuple<IQMPO,std::vector<double>, double>
rpara(SiteSet sites, double hz, std::default_random_engine e)
{
  //cf https://isocpp.org/files/papers/n3551.pdf
  std::uniform_real_distribution<double> d(0.0,1.0);
  std::vector<double> fields;
  auto L = sites.N();

  auto H_ampo = AutoMPO(sites);
  
  double sumhzj = 0;
  for(int b = 1; b <= L; ++b)
    {
      double hzj = hz * (2*d(e) - 1);
      sumhzj += std::abs(hzj);
      fields.push_back(hzj);
      H_ampo += hzj, "Sz",b;
    }
  auto H = IQMPO(H_ampo);

  //operator norm is bounded by sum of operator norms of terms
  double opnorm_bound = 3*(L-1)*0.25 + sumhzj*0.5;
  return std::make_tuple(H, fields, opnorm_bound);
}

// cf http://itensor.org/docs.cgi?page=formulas/entanglement_mps
// Note that this changes the location of the orthogonality center!
//
// TODO need to understand C++ type system well enough to make this an
// MPS fn
template <class Tensor>
Spectrum
entanglement_spectrum(MPOt<Tensor> psi, int b)
{
  psi.position(b); 

  //Here assuming an MPS of ITensors, but same code works
  //for IQMPS by replacing ITensor -> IQTensor

  //Compute two-site wavefunction for sites (b,b+1)
  auto wf = psi.A(b)*psi.A(b+1);

  //SVD this wavefunction to get the spectrum
  //of density-matrix eigenvalues
  auto U = psi.A(b);
  IQTensor S,V;
  auto spectrum = svd(wf,U,S,V);
  return spectrum;
}

template Spectrum entanglement_spectrum(MPOt<IQTensor> psi, int b);

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

template IQTensor single_mu(MPOt<IQTensor> const& Tn, MPOt<IQTensor> const& opp);

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

template IQTensor double_mu(MPOt<IQTensor> const& Tn, MPOt<IQTensor> const& Tm, MPOt<IQTensor> const& opp);

void write_singleKPM(IQTensor chtr,
		     std::ofstream& chtrre_file,
		     std::ofstream& chtrim_file)
  
{
  IQIndex i1;
  if (1 != rank(chtr)) Error("write_singleKPM: chtr has wrong rank");

  i1 = chtr.inds()[0];
  int N = i1.m();
  
  for(int n = 1; n <= N; n++){
    chtrre_file << real(chtr.cplx(i1(n))) << " " << std::flush;  
    chtrim_file << imag(chtr.cplx(i1(n))) << " " << std::flush;
  }
}

std::vector<std::vector<std::complex<double>>>
write_doubleKPM(IQTensor mu,
		std::ofstream& realmu_file,
		std::ofstream& imagmu_file )
{ 
  if (2 != rank(mu)) Error("write_doubleKPM: mu has wrong rank");

  IQIndex i1;
  IQIndex i2;
  i1 = mu.inds()[0];
  i2 = mu.inds()[1];
  int N = i1.m();
  if (i2.m() != N) Error("write_doubleKPM: mu not square");
  
  std::vector<std::vector<std::complex<double>>> vecmu;
  vecmu.reserve(N);
  //initialize mu to be 0
  for (int i = 0; i < N; i++) { vecmu.push_back(std::vector<std::complex<double>>(N,0)); }
  
  for(int n = 1; n <= N; n++){
    for(int m = 1; m <= N; m++){
      vecmu[n-1][m-1] = mu.cplx(i1(n),i2(m));
      realmu_file << real(mu.cplx(i1(n), i2(m)))  << std::flush;
      imagmu_file << imag(mu.cplx(i1(n), i2(m))) << " " << std::flush;  
      if (m < N) {
	realmu_file << " ";
	imagmu_file << " ";
      }
    }
    if (n < N) {
      realmu_file << "\n";
      imagmu_file << "\n";
    }
  }
  return vecmu;
}

template <class Tensor>
MPOt<Tensor>
commutator(MPOt<Tensor> H,
	   MPOt<Tensor> q)
{
  MPOt<Tensor> Hq, qH, Hqb;
  nmultMPO(H, q, Hq, {"Cutoff", 1e-15});
  nmultMPO(q, H, qH, {"Cutoff", 1e-15});
  Hqb = (Hq - qH).orthogonalize();
}

template <class Tensor>
std::vector<Tensor>
current(MPOt<Tensor> H,
	int p, //diameter of local conserved quantity
	std::vector<MPOt<Tensor>> qlist )
{
  //assume qlist[0] starts at left end of chain
  std::vector<Tensor> rightwards_currents;

  int N = H.N();

  //for verifying that it is in fact a conserved quantity
  MPOt<Tensor> Hqb = -im*commutator(H,qlist[0]);//list is zero-indexed
  MPOt<Tensor> HQb = Hqb;
  rightwards_currents[0] = -Hqb;

  //jleft is list index not site index (zero-indexing vs 1-index)
  for(int jleft = 1; jleft <= N-p; jleft++) {
    if(qlist[jleft].N() != N) { Error("current: qlist[j] wrong length");}
    Hqb = -im*commutator(H, qlist[jleft]);
    /* Hqb = rightwards + leftwards.
           = rightwards - (next rightwards)
       (next rightwards) = rightwards - Hqb
     */
    rightwards_currents.push_back((rightwards_currents[jleft - 1] - Hqb).orthogonalize());
    HQb += Hqb;
    HQb.orthogonalize();
  }

  //TODO HQb zero

  return rightwards_currents;
}

/* construct the identity-environments (partial trace
   environments) */
template <class Tensor>
std::vector<MPOt<Tensor>>
left_identity_environments(MPOt<Tensor> H)
{
  auto sites = H.sites();
  int N = H.N();
  std::vector<Tensor> left_id_environments;

  Tensor El = 0.5*H.A(1)*sites.op("Id",1);
  left_id_environments.push_back(El);
  
  //j is vector (0-) index, not site (1-) index
  for(int j = 1; j < N; j++){
    El = 0.5*(El*(H.A(j+1)*sites.op("Id", j+1)));
    left_id_environments.push_back(El);
  }
  return El;
}


template <class Tensor>
std::vector<MPOt<Tensor>>
right_identity_environments(MPOt<Tensor> H)
{
  auto sites = H.sites();
  std::vector<Tensor> right_id_environments;
  int N = H.A(N);
  Tensor Er = 0.5*H.A(N)*sites.op("Id",N);
  right_id_environments.push_back(Er);
  
  //j is vector (0-) index, not site (1-) index
  for(int j = 1; j < N; j++){
    Er = 0.5*(Er*(H.A(N-j-1)*sites.op("Id", N-j-1)));
    right_id_environments.push_back(Er);
  }

  Er.reverse();

  return Er;
}

/* Partial traces down to each contiguous p-site region, with factors
   of 1/2 . Represented as MPOs on whole system. */
template <class Tensor>
std::vector<MPOt<Tensor>>
psite_components(MPOt<Tensor> H, int p,
		 std::vector<Tensor> El,
		 std::vector<Tensor> Er)
{

  int N = H.N();
  std::vector<Tensor> cpts;

  cpts.push_back(H.A(1)*Er[0+1]);
  //jleft is a site index, not a vector index
  for(int jleft = 2; jleft <= N-p; jleft++){
    Tensor thiscpt = El[jleft];
    for(int k = jleft; k <= jleft + p; k++){ thiscpt *= H.A(k); }
    
    /* consider e.g. jleft = 1, p = 1: this is a single-site op on
       site 1 (leftmost site). We want the partial-trace up through
       site 2; that's Er[1] = Er[jleft + p - 1].
    */
    if(jleft != N-p) { thiscpt *= Er[jleft + p - 1 ]; }
    cpts.push_back(thiscpt);
  }

  cpts.push_back(El[N-2]*H.A(N));
  return cpts;
}

/* Take a tensor t that's an operator on sites [jleft, jright] and
   turn it into an operator on sites [kleft, kright] by tensoring on
   identities from SiteSet sites.

   Doesn't verify that t is what it claims to be. (Probably should...)

   jleft, jright are site indices, not vector indices.
 */
template <class Tensor>
Tensor
embed(Tensor t, SiteSet sites,
      int jleft, int jright,
      int kleft, int kright )
{
  assert(kleft <= jleft);
  assert(jright <= kright);
  for(int n = kleft; n < jright; n++)     { t *= sites.op("Id", n);}
  for(int n = jright+1; n <= kright; n++) { t *= sites.op("Id", n);}
  return t;
}
  
/* Energy densities of a p-local Hamiltonian as list of tensors. Tries
   to be reasonably symmetric, so more complicated than it could be. */

template <class Tensor>
std::vector<MPOt<Tensor>>
energy_density(MPOt<Tensor> H, int p)
{
  int N = H.N();
  auto El = left_identity_environments(H);
  auto Er = right_identity_environments(H);

  auto eps = psite_components(H, p, El, Er);
  for(int r = 1; r < p; r++) {
    auto hp = psite_components(H, r, El, Er);

    //j is vector index, not site index
    for(int j = 0; j <= N-r; j++) {
      for(int k = j; k <= k + p - r; k++) {
	//Take my r-site operator. How many energy density terms is it split over?
	int number_split = std::min(std::min(k+1, N-k), p-r+1);
	eps[j] += embed(hp[k]/number_split, k + 1, k+r, j+1, j+r);
      }
    } 
  }
}

#endif //#ifndef MPOKPM_UTIL
