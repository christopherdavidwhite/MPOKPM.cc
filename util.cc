// Miscellaneous utilities.

#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"

#ifndef MPOKPM_UTIL
#define MPOKPM_UTIL
using namespace itensor;

#define OPENE(fn, f) std::ofstream f(fn); if (!f) {error("Could not open file for writing");} f.precision(15);

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


template <class Tensor>
double
operatornorm(MPOt<Tensor> const& B)
{
  //TODO
}

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



void write_singleKPM(IQTensor chtr,
		     int N,
		     std::ofstream& chtrre_file,
		     std::ofstream& chtrim_file)
  
{
  IQIndex i1;

  //todo assert(1 == chtr.rank);
  //print the single-mu traces
  i1 = chtr.inds()[0];
  for(int n = 1; n <= N; n++){
    chtrre_file << real(chtr.cplx(i1(n))) << " ";  
    chtrim_file << imag(chtr.cplx(i1(n))) << " ";
  }
}

std::vector<std::vector<std::complex<double>>>
write_doubleKPM(IQTensor mu,
		int N,
		std::ofstream& realmu_file,
		std::ofstream& imagmu_file )
{ 
  //todo assert(2 == mu.rank);
  
  IQIndex i1;
  IQIndex i2;
  std::vector<std::vector<std::complex<double>>> vecmu;
  vecmu.reserve(N);
  //initialize mu to be 0
  for (int i = 0; i < N; i++) { vecmu.push_back(std::vector<std::complex<double>>(N,0)); }
  
  i1 = mu.inds()[0];
  i2 = mu.inds()[1];
  for(int n = 1; n <= N; n++){
    for(int m = 1; m <= N; m++){
      vecmu[n-1][m-1] = mu.cplx(i1(n),i2(m));
      realmu_file << real(mu.cplx(i1(n), i2(m))) << " ";  
      imagmu_file << imag(mu.cplx(i1(n), i2(m))) << " ";  
    }
    realmu_file << "\n";
    imagmu_file << "\n";
  }
  return vecmu;
}
#endif //#ifndef MPOKPM_UTIL
