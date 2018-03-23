// Miscellaneous utilities.

#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"

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


template <class Tensor>
double
operatornorm(MPOt<Tensor> const& B)
{
  //TODO
}

//construct our random-field Heisenberg Hamiltonian
std::tuple<IQMPO,std::vector<double>>
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
      H_ampo +=     "Sz",b,"Sz",b+1;
    }
  
  for(int b = 1; b <= L; ++b)
    {
      double hzj = hz * (2*d(e) - 1);
      fields.push_back(hzj);
      H_ampo += hzj, "Sz",b;
    }
  auto H = IQMPO(H_ampo);

  return std::make_tuple(H, fields);
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



#endif //#ifndef MPOKPM_UTIL
