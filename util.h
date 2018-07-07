#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"
#include "globals.h"

using namespace itensor;

#define OPENE(fn, f) std::ofstream f(fn); if (!f) {error("Could not open file for writing");} f.precision(15);

MPOt<IQTensor> eye(SiteSet sites);

template <class Tensor> Tensor full(MPOt<Tensor> B);

std::tuple<IQMPO,std::vector<double>, double>
XX(SiteSet sites, double hz, std::default_random_engine e);

//construct our random-field Heisenberg Hamiltonian
//returns a tuple of IQMPO, fields, and upper bound on operator norm
std::tuple<IQMPO,std::vector<double>, double>
rfheis(SiteSet sites, double hz, std::default_random_engine e);

std::tuple<IQMPO,std::vector<double>, double>
rf_2NJW(SiteSet sites, double W, std::default_random_engine e);

std::tuple<IQMPO,std::vector<double>, double>
rpara(SiteSet sites, double hz, std::default_random_engine e);

template <class Tensor>
Spectrum
entanglement_spectrum(MPOt<Tensor> psi, int b);

template <class Tensor>
Tensor
single_mu(MPOt<Tensor> const& Tn,
	  MPOt<Tensor> const& opp );
  
template <class Tensor>
Tensor
double_mu(MPOt<Tensor> const& Tn,
	  MPOt<Tensor> const& Tm,
	  MPOt<Tensor> const& j);

template <class Tensor>
Tensor
double_mu(MPOt<Tensor> Tn,
	  MPOt<Tensor> j1,
	  MPOt<Tensor> Tm,
	  MPOt<Tensor> j2);

void write_singleKPM(IQTensor chtr,
		     std::ofstream& chtrre_file,
		     std::ofstream& chtrim_file);

std::vector<std::vector<std::complex<double>>>
write_doubleKPM(IQTensor mu,
		std::ofstream& realmu_file,
		std::ofstream& imagmu_file );

int export_hdf5(ITensor A, std::string filename, int flags);
int export_hdf5(ITensor A, std::string filename);
int export_hdf5(IQTensor A, std::string filename);
int export_hdf5(IQTensor A, std::string filename, int flags);

template <class Tensor>
MPOt<Tensor>
commutator(MPOt<Tensor> H,
	   MPOt<Tensor> q);

template <class Tensor>
std::vector<Tensor>
current(MPOt<Tensor> H,
	int p, //diameter of local conserved quantity
	std::vector<MPOt<Tensor>> qlist );

template <class Tensor>
std::vector<MPOt<Tensor>>
left_identity_environments(MPOt<Tensor> H);

template <class Tensor>
std::vector<MPOt<Tensor>>
right_identity_environments(MPOt<Tensor> H);

template <class Tensor>
std::vector<MPOt<Tensor>>
psite_components(MPOt<Tensor> H, int p,
		 std::vector<Tensor> El,
		 std::vector<Tensor> Er);

template <class Tensor>
std::vector<MPOt<Tensor>>
energy_density(MPOt<Tensor> H, int p);
