#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <random>
#include <complex>
#include <gtest/gtest.h>
#include <tgmath.h>
#define TESTING true
 
#include "chebyshev.cc"
#include "util.h"


namespace {

// The fixture for testing class Foo.
  class ConductivityTest: public ::testing::Test {
  protected:

  };
}

// Check double-KPM coefficients mu on random paramagnet Note that
// this isn't a real double KPM, in the sense that it doesn't use the
// current operator for the random paramagnet---that operator is
// identically 0. It does, however, approximate the MBL limit.
TEST_F(ConductivityTest, RandomParamagnetMu) {

  int L = 8;
  SiteSet sites   = SpinHalf(L);
  std::vector<double>  hz;
  std::complex<double> im(0,1); //1i
  int s = 0;
  std::default_random_engine e{s};
  std::uniform_real_distribution<double> d(0.0,1.0);
  
  auto rp_ampo = AutoMPO(sites);
  for(int b = 1; b <= L; b++)
    {
      double thz = d(e);
      hz.push_back(thz);
      rp_ampo += 2*thz,"Sz",b; //Pauli sigma not S
    }
  
  auto random_paramagnet = IQMPO(rp_ampo);
  auto j_ampo = AutoMPO(sites);
  for(int b = 1; b < L; ++b)
    {
      auto coeff = im; //sigma
      j_ampo += coeff,"S+",b,"S-",b+1;
      j_ampo += -coeff,"S-",b,"S+",b+1;
    }

  auto j = IQMPO(j_ampo);
  auto I = eye(sites);

  std::string filename = "conductivity.test";
  double cutoff=1e-14;
  
  int Maxm   = 1000;
  int N      = 4;
  double ups = 1;

  double shzl2 = 0;     for (double& thz : hz) { shzl2 += pow(thz,2); }
  double shzlhzlp1 = 0; for (int n = 0; n < L; n++) { shzlhzlp1 += hz[n]*hz[n+1]; }

  std::cout << "shzl2 " << shzl2 << "\n";
  std::cout << "shzlhzlp1 " << shzlhzlp1 << "\n";
  auto cheb = listandwrite_dangler(random_paramagnet, filename, N, Maxm, cutoff, true, 4, 16);
  auto mu = double_mu(cheb, cheb, j);
  
  ASSERT_EQ(mu.r(), 2);
  
  IQIndex i1 = mu.inds()[0];
  IQIndex i2 = mu.inds()[1];

  double mu00 = pow(ups,2) * pow(2, -1) * (L-1);
  double mu02 = pow(ups,2) * pow(2, -1) * ((L-1) * (2*shzl2 - 1) - 4*shzlhzlp1);

  //note that mu00 has the 2^-L I put in the conductivity calculation
  EXPECT_NEAR(double_mu(I,I,j).real(), pow(2, L)* mu00 , 1e-10);
  EXPECT_NEAR(real(mu.cplx(i1(0+1), i2(0+1))), mu00, 1e-10);
  EXPECT_NEAR(imag(mu.cplx(i1(0+1),i2(0+1))), 0, 1e-10);
  
  EXPECT_NEAR(real(mu.cplx(i1(0+1),i2(1+1))), 0, 1e-10);
  EXPECT_NEAR(imag(mu.cplx(i1(0+1),i2(1+1))), 0, 1e-10);
  
  EXPECT_NEAR(real(mu.cplx(i1(1+1),i2(0+1))), 0, 1e-10);
  EXPECT_NEAR(imag(mu.cplx(i1(1+1),i2(0+1))), 0, 1e-10);
  
  EXPECT_NEAR(real(mu.cplx(i1(0+1),i2(2+1))), mu02, 1e-10);
  EXPECT_NEAR(imag(mu.cplx(i1(0+1),i2(2+1))), 0, 1e-10);

  EXPECT_NEAR(real(mu.cplx(i1(2+1),i2(0+1))), mu02, 1e-10);
  EXPECT_NEAR(imag(mu.cplx(i1(2+1),i2(0+1))), 0, 1e-10);
}


//check correlation function
TEST(twopoint_correlation, allI_II) {
  int L = 8;
  SiteSet sites = SpinHalf(L);
  auto I = eye(sites);
  I.position(1);

  auto mu = twopoint_correlation(I, "Id", "Id");
  for(int j1 = 1; j1 <= L; j1++) {
    for(int j2 = 1; j2 <= L; j2++) {
      EXPECT_NEAR(real(mu.cplx(j1,j2)), pow(2,L), 1e-10);
      EXPECT_NEAR(imag(mu.cplx(j1,j2)), 0       , 1e-10);
    }
  }
}
  
TEST(twopoint_correlation, allI_zz) {
  int L = 8;
  SiteSet sites = SpinHalf(L);
  auto I = eye(sites);
  I.position(1);

  auto mu = twopoint_correlation(I, "Sz", "Sz");
  for(int j1 = 1; j1 <= L; j1++) {
    for(int j2 = 1; j2 <= L; j2++) {
      int expctval = j1 == j2 ? (pow(2,L)/4) : 0;
      EXPECT_NEAR(real(mu.cplx(j1,j2)), expctval, 1e-10);
      EXPECT_NEAR(imag(mu.cplx(j1,j2)), 0       , 1e-10);
    }
  }
}


TEST(twopoint_correlation, zz_II) {
  int L = 8;
  SiteSet sites = SpinHalf(L);
  auto zz_ampo = AutoMPO(sites);
  int k1 = 2;
  int k2 = L-3;
    
  zz_ampo += 1.0,"Sz",k1;
  zz_ampo += 1.0,"Sz",k2;
  zz_ampo += 1.0,"Id",1;
  
  IQMPO zz = IQMPO(zz_ampo);
  zz.position(1);
  std::cout << "position done\n"<<std::flush;
  auto mu = twopoint_correlation(zz, "Sz", "Sz");
  std::cout << "twpoint_correlation done\n"<<std::flush;
  for(int j1 = 1; j1 <= L; j1++) {
    for(int j2 = 1; j2 <= L; j2++) {
      int expctval = 0;
      expctval += (j1 == j2) ? pow(2,L)*3/8 : 0;
      bool match = (j1 == k1 && j2  == k2) || (j1 == k2 && j2  == k1);
      expctval += match ? pow(2,L)/8 : 0 ;
      //std::cout << j1 << " " << j2 << " " << (j1 == j2) << " " << match << " " << expctval << " " << real(mu[j1][j2].cplx()) << "\n";
      EXPECT_NEAR(real(mu.cplx(j1,j2)), expctval, 1e-10);
      EXPECT_NEAR(imag(mu.cplx(j1,j2)), 0       , 1e-10);
    }
  }
}

TEST(IOTest, writetohdf5) {
  int N = 3;
  
  /* create a tensor to write */
  Index ini = Index("i", N, Dangler, 0);
  Index inj = Index("j", N, Dangler, 0);
  Index ink = Index("k", N, Dangler, 0);

  auto tensor = ITensor(ini, inj, ink);

  /* populate this tensor in such a way that elements are unique and
     predictable */
  for(int i = 1; i <= N; i++){
    for(int j = 1; j <= N; j++){
      for(int k = 1; k <= N; k++) {
	std::cout << i << " " << j << " " << k << " " << (i-1)*pow(N, 2) + (j-1)*N + (k-1) << "\n";
	tensor.set(ini(i), inj(j), ink(k), (i-1)*pow(N, 2) + (j-1)*N + (k-1));
      }
    }
  }

  /* write the tensor */
  export_hdf5(tensor, "test-tensor.h5");
  tensor *= im;
  export_hdf5(tensor, "test-tensor.cplx.h5");

}

/* Julia
   include("analysis/post-hoc-verification.jl")
   H = rf_2NJW(zeros(L)); X,Y,Z,P,M = pauli_matrices_sparse(L); trace(Z[1]*H*Z[1]*H) * 2.0^(-L)/4
   etc.
*/
TEST(Verification, double_mu) {
  double hz = 0;
  int L = 4;
  SiteSet sites = SpinHalf(L);
  std::string filename = "/tmp/conductivity";
  
  std::vector<double> hzs;
  IQMPO H;
  double opnorm_bound;
  double cutoff = 1e-14;
  bool sweep_flag = false;
  int nsweeps = 0;
  int Maxm = 500;
  std::default_random_engine e{0};
  
  std::tie(H, hzs, opnorm_bound) = rf_2NJW(sites, hz, e);
  H *= 1/opnorm_bound;

  const double Sz1_H_Sz1_H = 0.00739644970414201;
  const double Sz1_H_Sz2_H = 0.002958579881656805;
  const double Sz1_Sz1_H   = 0;
  const double Sz1_HH      = 0;
  const double Sz1_Sz2_H   = 0.019230769230769232;
  //I am surprised that this is the same as the above
  const double HH          = 0.07692307692307693;
  /* different from Sz1_H_Sz1_H because commutes nontrivially through
     more terms: not PBC!*/
  const double Sz2_H_Sz2_H = 0.0014792899408284025;
  IQMPO Tn = listandwrite_dangler(H, filename, 0, Maxm, cutoff, sweep_flag, nsweeps, 32);
  IQTensor TnTn = double_mu(Tn, eye(sites), Tn, eye(sites));
  EXPECT_NEAR(TnTn.real(1,1), 1.0, 1e-10);
  EXPECT_NEAR(TnTn.real(1,2), 0.0, 1e-10);
  EXPECT_NEAR(TnTn.real(2,1), 0.0, 1e-10);
  EXPECT_NEAR(TnTn.real(2,2), HH,  1e-10);

  
  IQMPO Sz1 = singlesite_IQMPO("Sz", 1, sites);
  IQMPO Sz2 = singlesite_IQMPO("Sz", 2, sites);
  IQTensor mu11 = double_mu(Tn, Sz1, Tn, Sz1);
  EXPECT_NEAR(mu11.real(1,1), 0.25, 1e-10);
  EXPECT_NEAR(mu11.real(1,2), 0.00, 1e-10);
  EXPECT_NEAR(mu11.real(2,1), 0.00, 1e-10);
  EXPECT_NEAR(mu11.real(2,2), Sz1_H_Sz1_H, 1e-10);
    
  IQTensor mu12 = double_mu(Tn, Sz1, Tn, Sz2);
  EXPECT_NEAR(mu12.real(1,1), 0, 1e-10);
  EXPECT_NEAR(mu12.real(1,2), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mu12.real(2,1), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mu12.real(2,2), Sz1_H_Sz2_H, 1e-10);

  IQMPO alg = oplus(Sz1, Sz2, "cor_loc");
  //PrintData(alg.A(1));
  //PrintData(alg.A(2));
  IQTensor FSzSz = double_mu(alg, eye(sites), alg, eye(sites));
  EXPECT_NEAR(FSzSz.real(1,1), pow(2,L)/4, 1e-10);
  EXPECT_NEAR(FSzSz.real(1,2), 0, 1e-10);
  EXPECT_NEAR(FSzSz.real(2,1), 0, 1e-10);
  EXPECT_NEAR(FSzSz.real(2,2), pow(2,L)/4, 1e-10);


  IQMPO Sz1Sz2I = oplus( eye(sites), alg);
  IQTensor triple = double_mu(Sz1Sz2I, Sz1Sz2I, Sz1Sz2I, eye(sites));
  EXPECT_NEAR(triple.real(1,1,1), pow(2,L), 1e-10);
  EXPECT_NEAR(triple.real(2,2,1), pow(2,L)/4, 1e-10);
  EXPECT_NEAR(triple.real(2,1,2), pow(2,L)/4, 1e-10);
  EXPECT_NEAR(triple.real(1,2,2), pow(2,L)/4, 1e-10);
  EXPECT_NEAR(triple.real(3,3,1), pow(2,L)/4, 1e-10);
  EXPECT_NEAR(triple.real(3,1,3), pow(2,L)/4, 1e-10);
  
  EXPECT_NEAR(triple.real(1,1,3), 0, 1e-10);
  EXPECT_NEAR(triple.real(1,2,1), 0, 1e-10);
  EXPECT_NEAR(triple.real(3,3,3), 0, 1e-10);
  EXPECT_NEAR(triple.real(2,1,3), 0, 1e-10);

  IQTensor mualg = double_mu(alg, Tn, alg, Tn);
  EXPECT_NEAR(mualg.real(1,1,1,1), 0.25, 1e-10);
  EXPECT_NEAR(mualg.real(1,1,1,2), 0.0, 1e-10);
  EXPECT_NEAR(mualg.real(1,1,2,1), 0.0, 1e-10);
  EXPECT_NEAR(mualg.real(1,1,2,2), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mualg.real(1,2,1,1), 0.0, 1e-10);
  EXPECT_NEAR(mualg.real(1,2,1,2), Sz1_H_Sz1_H, 1e-10);
  EXPECT_NEAR(mualg.real(1,2,2,1), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mualg.real(1,2,2,2), Sz1_H_Sz2_H, 1e-10);

  EXPECT_NEAR(mualg.real(2,1,1,1), 0, 1e-10);
  EXPECT_NEAR(mualg.real(2,1,1,2), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mualg.real(2,1,2,1), 0.25, 1e-10);
  EXPECT_NEAR(mualg.real(2,1,2,2), 0.0, 1e-10);
  EXPECT_NEAR(mualg.real(2,2,1,1), Sz1_Sz2_H, 1e-10);
  EXPECT_NEAR(mualg.real(2,2,1,2), Sz1_H_Sz2_H, 1e-10);
  EXPECT_NEAR(mualg.real(2,2,2,1), 0, 1e-10);
  EXPECT_NEAR(mualg.real(2,2,2,2), Sz2_H_Sz2_H, 1e-10);
  
  /* (/ 0.0073964497041420097 0.0014792899408283937)
     this is really close to five (of all things)
  */
  
}
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
