#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"

#include <random>
#include <complex>
#include <gtest/gtest.h>
#include <tgmath.h>
#define TESTING true 
#include "conductivity.cc" //hack


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
    
  std::ofstream realmu_file("conductivity.test.re");
  std::ofstream imagmu_file("conductivity.test.im");
  
  int Maxm   = 1000;
  int N      = 4;
  double ups = 1;

  double shzl2 = 0;     for (double& thz : hz) { shzl2 += pow(thz,2); }
  double shzlhzlp1 = 0; for (int n = 0; n < L; n++) { shzlhzlp1 += hz[n]*hz[n+1]; }
  auto mu = all_mu(random_paramagnet, j, realmu_file, imagmu_file, N, Maxm, 1);
  double mu00 = pow(ups,2) * pow(2, L - 1) * (L-1);
  double mu02 = pow(ups,2) * pow(2, L-1) * ((L-1) * (2*shzl2 - 1) - 4*shzlhzlp1);
  
  EXPECT_NEAR(real(single_mu(I,I,j)), mu00 , 1e-10);
  EXPECT_NEAR(imag(single_mu(I,I,j)), 0, 1e-10);
  
  EXPECT_NEAR(real(mu[0][0]), mu00, 1e-10);
  EXPECT_NEAR(imag(mu[0][0]), 0, 1e-10);
  
  EXPECT_NEAR(real(mu[0][1]), 0, 1e-10);
  EXPECT_NEAR(imag(mu[0][1]), 0, 1e-10);
  
  EXPECT_NEAR(real(mu[1][0]), 0, 1e-10);
  EXPECT_NEAR(imag(mu[1][0]), 0, 1e-10);
  
  EXPECT_NEAR(real(mu[0][2]), mu02, 1e-10);
  EXPECT_NEAR(imag(mu[0][2]), 0, 1e-10);

  EXPECT_NEAR(real(mu[2][0]), mu02, 1e-10);
  EXPECT_NEAR(imag(mu[2][0]), 0, 1e-10);
}


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
