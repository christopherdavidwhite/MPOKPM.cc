#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"
#include "util.h"

#ifndef MPOKPM_ALGEBRA
#define MPOKPM_ALGEBRA
using namespace itensor;

const auto Dangler = IndexType("Dangler");

/* This is fundamentally broken: nominally uses IQ tensors, but
   doesn't even try to do the right thing wrt any nontrivial
   conservation. */
IQMPO
oplus(IQMPO A, IQMPO B)
{
  auto out = A;
  auto Bp = B;

  /* slightly obtuse: if either doesn't have a dangler index, add a
     trivial one, then proceed on the assumption that it does. */
  if (!findtype(out.A(1), Dangler)) {
    //add dangler index to A
    auto i1 = IQIndex("iq", Index("i",1, Dangler, 0), QN(0));
    auto S1 = setElt(i1(1));
    out.setA(1, out.A(1)*S1);
  }

  if (!findtype(Bp.A(1), Dangler)) {
    //add dangler index to B
    auto i2 = IQIndex("iq", Index("i",1, Dangler, 0), QN(0));
    auto S2 = setElt(i2(1));
    Bp.setA(1, Bp.A(1)*S2);
  }

  IQIndex iA = findtype(out.A(1), Dangler).dag();
  IQIndex iB = findtype(Bp.A(1), Dangler).dag();
  IQIndex iout = IQIndex("iq", Index("iout", iA.m() + iB.m(), Dangler, 0), QN(0));

  IQTensor Acap = IQTensor(iout, iA);
  for(int l = 1; l <= iA.m(); l++) { Acap.set(iout(l), iA(l), 1.0); }
  
  IQTensor Bcap = IQTensor(iout, iB);
  for(int l = 1; l <= iB.m(); l++) { Bcap.set(iout(iA.m() + l), iB(l), 1.0); }
  out.setA(1,Acap*out.A(1));
  Bp.setA(1, Bcap*Bp.A(1));

  std::cout << std::flush;
  
  //add for {A, B}
  out.plusEq(Bp);
  return out;
}


/* This guy behaves differently from Miles's nmultMPO in a
   user-relevant way in how it does truncation. In particular:
 
   - the denmatDecomp uses a hard-coded 1e-14 cutoff and no maxM
   - the Args args goes to the orthogonalize at the end.

   I really don't understand why Miles was doing things the way he
   did. Maybe he wanted something like zip-up---in this case, sure,
   apply one Args to the denmatDecomp, but you want to pass in another
   Args to orthogonalize. Don't make the user do that.

*/
template<class Tensor>
void 
nmultMPAlgebra(MPOt<Tensor> const& Aorig, 
	       MPOt<Tensor>const& Borig,
	       Tensor const& lbc, //left boundary condition
	       MPOt<Tensor>& res,
	       Args args)
{
  if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

  if(Aorig.N() != Borig.N()) Error("nmultMPO(MPOType): Mismatched N");
  const int N = Borig.N();

  auto A = Aorig;
  A.position(1);

  MPOt<Tensor> B;
  if(&Borig == &Aorig) { B = A; }
  else                 { B = Borig; B.position(1); }

  B.primeall();

  res=A;
  res.primelinks(0,4);
  res.mapprime(1,2,Site);

      
  Tensor clust,nfork;
  for(int i = N; i > 1; --i)
    {
      if(i == N) { clust = A.A(i) * B.A(i); }
      else       { clust = A.A(i) * (B.A(i) * nfork); }

      nfork = Tensor(leftLinkInd(A,i),leftLinkInd(B,i),leftLinkInd(res,i).dag());

      Tensor C;
      denmatDecomp(clust,nfork,C ,Fromright,{"Cutoff", 1e-14});
      res.setA(i, C);
    }
  clust = (lbc * A.A(1))*(B.A(1)*nfork);
  res.setA(1, clust);
  Tensor theta = clust * res.A(2);

  // res.A(1) needs
  //
  //  - common physical indices of A.A(1), B.A(1). Get these from the
  //    siteset.
  //
  //  - left dangler index of lbc (*not* either of the dangler indices
  //    of A,B). This is the dangler of clust.
  //
  //  - right link index for site 1 (bond 1)
  //
  // (This is so svdBond will group indices in the right way.)
  
  res.svdBond(1,theta,Fromleft);
  res.noprimelink();
  res.mapprime(2,1,Site);
  res.orthogonalize(args);

}//void nmultMPO(const MPOType& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)

template<class Tensor>
void
fitmultMPAlgebra(MPOt<Tensor> const& psi,
		 MPOt<Tensor> const& origK,
		 Tensor const& lbc, //left boundary condition
		 MPOt<Tensor>& res,
		 Sweeps const& sweeps,
		 Args args)
{
  auto N = psi.N();
  auto verbose = args.getBool("Verbose",false);
  auto normalize = args.getBool("Normalize",false);
  
  auto origPsi = psi;
  
  if(not res) res = origPsi;
  res.mapprime(1,2, Site);
  res.position(1);
  
  auto K = origK;
  
  K.primeall();
  origPsi.setA(1, origPsi.A(1)*lbc);

  /* Need to make the dangler indices match up. This means sticking
     something like lbc onto res. Can't just be lbc, because that'd
     leave an extra dangler index (the one that's supposed to hook up
     to K) dangling---we need to tensor on a vector, call it vec,
     such that lbc*vec is full-rank. Right way to do this
     *generically* is to choose vec random, but I know things about
     lbc and I'm lazy.*/

  //TODO assert that there's only one
  auto i = findtype(K.A(1), Dangler);
  auto vec = Tensor(i);
  vec.set(i(1), 1);
  res.setA(1, res.A(1)*lbc*vec);
  
  auto BK = std::vector<Tensor>(N+2);
  BK.at(N) = origPsi.A(N)*K.A(N)*dag(res.A(N));
  for(auto n = N-1; n > 2; --n) { BK.at(n) = BK.at(n+1)*origPsi.A(n)*K.A(n)*dag(res.A(n)); }
  
  for(auto sw : range1(sweeps.nsweep())) {
    args.add("Sweep",sw);
    args.add("Cutoff",sweeps.cutoff(sw));
    args.add("Minm",sweeps.minm(sw));
    args.add("Maxm",sweeps.maxm(sw));
    args.add("Noise",sweeps.noise(sw));
    
    for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
      {
	  
	if(verbose) { printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1); }
	auto lwfK = (BK.at(b-1) ? BK.at(b-1)*origPsi.A(b) : origPsi.A(b));
	lwfK *= K.A(b);
	auto rwfK = (BK.at(b+2) ? BK.at(b+2)*origPsi.A(b+1) : origPsi.A(b+1));
	rwfK *= K.A(b+1);
	  
	auto wfK = lwfK*rwfK;
	  
	if(normalize) wfK /= norm(wfK);
	//auto spec = res.svdBond(b,wfK,(ha==1?Fromleft:Fromright),PH,args);
	res.svdBond(b,wfK,(ha==1?Fromleft:Fromright),args);
	    
	if(ha == 1) { BK.at(b) = lwfK * dag(res.A(b)); }
	else { BK.at(b+1) = rwfK * dag(res.A(b+1)); }
      }
  }
  res.mapprime(2,1, Site);
}


template<class Tensor>
void
fitmultMPAlgebra( MPOt<Tensor> const& A,
		  MPOt<Tensor> const& B,
		  Tensor const& lbc,
		  MPOt<Tensor>& res,
		  Args const& args,
int nsweep )
{
  Sweeps sweeps(nsweep);
  auto cutoff = args.getReal("Cutoff",-1);
  if(cutoff >= 0) sweeps.cutoff() = cutoff;
  auto maxm = args.getInt("Maxm",-1);
  if(maxm >= 1) sweeps.maxm() = maxm;
  fitmultMPAlgebra(A,B,lbc, res,sweeps,args);
}


/*
  Compute tr[ Tn t1 Tn^{\dag} t2]

  Returns a guy with two danngler indices.

  Assumes orthogonality center left of jleft + len(tens).
*/
template<class Tensor>
Tensor
correlation(MPOt<Tensor> const& Tn,
	    int j1, int j2,
	    Tensor t1, Tensor t2)
{
  auto sites = Tn.sites(); 
  
  assert(Tn.orthoCenter() <= std::min(j1,j2) + 1);
  int jleft  = std::min(j1,j2);
  int jright = std::max(j1,j2);

  Tensor C;
  
  Tensor TAdag;
  Tensor tf, ts;

  //rightmost site
  if (j1 == j2)    { tf = t1;                    ts = t2; }
  else if(j1 < j2) { tf = sites.op("Id",jright); ts = t2; }
  else if(j2 < j1) { tf = t1;                    ts = sites.op("Id",jright); }

  TAdag = dag(Tn.A(jright));
  TAdag = swapPrime(TAdag, 0,1,Site);
  if(jright == 1) { TAdag.prime(Dangler); }
  else            { TAdag.prime(leftLinkInd(Tn,jright)); }
  
  C = Tn.A(jright);
  C *= prime(tf);
  C.mapprime(2,1,Site);
  C *= prime(TAdag, Site);
  C.mapprime(3,1, Site);
  C *= mapprime(dag(ts), 1,2);

  if (j1 == j2) { return C; }
  //between the two sites
  if(std::abs(j1 - j2) > 1) {
    for(int k = jright-1; k > jleft; k--){
      C *= Tn.A(k);
      C *= prime(dag(Tn.A(k)), Link);
    }
  }

  if(j1 < j2)      { tf = t1;                   ts = sites.op("Id",jleft); }
  else if(j2 < j1) { tf = sites.op("Id",jleft); ts = t2; }

  TAdag = dag(Tn.A(jleft));
  TAdag = swapPrime(TAdag, 0,1,Site);
  if(jleft == 1) { TAdag.prime(Dangler); }
  else           { TAdag.prime(leftLinkInd(Tn,jleft)); }

  C *= Tn.A(jleft);
  C *= prime(tf);
  C.mapprime(2,1,Site);
  C *= prime(TAdag);
  C.mapprime(3,1, Site);
  C *= mapprime(dag(ts), 1,2);
  C.mapprime(2,1);
  assert(2==C.r());
  return C ;
}

template<class Tensor>
std::vector<std::vector<Tensor>>
twopoint_correlation(MPOt<Tensor> const& Tn,
		     std::string q1,
		     std::string q2)
{
  assert(Tn.orthoCenter() == 1);
  int N = Tn.N();
  
  auto sites = Tn.sites();
  Tensor E; 
  std::vector<std::vector<Tensor>> correlations(N, std::vector<Tensor>(N));

  for(int jleft = 1; jleft <= N; jleft++) {
    for(int jright = jleft; jright <= N; jright++) {
      /*does twice as much as it has to for jleft=jright*/
      
      int j1, j2;
      Tensor C;
      j1 = jleft; j2 = jright;
      C = correlation(Tn, j1, j2, sites.op(q1,j1), sites.op(q2,j2));
      if(jleft != 1) { C *= E; }
      correlations[j1-1][j2-1] = C;

      j1 = jright; j2 = jleft;
      C = correlation(Tn, j1, j2, sites.op(q1,j1), sites.op(q2,j2));
      if(jleft != 1) {C *= E;}
      correlations[j1-1][j2-1] = C;
    }
    auto TAdag = dag(Tn.A(jleft));
    TAdag.prime(Dangler);
    TAdag.prime(Link);
    if (1 == jleft) { E  = Tn.A(jleft)*TAdag; }
    //suboptimal contraction ordering
    else            { E *= Tn.A(jleft)*TAdag; }
  } //jleft
  
  return correlations;
}
		     



#endif //#ifndef MPOKPM_ALGEBRA
