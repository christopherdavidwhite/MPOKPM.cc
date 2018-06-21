#include "itensor/mps/autompo.h"
#include "itensor/mps/sites/spinhalf.h"
#include "itensor/decomp.h"

#ifndef MPOKPM_ALGEBRA
#define MPOKPM_ALGEBRA
using namespace itensor;

const auto Dangler = IndexType("Dangler");
template <class Tensor>
MPOt<Tensor>
oplus(MPOt<Tensor> A, MPOt<Tensor> B)
{
  auto out = A;
  auto Bp = B;
  
  //dangler index
  auto i = IQIndex("iq", Index("i",2, Dangler, 0), QN(0));
  
  //add dangler index to A
  auto S1 = setElt(i(1));
  out.setA(1, out.A(1)*S1);
  
  //add dangler index to B
  auto S2 = setElt(i(2));
  Bp.setA(1, Bp.A(1)*S2);
  
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


/* L is a four-index guy: two dangler indices and two link indices
   (should be the left link index at site jleft).

   Assumes orthogonality center left of jleft + len(tens).

   Probably want tens to be std::vector<const &Tensor> or so: I find
   it aesthetically displeasing to be copying tensors all over the
   place like this. 
*/

template<class Tensor>
void
correlation(MPOt<Tensor> const& Tn,
	    Tensor const& E,
	    int jleft,
	    std::vector<Tensor> tens)
{
  int l = tens.size();
  int jend = jleft + l;
  assert(Tn.orthoCenter() <= j + l);

  //index linking i to i+1:
  auto ir = commonIndex(Tn.A(jend),Tn.A(jend+1),Link);
  
  auto C = Tn.A(jend)*tens(l-1)*dag(prime(Tn.A(jend),Site,ir));
  
  /* n counts leftwards. Already did n=0 just above. */
  for(int n = 1; n < l; n++){
    auto t = tens[l - n];
    C *= Tn.A(jleft + l -n);
    C *= t;
    C *= dag(prime(Tn.A(jleft + l - n),Link));
  }

  return E * C ;
}

template<class Tensor>
std::vector<std::vector<Tensor>>
twopoint_correlation(MPOt<Tensor> const& Tn,
		     std::vector<std::string> q1,
		     std::vector<std::string> q2)
{
  assert(Tn.orthoCenter() == 1);
  int N = Tn.N();
  
  auto sites = Tn.sites();
  int p1 = q1.size();
  int p2 = q2.size();

  std::vector<std::vector<Tensor>> correlations;

  /* slightly strange because jleft is a site, which is 1-indexed, but
     my output matrix is 0-indexed */

  auto i = findtype(Tn.A(1), Dangler);
  Tensor E = diag(prime(i), dag(i)) * prime(diag(prime(i),dag(i)));
  
  for(int jleft = 1; jleft <= N-fmax(p1,p2); jleft++) {
    for(int jright = 1; jright <= N-fmin(p1,p2); jright++) {
      for(int order = 0 ; order <= 1; order++) {
	std::vector<Tensor> tens;
	
	//for scoping
	int pleft, pright;
	std::vector<std::string> qleft, qright;
	if(0 == order) {
	  pleft  = p1;
	  pright = p2;
	  qleft  = q1;
	  qright = q2;
	} else if(1 == order) {
	  pleft  = p2;
	  pright = p1;
	  qleft  = q2;
	  qright = q1;
	} else {
	  Error("twopoint_correlation bad order");
	}

	/* if the two operators overlap, multiply them.
	   Some redundancy with non-overlap case, but I don't feel
	   compelled to be clever
	*/
	
	if(jright >= jleft + pleft - 1) {
	  
	  //bit that's just qleft
	  for(int k = jleft; k < jright; k++)
	    { tens.push_back(sites.Op(qleft[k-jleft], k)); }
	  
	  //overlap
	  for(int k = jright; k < jleft + pleft; k++) {
	    auto t = sites.Op(qleft[k-jleft], k);
	    t *= sites.Op(qright[k - jright], k);
	    t.mapprime(2, 1, Site);
	    tens.push_back(t);
	  }
	  
	  //bit that's just qright
	  for(int k = jright; k <= jright + pright - 1; k++) 
	    { tens.push_back(sites.Op(qright[k-jright], k)); }
	  
	} else {

	  //case where strings don't overlap
	  
	  for(int k = jleft; k < jleft + pleft; k++)
	    { tens.push_back(sites.op(qleft[k-jleft], k)); }
	  
	  for(int k = jleft + pleft; k < jright; k++)
	    { tens.push_back(sites.op("Id", k)); }
	  
	  for(int k = jright; k <= jright + pleft - 1; k++)
	    { tens.push_back(sites.op(qright[k-jright], k)); }
	} 	
	correlations[jright-1].push_back(correlation(Tn, E, jleft, tens));
      }//order
    } //jright
    //Update the left environment tensor
    E *= Tn.A(jleft)*prime(Tn.A(jleft));
  } //jleft
  
  return correlations;
}
		     



#endif //#ifndef MPOKPM_ALGEBRA
