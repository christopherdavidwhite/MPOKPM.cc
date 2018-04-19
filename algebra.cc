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
      else       { clust = nfork * A.A(i) * B.A(i); }

      nfork = Tensor(linkInd(A,i-1),linkInd(B,i-1),linkInd(res,i-1));

      denmatDecomp(clust,nfork,res.Anc(i) ,Fromright,{"Cutoff", 1e-14});

      auto mid = commonIndex(res.A(i),nfork,Link);
      mid.dag();
      res.Anc(i-1) = Tensor(mid,dag(res.sites()(i-1)),prime(res.sites()(i-1),2),rightLinkInd(res,i-1));
    }

  nfork = lbc * clust * A.A(1) * B.A(1);

  // res.A(1) needs
  //
  //  - common physical indices of A.A(1), B.A(1). Get these from the
  //    siteset.
  //
  //  - left dangler index of lbc (*not* either of the dangler indices
  //    of A,B). This is the dangler of nfork.
  // 
  //  - some kind of mid link index. Not entirely sure what to do
  //    here, but I'll copy Miles above.
  
  auto mid = commonIndex(res.A(1),nfork,Link);
  mid.dag();
  auto output_dangler = findtype(nfork, Dangler);
  res.Anc(1) = Tensor(output_dangler, mid,dag(res.sites()(1)),prime(res.sites()(1),2),rightLinkInd(res,1));
  
  res.svdBond(1,nfork,Fromleft);
  res.noprimelink();
  res.mapprime(2,1,Site);
  res.orthogonalize(args);

}//void nmultMPO(const MPOType& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)

#endif //#ifndef MPOKPM_ALGEBRA
