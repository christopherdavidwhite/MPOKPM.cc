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

template<class MPOType>
void 
nmultMPAlgebra(MPOType const& Aorig, 
         MPOType const& Borig, 
         MPOType& res,
         Args args)
{
  using Tensor = typename MPOType::TensorT;

  if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

  if(Aorig.N() != Borig.N()) Error("nmultMPO(MPOType): Mismatched N");
  const int N = Borig.N();

  auto A = Aorig;
  A.position(1);

  MPOType B;
  if(&Borig == &Aorig)
    {
      B = A;
    }
  else
    {
      B = Borig;
      B.position(1);
    }

  B.primeall();

  res=A;
  res.primelinks(0,4);
  res.mapprime(1,2,Site);

  // need to give res.A(1) B's dangler index as well as A's so that
  // denmatDecomp (later) deduces the indices correctly
  auto Bd = findtype(B.A(1), Dangler);
  auto S = setElt(Bd(1));
  res.setA(1, S*res.A(1));
      
  Tensor clust,nfork;
  for(int i = 1; i < N; ++i)
    {
      if(i == 1) 
	{ 
	  clust = A.A(i) * B.A(i); 
	}
      else       
	{ 
	  clust = nfork * A.A(i) * B.A(i); 
	}

      if(i == N-1) break;

      nfork = Tensor(linkInd(A,i),linkInd(B,i),linkInd(res,i));

      denmatDecomp(clust,res.Anc(i),nfork,Fromleft,args);

      auto mid = commonIndex(res.A(i),nfork,Link);
      mid.dag();
      res.Anc(i+1) = Tensor(mid,dag(res.sites()(i+1)),prime(res.sites()(i+1),2),rightLinkInd(res,i+1));
    }

  nfork = clust * A.A(N) * B.A(N);

  res.svdBond(N-1,nfork,Fromright);
  res.noprimelink();
  res.mapprime(2,1,Site);
  res.orthogonalize();

}//void nmultMPO(const MPOType& Aorig, const IQMPO& Borig, IQMPO& res,Real cut, int maxm)


#endif //#ifndef MPOKPM_ALGEBRA
