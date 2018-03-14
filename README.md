# Introduction

Compute conductivity (soon other quantities...) using a variant of the
kernel polynomial method in which operators are represented by matrix
product operators.

Computing dynamical quantities of large quantum systems is hard. One
can exactly diagonalize the Hamiltonian, at which point one
essentially knows everything, but computation times scale like N^3
(where N is the dimension of the Hilbert space) and memory
requirements like N^2. The Kernel polynomial method (KPM) brings both
of these down by a factor of N, for computation times like N^2 and
memory requirements like N.

For single-particle systems this is enough. In that case N ~ L^d,
where L is the system size and d the spatial dimension, but for
many-body systems (e.g. 3 for a chunk of silicon), so you can get to
fairly large sizes. For many-body systems, on the other hand, one has
to be more clever. The Hilbert space dimension goes as N ~ 2^(L^d), so
even after the speedup KPM buys us we can't proceed straightforwardly.

Enter the matrix product operator. All you really need for KPM is to
be able to multiply and add operators; matrix product operators (MPOs)
provide a way to lossily compress operators, dramatically reducing
memory requirements and computation times. The downside is that after
every multiplication or addition one has to re-compress; this repeated
lossy compression may (or may not!) destroy the physics you're trying
to understand.

# Method sketch and code layout
We use Miles Stoudenmire's ITensor (http://itensor.org/) for MPO
representations/operations. (It's a pretty great library, from what
I've seen: I've used a quasi-proprietary competitor and written my
own, and on a few days acquaintance ITensor seems way, way better. If
Miles is ever looking for a job, go hire him.)

KPM rests on the fact that one can write all kinds of quantities of
interest in terms of coefficients like

    mu[n,m] = trace(T[n](H) j T[m](H) j)

where `T[n](H)` and `T[m](H)` are the nth and mth Chebyshev
polynomials of the given Hamiltonian H and j is some operator
(e.g. the current operator). The Chebyshev polynomials are defined by
the base case + recursion relation

    T[0](H) = identity operator
    T[1](H) = H
    T[n](H) = 2*H*T[n-1](H) - T[n-2](H) ;

Note that despite the lossy compression a single Chebyshev may gake
gigabytes of memory to store, so we don't want to have more than the
bare minimum lying around.

Functions of interest:

  - `single_mu` computes a single coefficient `mu[n,m]` for given
      Chebyshevs `T[n]`,`T['m]`
  - `advance_chebyshevs` implements one sets `T[n]`, `T[n-1]` to
    `T[n+1]`, `T[n]` respectively
  - `all_mu` is a driver function using `advance_chebyshevs` to walk
    through the Chebyshev polynomials and `single_mu` to compute the
    relevant coefficients at each step.
	
You'll see a parameter `Maxm` in various places; this controls how
vigorously the lossy compression works.
