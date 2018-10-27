# Status and usage
This is methods-research code. We don't have a solid understanding of
where the methods we implement will be reliable and useful, though we
have some notions; the whole point is to try that and see. If you try
to use this for cool physics, you'll be a guinea pig for the
methods. Which is great! We love methodological guinea pigs! But you
should be very aware of what you're signing up for.

If you do use the code (or just the methods it implements) and get
good results, let us know---we'd like to put out a methods paper with
you, in addition to your physics paper.

# Introduction

Compute conductivity, density of states, S(q,ω) using a variant of the
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
where L is the system size and d the spatial dimension (e.g. 3 for a
chunk of silicon), but for many-body systems, so you can get to fairly
large sizes. For many-body systems, on the other hand, one has to be
more clever. The Hilbert space dimension goes as N ~ 2^(L^d), so even
after the speedup KPM buys us we can't proceed straightforwardly.

Enter the matrix product operator. All you really need for KPM is to
be able to multiply and add operators; matrix product operators (MPOs)
provide a way to lossily compress operators, dramatically reducing
memory requirements and computation times. The downside is that after
every multiplication or addition one has to re-compress; this repeated
lossy compression may (or may not!) destroy the physics you're trying
to understand.

We do one better: we represent the whole space of powers of the
Hamiltonian in a single compact datastructure. Computing that is the
expensive part; we can write it to disk and use it at our leisure to
compute quantities of physical interest.

# Method sketch

KPM rests on the fact that one can write all kinds of quantities of
interest in terms of coefficients like

    μ[n,m] = trace(T[n](H) j T[m](H) j)

where `T[n](H)` and `T[m](H)` are the nth and mth Chebyshev
polynomials of the given Hamiltonian H and j is some operator
(e.g. the current operator). The Chebyshev polynomials are defined by
the base case + recursion relation

    T[0](H) = identity operator
    T[1](H) = H
    T[n](H) = 2*H*T[n-1](H) - T[n-2](H) ;

We compute the whole set of `T[n](H)`, they share a substantial amount of information, so it is better to store them

# Code structure 

We use Miles Stoudenmire's ITensor (http://itensor.org/) for MPO
representations/operations. (It's a pretty great library, from what
I've seen: I've used a quasi-proprietary competitor and written my
own, and ITensor is way, way better.)

`construct_algebra.cc` compiles to `construct_algebra`, which is the
core MPOKPM utility: it constructs the set of Chebyshev polynomials
`T[n](H)` and writes them to disk using in the ITensor native format.

A series of utilities read this structure and use them to compute quantities of physical interest. Those utilities are
 - `conductivity <-- conductivity.cc` computes conductivity KPM coefficients
   `tr( T[n](H) j T[m](H) j)`, where `j` is a current operator
 - `fourier <-- fourier.cc`, which computes
   `tr( T[n](H) Sz[q] T[m](H) Sz[q])`, where `Sz[q]` is the Fourier
   cosine transform of the onsite spins `Sz[j]`, `Sz[q] = Σ cos(πj/L)
   Sz[j]`.
 - `dos <-- dos.cc` computes the density-of-states KPM coefficients `tr( T[n](H) )`
 - `twopoint-correlation <-- twopoint-correlation.cc` computes all the two-point correlation
   functions `tr( T[n](H) Sz[j1] T[n](H) Sz[j2]` .
   
We test with `test <-- test.cc` and `post-hoc-verification.jl`. `test`
runs unit tests, while `post-hoc-verification.jl` compares with a
small-system reimplementation in Julia using exact diagonalization. If
you want to understand exactly what computations we're doing,
`post-hoc-verification` is a pretty good place to start.

# Testing
Do `make test` or `make verification`.
   
# Code printouts
If you want pdfs of the code for printing/reviewing in
Notability/etc., do `make ps`.

# Documentation for individual utilities
(at the moment you're going to have to either read the code or talk to
me (CDW) for this.)
## `conductivity`
### options
### example
## `construct_algebra`
### options
 - `--sweep` Use fit-multiply (sweeping multiplication
   algorithm). Don't do this.
 - `-e e` Error cutoff in multiplication
 - `-f FILE` Output to file `FILE`
 - `-h h` Width of disorder in onsite `Sz` fields 
 - `-Q Q` artificially rescale Hamiltonian by `Q`
 - `-s s` PRNG seed `s`
 - `-w, --nsweeps w` Use `w` sweeps if using fit-multiplication. Don't
   do this.
 - ` -L, --system-size L` System size
 - `-m, --model` Model to simulate; one of 
   + `rfheis`: random-field Heisenberg
   + `xx`: (isotropic) XY model
   + `rpara` : random paramagnets
   + `2NJW`: Heisenberg with second-neighbor Fermion hopping
 - `-M, --bond-dimension M` Maximum bond dimension (program exits when
   `T[0...n](H)` reaches this bond dimension)
 - `-N, --chebyshev-order N` Maximum Chebyshev order (program exits when
   it has computed `T[0..N]`)
### example
## `dos`
### options
### example
## `fourier`
### options
### example
## `twopoint-correlation`
### options
### example


# Dependencies

Core dependencies:
 - libpthread
 - libhdf5, libhdf5_cpp
 - git (the results of git are only used for filing the results of
   `make verification`, but at Present I expect every `make` command
   to chooke if it's not installed)

Additional dependencies for tests (`make test`):
 - Google Test
 
Additional dependencies for verification (`make verification`)
 - Julia < v0.7.0
 - Julia packages:
   + ArgParse
   + HDF5
   + PyCall
 - Python packages
   + h5py
