# Some notes on runtime checks
## Requirements for runtime checks
 - Need to turn on and off via compile-time flag. Macro? Probably just #IF.

    template <class T>
	void 
	check(T val, T target, string name, fail = None, desc = None, what="each")

 - Need a "timestep" etc.
## Checking specific functions
What invariants do we have for the Chebyshev polynomials?
 - bandwidth of `Tn` < 1 (this requires diagonalization)
 - `T[n+1] - T[n-1] = 2*H*Tn`
 - `T[n]^2 + (1 - H^2) U[n-1]^2 = 1` (this requires us to track the `U[n]`s)
 - 
