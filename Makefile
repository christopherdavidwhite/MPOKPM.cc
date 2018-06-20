LIBRARY_DIR=$(ITENSOR_DIR)

COMMIT=`git log | head -1 | cut -d ' ' -f 2 | cut -b 1-7`
DATE=`date +%Y-%m-%d`

ifndef JULIA
JULIA=julia
endif


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

.PHONY: all
all: construct-algebra conductivity dos


%.o: %.cc algebra.cc util.cc chebyshev.cc
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

conductivity: conductivity.o 
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) -lpthread

dos: dos.o 
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) -lpthread

construct-algebra: construct-algebra.o 
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) -lpthread

clean:
	rm -fr .debug_objs *.o *-g *.ps conductivity dos construct-algebra

verification: conductivity construct-algebra dos
	rm -rf $(VDIR)
	mkdir -p $(VDIR)
	./construct-algebra -L $(vL) -N $(vN) -M $(vM) -f $(vfn).rfheis -m rfheis
	./construct-algebra -L $(vL) -N $(vN) -M $(vM) -f $(vfn).2NJW  -m 2NJW
	./conductivity --sites $(vfn).rfheis.sites --dangler $(vfn).rfheis.chMPA -o $(vfn).rfheis
	./conductivity --sites $(vfn).2NJW.sites --dangler $(vfn).2NJW.chMPA -o $(vfn).2NJW #this will be wrong
	./dos          --sites $(vfn).rfheis.sites --dangler $(vfn).rfheis.chMPA -o $(vfn).rfheis
	./dos          --sites $(vfn).2NJW.sites --dangler $(vfn).2NJW.chMPA -o $(vfn).2NJW
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(vfn).rfheis -o $(vfn).rfheis -m rfheis
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(vfn).2NJW -o $(vfn).2NJW -m 2NJW

.PHONY: ps
ps: construct-algebra.cc.ps chebyshev.cc.ps  dos.cc.ps conductivity.cc.ps util.cc.ps algebra.cc.ps Makefile.ps

construct-algebra.cc.ps: construct-algebra.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

chebyshev.cc.ps: chebyshev.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

dos.cc.ps: dos.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

conductivity.cc.ps: conductivity.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

util.cc.ps: util.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

algebra.cc.ps: algebra.cc
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@

Makefile.ps: Makefile
	enscript --line-numbers -2 -E -q -Z -p - -f Courier6 $< > $@
