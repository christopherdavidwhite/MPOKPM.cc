LIBRARY_DIR=$(ITENSOR_DIR)

COMMIT=`git log | head -1 | cut -d ' ' -f 2 | cut -b 1-7`
DATE=`date +%Y-%m-%d`

ifndef JULIA
JULIA=julia
endif

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

ifdef DEBUG
$(info DEBUG)
LIBFLAGS=$(LIBGFLAGS)
CCFLAGS=$(CCGFLAGS)
endif

LIBFLAGS += -lhdf5 -lhdf5_cpp

.PHONY: all
all: construct-algebra conductivity dos twopoint-correlation fourier

.PHONY: test
test: test.o util.o
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) $(OBJFILES) -lpthread -lgtest
	./test


OBJFILES=util.o
HFILES=util.h

%.o: %.cc algebra.cc util.cc chebyshev.cc globals.h util.h
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

fourier: fourier.o $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) $(OBJFILES) -lpthread

conductivity: conductivity.o $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< -o $@ $(LIBFLAGS) $(OBJFILES) -lpthread

twopoint-correlation: twopoint-correlation.o $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< $(OBJFILES) -o $@ $(LIBFLAGS) -lpthread

dos: dos.o $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< $(OBJFILES) -o $@ $(LIBFLAGS) -lpthread

construct-algebra: construct-algebra.o $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< $(OBJFILES) -o $@ $(LIBFLAGS) -lpthread

scratch: scratch.o  $(OBJFILES) $(HFILES)
	$(CCCOM) $(CCFLAGS) $< $(OBJFILES) -o $@ $(LIBFLAGS) -lpthread

clean:
	rm -fr .debug_objs *.o *-g *.ps conductivity dos construct-algebra twopoint-correlation scratch test

#parameters for verification
vL=8
vN=4
vM=512
VDIR=verification-$(DATE)-$(COMMIT)
vfn = $(VDIR)/L$(vL)-N$(vN)-M$(vM)

#there's almost certainly a cleaner way to do this
#in particular: separate target for each model
.PHONY: verification
verification: verification-setup rfheis 2NJW

.PHONY: verification-setup
verification-setup:

	rm -rf $(VDIR)
	mkdir -p $(VDIR)

.PHONY: rfheis
rfheis:
	./construct-algebra -L $(vL) -N $(vN) -M $(vM) -f $(vfn).$@ -m $@
	./fourier --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn)$@.
	./twopoint-correlation --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@ Sz
	./conductivity --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@ --model $@
	./dos          --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(vfn).$@ -o $(vfn).$@ -m $@

.PHONY: 2NJW
2NJW:
	./construct-algebra -L $(vL) -N $(vN) -M $(vM) -f $(vfn).$@ -m r$@
	./fourier --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn)$@.
	./twopoint-correlation --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@Sz
	./conductivity --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@ --model $@
	./dos          --sites $(vfn).$@.sites --dangler $(vfn).$@.chMPA -o $(vfn).$@
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(vfn).$@ -o $(vfn).$@ -m $@
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
