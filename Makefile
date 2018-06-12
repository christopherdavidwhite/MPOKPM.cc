# copied from Miles' example

LIBRARY_DIR=$(ITENSOR_DIR)

COMMIT=`git log | head -1 | cut -d ' ' -f 2 | cut -b 1-7`
DATE=`date +%Y-%m-%d`

ifdef app
APP=$(app)
else
APP=conductivity
endif

ifndef JULIA
JULIA=julia
endif

CCFILES=$(APP).cc
#################################################################
#################################################################
#################################################################
#################################################################


include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/all.h

#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
#	$(CCCOM) -c $(CCGFLAGS) --define-macro CHECK=true -o $@ $<
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<
#Targets -----------------

build: $(APP)
debug: $(APP)-g
check: $(APP)-c

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS) -lpthread

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)  -lpthread 

clean:
	rm -fr .debug_objs *.o *.ps $(APP) $(APP)-g $(APP)-c

mkdebugdir:
	mkdir -p .debug_objs

#parameters for verification
vL=6
vN=30
vM=512
VDIR=verification-$(DATE)-$(COMMIT)
vfn = $(VDIR)/L$(vL)-N$(vN)-M$(vM)

#there's almost certainly a cleaner way to do this
#in particular: separate target for each model
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
