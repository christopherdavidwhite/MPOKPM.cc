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
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)  -lpthread -lgtest

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g $(APP)-c

mkdebugdir:
	mkdir -p .debug_objs

#parameters for verification
vL=6
vN=8
vM=512
VDIR=verification-$(DATE)-$(COMMIT)
vfn = $(VDIR)/L$(vL)-N$(vN)-M$(vM)

#there's almost certainly a cleaner way to do this
verification: conductivity
	rm -rf $(VDIR)
	mkdir -p $(VDIR)
	./conductivity -L $(vL) -N $(vN) -M $(vM) -f $(vfn) -d
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(vfn) -o $(vfn)
	evince $(vfn)-plt-trTnerr.pdf &
	evince $(vfn)-plt-trTnjTnjerr.pdf &
