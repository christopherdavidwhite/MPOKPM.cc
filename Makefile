# copied from Miles' example

LIBRARY_DIR=$(ITENSOR_DIR)

COMMIT=`git log | head -1 | cut -d ' ' -f 2 | cut -b 1-7`
DATE=`date +%Y-%m-%d`
VDIR=verification-$(DATE)-$(COMMIT)

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

#there's almost certainly a cleaner way to do this
verification: dos
	mkdir -p verification-$(DATE)-$(COMMIT)
	./dos -L 12 -N 32 -M 512 -f $(VDIR)/L12-N32-M512
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(VDIR)/L12-N32-M512 -o $(VDIR)/L12-N32-M512
	./dos -L 32 -N 32 -M 512 -f $(VDIR)/L32-N32-M512
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(VDIR)/L32-N32-M512 -o $(VDIR)/L32-N32-M512
	./dos -L 64 -N 32 -M 512 -f $(VDIR)/L64-N32-M512
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(VDIR)/L64-N32-M512 -o $(VDIR)/L64-N32-M512
	./dos -L 128 -N 32 -M 512 -f $(VDIR)/L128-N32-M512
	$(JULIA) ./analysis/post-hoc-verification.jl -i  $(VDIR)/L128-N32-M512 -o $(VDIR)/L128-N32-M512
	./dos -L 256 -N 32 -M 512 -f $(VDIR)/L256-N32-M512
	$(JULIA) ./analysis/post-hoc-verification.jl -i $(VDIR)/L256-N32-M512 -o $(VDIR)/L256-N32-M512
