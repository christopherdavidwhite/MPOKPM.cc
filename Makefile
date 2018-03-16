# copied from Miles' example

LIBRARY_DIR=/home/christopher/itensor

ifdef app
APP=$(app)
else
APP=conductivity
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

