PROGRAM       = DrawingClasses
version       = TrackPairEfficiency
CXX           = g++
CXXFLAGS      = -O -Wall -Wno-bitwise-instead-of-logical -Wno-sizeof-pointer-div -fPIC -D$(version) # Add -g for debugger
LD            = g++
LDFLAGS       = -O
SOFLAGS       = -shared
#############################################
# -bind_at_load helps to remove linker error
############################################
CXXFLAGS += $(shell root-config --cflags)
LDFLAGS  = $(shell root-config --libs)

CXXFLAGS  += $(INCS)

HDRS	+= TrackPairEfficiencyHistogramManager.h TrackPairEfficiencyCard.h

SRCS = $(HDRS:.h=.cxx)
OBJS = $(HDRS:.h=.o)

all:            $(PROGRAM)

$(PROGRAM) :    $(OBJS)
		@echo "Linking $(PROGRAM) ..."
		$(CXX) -shared -lEG -lPhysics -L$(PWD) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $(PROGRAM).so
		@echo "done"

%.cxx:

%: %.cxx
#  commands to execute (built-in):
	$(LINK.cc) $^ $(CXXFLAGS) $(LOADLIBES) $(LDLIBS) -o $@

%.o: %.cxx %.h
#  commands to execute (built-in):
	$(COMPILE.cc) $(OUTPUT_OPTION) $<


clean:
		rm -f $(OBJS) $(PROGRAM).o $(PROGRAM).so

cl:  clean $(PROGRAM)
