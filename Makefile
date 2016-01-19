# Makefile for ATLAS searches

# Handle GZIP support.



PYTHIA8=/Users/nishita/Code/pythia82/
FASTJET3=/Users/nishita/Code/
PYTHIALIB=$(PYTHIA8)/lib
PYTHIAINC=$(PYTHIA8)/include/
FJINC=$(FASTJET3)/include
FJLIB=$(FASTJET3)/lib
OBJECTS=DecayChain.o
INCLUDE=-I$(PYTHIAINC) -I./ -I$(FJINC) -DGZIPSUPPORT 
LIB=-L$(PYTHIALIB) -lpythia8 -L$(FJLIB) -lfastjet -lz
CXXFLAGS=-g -O0

decayChain: $(OBJECTS)
	$(CXX) -o decayChain $(OBJECTS) -L$(PYTHIALIB) -lpythia8

%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $^ -I$(PYTHIAINC)
