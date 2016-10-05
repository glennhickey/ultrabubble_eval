# this Makefile is largely derived from https://github.com/adamnovak/corg/blob/master/Makefile
.PHONY: all clean

# point this to wherever vg has been installed
VGDIR=../vg

VGLIBDIR=$(VGDIR)/lib
CXX=g++
INCLUDES=-I$(VGDIR)/ -I$(VGDIR)/include
CXXFLAGS:=-O3 -msse4.1 -fopenmp -std=c++11 -ggdb -g $(INCLUDES)
LDSEARCH=-L$(VGDIR)/src -L$(VGLIBDIR)
LDFLAGS=-lm -lpthread -lz -lbz2 -lsnappy -ldivsufsort -ldivsufsort64 -ljansson $(LDSEARCH)
LIBVG=$(VGLIBDIR)/libvg.a
LIBXG=$(VGLIBDIR)/libxg.a
LIBPROTOBUF=$(VGLIBDIR)/libprotobuf.a
LIBSDSL=$(VGLIBDIR)/libsdsl.a
LIBGSSW=$(VGLIBDIR)/libgssw.a
LIBSNAPPY=$(VGLIBDIR)/libsnappy.a
LIBROCKSDB=$(VGLIBDIR)/librocksdb.a
LIBHTS=$(VGLIBDIR)/libhts.a
LIBGCSA2=$(VGLIBDIR)/libgcsa2.a
LIBVCFLIB=$(VGLIBDIR)/libvcflib.a
LIBRAPTOR=$(VGLIBDIR)/libraptor2.a
LIBGFAK=$(VGLIBDIR)/libgfakluge.a
LIBSUPBUB=$(VGLIBDIR)/libsupbub.a
LIBSONLIB=$(VGLIBDIR)/libsonlib.a
LIBPINCHESANDCACTI=$(VGLIBDIR)/libpinchesandcacti.a
LIB3EDGECONNECTED=$(VGLIBDIR)/lib3edgeconnected.a
VGLIBS=$(LIBVG) $(LIBXG) $(LIBVCFLIB) $(LIBGSSW) $(LIBSNAPPY) $(LIBROCKSDB) $(LIBHTS) $(LIBGCSA2) $(LIBPROTOBUF) $(LIBRAPTOR) $(LIBGFAK) $(LIBSUPBUB) $(LIBSDSL) $(LIBSONLIB) $(LIBPINCHESANDCACTI) $(LIB3EDGECONNECTED)

#Some little adjustments to build on OSX
#(tested with gcc4.9 and jansson installed from MacPorts)
SYS=$(shell uname -s)
ifeq (${SYS},Darwin)
	LDFLAGS:=$(LDFLAGS) -L/opt/local/lib/ # needed for macports jansson
else
	LDFLAGS:=$(LDFLAGS) -lrt
endif

all: ub_eval

ub_eval.o: ub_eval.cpp
	$(CXX) ub_eval.cpp -c $(CXXFLAGS)

ub_eval: ub_eval.o $(VGLIBS)
	$(CXX) ub_eval.o $(VGLIBS) -o ub_eval $(CXXFLAGS) $(LDFLAGS)

clean:
	rm -f ub_eval.o ub_eval
