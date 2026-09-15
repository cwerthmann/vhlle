
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

CXX           = g++
# -march=native matches the optimised vHLLE build.  Drop it if you compile on a
# login node and run on compute nodes of a different microarchitecture -- on a
# heterogeneous cluster it produces binaries that die with SIGILL.
# -fno-math-errno lets sqrt() inline to a single instruction; it relaxes no
# IEEE rounding rule, unlike -ffast-math, which must NOT be used here (the
# Riemann solver and the primitive-variable recovery rely on well-defined
# behaviour near zero and on NaN comparisons).
CXXFLAGS      = -Wall -O3 -march=native -flto=auto -fno-math-errno -fopenmp -fPIC
LD            = g++
LDFLAGS       = -O3 -march=native -flto=auto -fno-math-errno -fopenmp

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoCMF.cpp eoCMFe.cpp eoHadron.cpp eoHadronPH.cpp eoAZH.cpp eoSmash.cpp eo2DTExS.cpp trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp ic.cpp ic3F.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp main.cpp rmn.cpp cornelius.cpp \
             icGlauber.cpp icGubser.cpp icGlissando.cpp icTrento.cpp multiHydro.cpp xsect.cpp EfI.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = hlle_visc
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
