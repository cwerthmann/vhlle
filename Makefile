
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

CXX           = g++
CXXFLAGS      = -Wall -fopenmp -fPIC -O3
LD            = g++
LDFLAGS       = -O3

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoHadron.cpp eoAZH.cpp eoSmash.cpp trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp ic.cpp ic3F.cpp ickw.cpp icPartUrqmd.cpp icPartSMASH.cpp main.cpp rmn.cpp cornelius.cpp \
             icGlauber.cpp icGubser.cpp icGlissando.cpp icTrento.cpp multiHydro.cpp xsect.cpp EfI.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC))

TARGET	   = hlle_visc
#------------------------------------------------------------------------------
$(TARGET):       $(OBJS)
		$(LD)  $(LDFLAGS) $^ -o $@ -fopenmp $(LIBS)
		@echo "$@ done"
clean:
		@rm -f $(OBJS) $(TARGET)

$(OBJS): | $(objdir)

$(objdir):
	@mkdir -p $(objdir)

obj/%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
