
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)

CXX           = g++
CXXFLAGS      = -Wall -fPIC -g -march=native
LD            = g++
LDFLAGS       = -g -march=native

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)

vpath %.cpp src
objdir     = obj

SRC        = cll.cpp eos.cpp eo3.cpp eo1.cpp eoChiral.cpp eoHadron.cpp eoAZH.cpp trancoeff.cpp fld.cpp hdo.cpp s95p.cpp icurqmd.cpp icGlauber.cpp ickw.cpp icPartUrqmd.cpp main.cpp rmn.cpp cornelius.cpp jets.cpp
OBJS       = $(patsubst %.cpp,$(objdir)/%.o,$(SRC)) 
              
TARGET	   = hlle_visc.14.3
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
