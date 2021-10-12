PROGRAMS =	rigel
LD =            ld
OPT =           -O2
CXX =           g++
SYMBOLS =       -g
LDFLAGS =	-lm
CXXFLAGS =      $(SYMBOLS) $(GPROF) $(OPT) $(STLPORT_FLAGS) $(INC) $(GPROF) $(LIBS)
CFLAGS =      $(SYMBOLS) $(GPROF) $(OPT) $(INC) $(GPROF) $(LIBS)
MAKEFILE =      Makefile


NODE_OBJ = node.o
RIGEL_OBJ = node.o rigel.o
RIGEL_SRC = node.cc node.h \
	      rigel.cc rigel.h Makefile

all:            $(PROGRAMS)

rigel:	$(RIGEL_OBJ)
	$(CXX) -o $@ $(RIGEL_OBJ) $(LDFLAGS) $(LIBS) $(GPROF)

nn:	$(NODE_OBJ)
	$(CXX) -o $@ $(NODE_OBJ) $(LDFLAGS) $(LIBS) $(GPROF)

.cc.o:
	$(CXX) $(CXXFLAGS) -c $< -o $@

.c.o:
	$(C) $(CFLAGS) -c $< -o $@

clean:
	rm -f core* *~ *.o $(PROGRAMS) *.o
