UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CXX=clang++-18
else
	CXX=clang++-mp-18
endif
all: avoid lattice

CXXFLAGS=-g -O3 -Wall -Wextra -std=c++20 -Wshadow -I/opt/local/include/

lattice.o: lattice.hex version.h
avoid.o: avoid.hex version.h

lattice.hex: lattice.cc
	xxd  -i lattice.cc  >lattice.hex
avoid.hex: avoid.cc
	xxd  -i avoid.cc  >avoid.hex

avoid: avoid.o  graphics_avoid.o
	$(CXX) $(CXXFLAGS) avoid.o  graphics_avoid.o -o $@  -L/opt/local/lib -lcairo   -lX11 -lm

lattice: lattice.o graphics.o 
	$(CXX) $(CXXFLAGS) lattice.o graphics.o -o $@ -L/opt/local/lib -lcairo    -lX11 -lm

version.h:
	git log -n 1 --format=format:"#define GIT_COMMIT \"%H\"%n" HEAD > $@
lint:
	scan-build-mp-18  -o html make
lint2:
	cppcheck --enable=all  --suppress=missingIncludeSystem avoid.cc
	cppcheck --enable=all  --suppress=missingIncludeSystem lattice.cc
clean:
	$(RM)  *.o avoid lattice *log *.dat *synctex.gz *~ *docx *out *.blg *.fls *.hex version.h
