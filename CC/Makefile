UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CXX=clang++-15
else
	CXX=clang++-mp-17
endif
all: avoid lattice

CXXFLAGS=-g -O3 -Wall -Wextra -std=c++20 -Wshadow -I/opt/local/include/


avoid: avoid.o
	$(CXX) $(CXXFLAGS) avoid.o -o $@
lattice: lattice.o graphics.o
	$(CXX) $(CXXFLAGS) lattice.o graphics.o -o $@ -L/opt/local/lib -lcairo  -L/opt/X11/lib  -lX11 -lm

lint:
	scan-build-mp-17  -o html make
lint2:
	cppcheck --enable=all  --suppress=missingIncludeSystem avoid.cc
	cppcheck --enable=all  --suppress=missingIncludeSystem lattice.cc
clean:
	$(RM)  *.o avoid lattice *log *.dat *synctex.gz *~ *docx *out *.blg *.fls
