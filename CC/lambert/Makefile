CPPFLAGS=-I/opt/local/libexec/boost/1.81/include/
CPPFLAGS=-I$(HOME)/boost/boost_1_86_0/
CXXFLAGS=-std=c++23  -g -O3 -Wall  -Wextra -Wshadow # -Wno-deprecated
CXX=clang++-mp-19

all:lambert

lambert.o:version.h timer.h lambert.hex
	$(CXX)  $(CXXFLAGS) $(CPPFLAGS) -c lambert.cc -o $@
lambert:lambert.o
	$(CXX) $(CXXFLAGS) lambert.o -o $@

version.h:
	git log -n 1 --format=format:"#define GIT_COMMIT \"%H\"%n" HEAD > $@

lambert.hex: lambert.cc
	xxd  -i lambert.cc  > lambert.hex

clean:
	$(RM) -r *.o wilson lambert *.log *.hex  html/* *.dat source.code foo.html

lint:
	scan-build-mp-19  -o html make
lint2:
	cppcheck --enable=all  --suppress=missingIncludeSystem lambert.cc
print:
	 enscript -E --color  --fancy-header --language=html -fCourier10 -pfoo.html lambert.cc
tags:
	 /opt/local/bin/ctags -x --c++-kinds=pf --language-force=c++ lambert.cc  
doc:
	doxygen
	open html/index.html
