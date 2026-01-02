CXX=clang++-mp-21  
CXX=g++-mp-15
CXX=g++
CXXFLAGS=-g -Wall -Wshadow -Wextra -O3  -std=c++23 -I/opt/local/include/ -I$(HOME)/boost/boost_1_87_0/

all:bricks 

bricks:bricks.o graphics.o
	$(CXX) $(CXXFLAGS)  graphics.o bricks.o -o $@  -L/opt/local/lib -lcairo   -lX11 

single:single.o graphics.o
	$(CXX) $(CXXFLAGS)  graphics.o single.o -o $@  -L/opt/local/lib -lcairo   -lX11 


bricks.hex:bricks.cc
	xxd  -i bricks.cc  > bricks.hex

bricks.o:bricks.cc version.h bricks.hex

version.h:
	git log -n 1 --format=format:"#define GIT_COMMIT \"%H\"%n" HEAD > $@

clean:
	$(RM) *.o *hex bricks version.h *.dat 

lint:
	scan-build-mp-21  -o html make

lint2:
	cppcheck --enable=all  --suppress=missingIncludeSystem bricks.cc

print:
	enscript -E --color  --fancy-header --language=html -fCourier10 -pfoo.html bricks.cc
format:
	clang-format-mp-21 bricks.cc >j.cc
	mv -f bricks.cc bricks.cc~
	mv -f j.cc bricks.cc
