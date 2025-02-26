UNAME=$(shell uname -s)

ifeq ($(UNAME), Darwin)
        CXX=clang++-mp-19
        CPPFLAGS=-I/opt/X11/include -I/opt/local//include  -I$(HOME)/github/HighFive/include -I$(HOME)/boost/boost_1_87_0/  -I/opt/X11/include -I/opt/local//include
        LDFLAGS=-L/opt/local/lib -lhdf5  -lm
else
        CPPFLAGS=-I/usr/include/hdf5/serial -I$(HOME)/github/HighFive/include
        LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial   -L/usr/include/hdf5/serial  -lhdf5  -lm
        CXX=clang++-18
endif
CXXFLAGS= -Wall -Wshadow -Wextra  -std=c++23  -g -O3
CXX=clang++

OBJ=roots.o graphics_avoid.o
all: simul   # testi  testlj


prepare:
	makedepend *.cc *.h 2>/dev/null



roots.o: roots.hex version.h roots.hex

roots.hex: roots.cc
	xxd  -i roots.cc  >roots.hex

testlj: testlj.o ljforce.o
	$(CXX) $(CXXFLAGS) -o $@  testlj.o ljforce.o 
simul: $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@  $(OBJ) $(LDFLAGS) -lcairo  -L/opt/X11/lib  -lX11  

clean:
	$(RM) *.o simul */*.html *hex version.h source.code
clobber:clean
	$(RM) -rf *.dat simul t *.w *~ ./html/* *.log *.w
tar:
	tar cvfz ecmc.tar.gz *.cc *.h *.sh Makefile wmake *.m
	rsync -av ecmc.tar.gz tony@holbein.pct.espci.fr:
export:
	rsync -av  *cc *.h Makefile *.sh *.m wmake amaggs@oliver.pct.espci.fr:Psweep/
pull:clobberdata
	$(RM) files_w.tar.gz	
	/opt/local/bin/wget --cache=off https://turner.pct.espci.fr/~acm/files_w.tar.gz
	tar xfz files_w.tar.gz
	$(RM) files_w.tar.gz
web:
	scp files_w.tar.gz acm@turner:public_html/
webtrail:
	tar cvfz webtrail.tar.gz */*mat
	scp webtrail.tar.gz acm@turner:public_html/
gettrail:
	/opt/local/bin/wget https://turner.pct.espci.fr/~acm/webtrail.tar.gz
drop:
	cp -r *.cc *.h Makefile *.sh *.m wmake /Users/tony/Dropbox/LJW/
lint:
	 scan-build-mp-18 -o html make
check:
	cppcheck *.cc *.c
tidy:
	clang-tidy $(CXXFILES) *cc
print:
#	enscript -G2rE  --color  -whtml -phtml/code.html  *.h *.cc
	enscript -G2rE  --color  --fancy-header --language=html -fCourier10 -pfoo.html *.h *.cc
	open foo.html
testi:testi.o int.o
			$(CXX) $(CXXFLAGS) -o $@ testi.o int.o
clobberdata:
	#$(RM) -rf 1* 2* 3* 4* 5* 6* 7* 8* 9*
tags:
	/opt/local/bin/ctags -x --c++-kinds=pf --language-force=c++ roots.cc
version.h:
	git log -n 1 --format=format:"#define GIT_COMMIT \"%H\"%n" HEAD > $@
