#include <sys/time.h>
#include <sys/resource.h>
#include <string>

#include <fstream>

#ifndef _TIMER
#define _TIMER



class Timer{
 private:
  struct timeval t1,t2;
  rusage tr1, tr2;
 public:
  Timer(){
    gettimeofday( &t1, NULL);
    getrusage(RUSAGE_SELF,&tr1);
  }
  void stats(std::ofstream & log){
    gettimeofday( &t2, NULL);
    log << "Wall clock time \t"<<(t2.tv_sec-t1.tv_sec) + (t2.tv_usec-t1.tv_usec)/1000000. << std::endl<<std::flush;
    getrusage(RUSAGE_SELF,&tr2);
    std::cout<<"\u001B[32m"<<"CPU\t" << (tr2.ru_utime.tv_sec +  tr2.ru_utime.tv_usec/1000000.)<<"\t"<<
      tr2.ru_stime.tv_sec +  tr2.ru_utime.tv_usec/1000000.<<"\u001b[0m"<<std::endl;
  }
};

#endif
