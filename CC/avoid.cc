#include <vector>
#include <iostream>
#include <random>
#include <cassert>
#include <fstream>
#include <cmath>
#include <list>
#include "timer.h"
#include "version.h"
#include "graphics_avoid.h"

using namespace std;
// energy (x_i - x_{i+1})^2/2 event drive mc on gaussian chain, to compare with true self-avoiding walk

namespace Param{
  int N = 1024; //number of particles
  int iter = 1024*128; //number of samples
  const int cold=true;  //hot or cold start
  const double T=1; // temperature
  const int GRAPH=0;//animation on screen
}

Graphics *g;

static string reset="\u001b[0m";
static string yellow="\u001b[33m";
static string red="\u001b[31m";

mt19937 gen(random_device{}()); //seed random number from quality source
uniform_real_distribution<double> uniform(0.,1.);
normal_distribution<double> gauss(0.0,1.0);


void
build_cholesky(vector<double> &d, vector<double> &o, vector<double> &c){ //cholesky factor for 1d laplacian
  int n= d.size();
  for(int i=0;i<n-1 ;i++) {
     d[i] = sqrt( (i+2.)/(i+1.));//diagonal
     o[i] = -1/d[i]; // diagonal+1
  }
  d[n-1]= sqrt(8.);//arbitrary to have +ve matrix
  for(int i=0;i<n;i++) c[i] = -1./sqrt( (i+1)*(i+2)); //correction for periodic bc
  c[n-2]= -sqrt( n/(n-1.));
  c[n-1]=d[n-1];
  o[n-2]= c[n-2];
}

void
random_config(vector<double> &conf, vector<double> & d,vector<double> & o,vector<double> & c, int test){ //solve for a configuration by back substitution
  int n=d.size();
  vector<double> r(n);
  if(!test) {
    for (int i=0;i<n;i++) r[i]=gauss(gen);
  }
  else{
    for (int i=0;i<n;i++) r[i]=sqrt(i+1);
  }

  conf[n-1]= r[n-1]/d[n-1];
  conf[n-2]= (r[n-2]-conf[n-1]*o[n-2] )/ d[n-2];
  for(int i=0;i<n-2;i++) conf[n-3-i]= (r[n-3-i]-conf[n-2-i]*o[n-3-i] - c[n-3-i]*conf[n-1])/ d[n-3-i];
}


void code(){//print the code to this program
#include "avoid.hex"
  for(auto i=0U;i<avoid_cc_len;i++)
    cout<<avoid_cc[i];
  exit(0);
}

double
step(vector<double> & x, vector<int> & count, int & a, int &b){//single ECMC step
  int n= x.size();
  double r_step, l_step, energy;
  int right = (a+1) % n;
  int left= (a-1+n) % n;
  
  //collision to the right
  double sep=x[right]-x[a];
  energy =  -Param::T*log( uniform(gen));
  if ( sep < 0 )
    energy += sep*sep/2.;
  r_step= sqrt(2*energy)+sep;

  //collision to the left
  sep=x[a] - x[left];
  energy = -Param::T*log( uniform(gen));
  if( sep > 0 )
    energy += sep*sep/2.;
  l_step =  sqrt( energy*2) -sep;

  
  if (r_step<l_step){ //now choose between two candidates
    x[a] += r_step;
    a=right;
    b++;
  }
  else{
    x[a] += l_step;
    a=left;
    b--;
  }
  count[a]++;
  return min(r_step, l_step);
}

void
sweep(vector<double> & x, vector<int> &count, int & a, int &b){ //N ECMC steps
  for (int i=0; i<Param::N;i++){
    step(x, count, a, b);
  }  
  count[a]--;
}

void
build_harmonic(vector<double> & x){
  int n=x.size();
  vector<double> d(n),o(n), c(n);
  build_cholesky(d, o, c);
  random_config(x, d,o, c,0);
}

void
save(vector<double> const&  x, vector<int> const & master, list<int> const & lst, list<int> const & blst, list<int> const & zlst, list<double> const & xlst){
  ofstream cnt ("count.dat"); //average L as a function of position
  ofstream opos ("pos.dat"); //final configuration
  ofstream lout ("list.dat"); //H for rho_2
  ofstream bout ("blist.dat");// unwrapped position of active pointer, for rho_1
  ofstream zout ("zlist.dat"); //returns to origin
  ofstream xout ("xlist.dat"); //x[a] for last point in simulation
  
  for (double i : x) opos<<i<<"\n";
  for (int i : master) cnt<<i<<"\n";
  for (int n : lst) lout << n << "\n";
  for (int n : blst) bout << n << "\n";
  for (int n : zlst) zout << n << "\n";
  for (double xx : xlst) xout << xx <<"\n";
}

void
start(vector<double> & x){
  int a=0;
  int b=0;
  vector<int> count(Param::N); 
  for (int i=0; i<Param::N ;i++){ //equilibration steps
    sweep( x , count, a, b);
  }
  cout<<"Warm start"<<endl;
}

void 
chain(vector<double> &x, vector<int> & count, int & a, int & b){
  double t_chain=0;
  double sample_time= Param::N;
  if(Param::GRAPH) sample_time *= 100;
  int olda,oldb;
  while(true){
    olda=a;
    oldb=b;
    double s=step(x,count,a,b);
    if( t_chain+s>sample_time ){// we have gone slightly too far, step back
        a=olda;
        b=oldb;
        count[olda]--;
        x[olda] -=  (-sample_time+ t_chain +s);
        return;
    }
    else{
        t_chain += s;
    }
   if(Param::GRAPH)
        g->draw(x,  60 , a , a);
  }
}

void
mc(vector<double> &  x, vector<int>  & master, list<int> & lst, list<int>  & blst, list<int>  & zlst, list<double>  & xlst, ofstream & log_out){
  uniform_int_distribution<> ran_int(0,Param::N-1); // distribution in range [0, N-1]
  vector<int> count(Param::N);
 if(Param::GRAPH) {
    double dmin=0.02*Param::N;
    double dmax=Param::N*1.05;
    double diameter=0.2;
    g= new Graphics(Param::N, 1800, dmin, dmax, diameter);
  }
  for (int l=0; l<Param::iter; l++){ //main loop for data 
    int a;
    if( l % (1024*8) == 0) {
      cout<<"l "<<l/1024<<endl;
      log_out<<"l "<<l/1024<<" of "<<Param::iter/1024<<endl<<flush;
    }
    if( l % (1024*128) == 0)   save(x, master, lst, blst, zlst, xlst);
    
    if (Param::cold) {
      std::fill(x.begin(), x.end(), 0); // reset position
      if (Param::GRAPH)
          for (int i=0;i<Param::N;i++) x[i] = 30*sin( (4*2.*i*M_PI)/Param::N)-20;
    }
    else{
      build_harmonic(x);
    }

    a = 0;
    int a0 = a; //remember origin
    int b = 0; // "b" is almost "a" but unwrapped to follow activity
    std::fill(count.begin(), count.end(), 0); // reset count of site visits    
    chain( x , count, a, b); //chain() or replace by sweep()
    lst.push_back( count[a] ); // H curve, rho_2
    blst.push_back( b ); //displacement, rho_1
    xlst.push_back( x[a] ); 
    zlst.push_back( count[a0] ); //origin visits
    
    for(int i=0; i<Param::N; i++){//global count curve, not analysed
      int pos = (a0+i)%Param::N;
      master[i] += count[pos];
    }
  }
}

void
test() {
  cout<<"Test of cholesky code, to compare with matlab solution for debugging"<<endl;
  vector<double> x(7);
  int n=x.size();
  vector<double> d(n),o(n), c(n);
  build_cholesky(d, o, c);
  random_config(x, d,o, c,1);
  cout<<setprecision(5);
  cout<<"d \t";
  for(auto y:d) cout<<y<<" ";
  cout<<"\no \t";
  for(auto y:o) cout<<y<<" ";
  cout<<"\nc \t";
  for(auto y:c) cout<<y<<" ";
  cout<<"\nx \t";
  for(auto y:x) cout<<y<<" ";
  cout<<endl;
  cout<<"this last line should be"<<endl;
  cout<<"  3.9533    5.5570    6.1357    5.8691    4.8665    3.2032    0.9354"<<endl;
}


int
main(int argc, char * argv[]){
  if( argc ==2 && argv[1] == string("code"))//print out this c++ code
    code();
  if (argc  > 1) 
    Param::N = stoi(argv[1]);
  if (argc > 2)
    Param::iter = 1024*stoi(argv[2]);
  test();
  ofstream log_out ("harmonic.log");
  log_out<<"git commit string "<<GIT_COMMIT<<endl;
  list<int> zlst; // returns to origin, local time
  list<int> lst; //measures h, rho_2
  list<int> blst; //unwrapped final index, rho_1
  list<double> xlst; //final position in x of chain
 
  vector<double> x(Param::N);
  vector<int> master(Param::N);
  Timer t1;
  
  cout<<yellow<<"N= "<<Param::N<<"Iterations="<<Param::iter<<"\tcold="<<Param::cold<<endl<<flush;
  log_out<<"N="<<Param::N<<"\tIterations="<<Param::iter<<"\tcold="<<Param::cold<<endl<<flush;

  mc(x, master, lst, blst, zlst, xlst, log_out);

  save(x, master, lst, blst, zlst, xlst);
  log_out<<"Done"<<endl;
  cout<<red<<"Done"<<reset<<endl;
  t1.stats(log_out);
  return 0;
}
