/** @file
* **Factor field for Lennard-Jones.**
*/
#define DOGRAPH 0

#include <string>
#include <cassert>
#include <iostream>
#include <tuple>
#include <random>
#include <string>
#include <filesystem>
//#define BOOST_MATH_INSTRUMENT //debug
#include <boost/math/tools/roots.hpp>
#include  <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <fstream>
#include <list>
#include <string>

#include <highfive/H5File.hpp>
#include <highfive/H5Easy.hpp>
#include <highfive/H5DataSpace.hpp>

#include "timer.h"
#include "version.h"

#if DOGRAPH
#include "graphics_avoid.h"
#endif

namespace LJ{
  const double TAIL=1.; /*!< TAIL=tail amplitude in LJ, normally 1 */
  double Delta=(1+pow(2,1./6.))/2.; /*!< Delta: reset by integrate(), mean spacing in chain*/
  const double min_r = pow(2,1./6.); /*!< min_r position of LJ minimum */

  double p0 = 0.03; //fabs(-6*TAIL/pow(Delta,7) + 12/pow(Delta,13)); /*!< p0 pressure */
  double beta = 8.; /*!< \f[ \beta \f] inverse temperature*/

  double V0 = -1; /*!< V0: minimum of factor energy */
  double r0 = -1 ; /*!< r0: position of stationary point in potential*/
  double L=-1;
#if DOGRAPH
  Graphics *g;
#endif
}

namespace Control{
  int N=64;// chain size
  const int iter=128*8*8*2*8;
  double sample_time=128. ; //length of chain
  const int sweeps=10000 ;//for Sq-virial runs

  const bool load_samples=true;// for hot or cold runs
  const bool single_run=false;

  const bool Print=false; /*!< Control::Print detailed debugging information */
  bool DIST=false; /*!< DISR: print statistics on interactions/energy */
}

using namespace std;

ofstream dm ("debug.out");
ofstream em ("dist.out");
long int acounts[2]={0,0};
static string reset="\u001b[0m"; //color printing
static string yellow="\u001b[33m";
static string red="\u001b[31m";
static string green="\u001B[32m";


mt19937 gen(random_device{}()); //seed mersenne random numbers from quality source
uniform_real_distribution<double> uni;
exponential_distribution<> poisson(LJ::beta);
vector<float> samples;

void load_samples(){//load data generated by matlab program for distance between particles
  if (!filesystem::exists("delta.h5") ){
    cout<<red<<"No delta.h5\n"<<reset;
    exit(1);
  }
  H5Easy::File fd("delta.h5" , H5Easy::File::ReadOnly);
  samples = H5Easy::load< vector<float> > (fd,"/delta"); //whole energy array
  cout<<"first samples ";
  for(int i=0; i<10; i++) cout<<samples[i]<<" ";
  cout<<endl;
  double P =  H5Easy::load<double>(fd,"/P");
  double beta =  H5Easy::load<double>(fd,"/beta");
  cout<<"P "<<P<<" "<<beta<<endl;
  assert( fabs(P-LJ::p0) < 1.e-6 );
  assert( fabs(beta-LJ::beta) < 1.e-6 );
}

double
LJ_integrate(){/// find average pressure by numerical integration with weight \f[exp(-\beta(V(x)-px))\f]
  using namespace boost::math::quadrature;
  auto f1 = [](double t) { return exp(-LJ::beta* (1/pow(t,12) -LJ::TAIL/pow(t,6) +LJ::p0*t) ) ; }; /** Z */
  auto f2 = [](double t) { return t*exp(-LJ::beta* (1/pow(t,12) -LJ::TAIL/pow(t,6) +LJ::p0*t) );  };/// extension
  auto f3 = [](double t) { return (12./pow(t,12) -6.*LJ::TAIL/pow(t,6)) *exp(-LJ::beta* (1/pow(t,12) -LJ::TAIL/pow(t,6) +LJ::p0*t) );  };/// virial
  auto f4 = [](double t) { return (12./pow(t,13) -6.*LJ::TAIL/pow(t,7)) *exp(-LJ::beta* (1/pow(t,12) -LJ::TAIL/pow(t,6) +LJ::p0*t) );  };/// bond force

  double error1, error2, error3, error4;/**< errors */
  double Q1 = gauss_kronrod<double, 61>::integrate(f1, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error1);
  double Q2 = gauss_kronrod<double, 61>::integrate(f2, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error2);
  double Q3 = gauss_kronrod<double, 61>::integrate(f3, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error3);
  double Q4 = gauss_kronrod<double, 61>::integrate(f4, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error4);
  
  double l0=Q2/Q1;
  double W=Q3/Q1+1./LJ::beta;
  double av_force=Q4/Q1;
  LJ::Delta=l0;
  
  cout <<yellow<< "average virial= " <<W << " " <<reset<< endl;
  cout <<yellow<< "average force= " << av_force << " " <<reset<< endl;
  cout <<yellow<< "spring energy= " <<LJ::p0*l0 << " " <<reset<< endl;
  cout <<yellow<< "test " <<LJ::p0*l0 - W << " " <<reset<< endl;
  return l0;
}

inline
void prd_half(double & x){///< periodic boundary conditions
    if( x< -LJ::L/2. ) x+= LJ::L;
    if( x > LJ::L/2. ) x-= LJ::L;
}

inline
double LJ_energy(int l, int m, vector<double>& z){ ///< Potential energy
  double x = z[m]-z[l];
  prd_half(x);

  double x2= x * x;
  double x6 = 1./(x2* x2* x2);
  double x12= x6*x6;
  return (x12 - LJ::TAIL*x6) ;
}

/** LJ+linear force, derivatives for numerical root finder, called by boost halley solver, calculate function, first and second derivatives */
class rootLJ{
 private:
  double tail, force, E;
 public:
  /** Empty constructor */
  rootLJ   ( double _t, double _f, double _e ): tail(_t), force(_f), E(_e) {}
/** return function and derivates */
  tuple<double, double, double> operator()(double const& x){
    double x2=1./(x*x);
    double x6=x2*x2*x2;
    double x12=x6*x6;

    double fx= x12 - tail * x6 + force * x - E;
    double dx = (-12 * x12+ tail * 6 * x6)/x + force;
    double ddx = (12 * 13 * x12  - tail * 6 * 7 * x6) * x2;
    if(Control::Print) dm <<"Call "<< x<<" "<<fx<<" "<<dx<<" "<<ddx<<endl;

    return make_tuple(fx, dx, ddx);
  }
};

/** Find the minimum of potential */
double getstationary(){
//  cout<<green<<"Stationary "<<reset<<endl;
  double guess=1.;
  double min=.1;
  double max=10;
  int digits=40;
  unsigned long int it=100;
  rootLJ v(LJ::TAIL/2., -LJ::p0/12., 0.);
  double result = boost::math::tools::halley_iterate(  v , guess, min, max, digits, it);
//  cout<< green << "Stationary point\t" << result<<" "<< get<0>(v(result)) << reset << endl;
  return result;
}

/** Root finder using Halley iterates, implemented in boost. Delicate point is finding a good guess to start. */
double rightroot(double E, double x0, double force, double V0, double initial){  ///< root finder for LJ+constant force
  double min=.5;
  if(Control::Print) dm<<"x0 initial\t"<<x0<<" "<<initial<<endl;
  
  double mx=x0 < initial ? x0 : initial;
   assert(E-V0>=0);

  int digits=40;
  double stiff = (12*13/pow(LJ::min_r,14) - 6*7*LJ::TAIL/pow(LJ::min_r,8));
  double guess = x0 - sqrt((E-V0)/stiff);

  if(guess<min || E-V0 > 0.4){
    guess= pow(E+sqrt(E),-1./12.);
  }
  if(guess> mx){
     guess = (min+mx)/2.;
//     cout<<"yyy " <<guess <<" "<<min<<" "<<mx<<endl;
  }
  if(Control::Print) dm<<"right guess min max\t "<<guess<<"\t "<<min<<"\t "<<mx<<endl;
  
  unsigned long int it=1000;
  double result = boost::math::tools::halley_iterate(  rootLJ(LJ::TAIL, force, E), guess, min, mx, digits, it);
  if(Control::DIST) em <<"rit "<<it<<" "<<E-V0<<endl;

  if( initial - result <0 || it>15 ) {
    cout.precision(10);
    cout<<"right Iterations "<<it<<endl;
    rootLJ t(LJ::TAIL,force,E);
    tuple<double,double,double> tt=t(result);
    cout<<"rTuple "<< get<0>(tt)<<" "<<get<1>(tt)<<" "<<get<2>(tt)<<endl<<endl;
    cout<<"initial result:\t"<<initial <<" "<<result<<endl;

    tt=t(guess); 
    cout<<"guess rTuple "<< get<0>(tt)<<" "<<get<1>(tt)<<" "<<get<2>(tt)<<endl<<endl;
    
    sleep(1);
    cout.precision(6);
    abort();
  }

  return result;  
}


/** Root finder using Halley iterates, implemented in boost. Delicate point is finding a good guess to start. */
double leftroot(double E, double x0, double force, double V0, double initial){
  if(Control::Print) dm<<"left root "<<endl;
  
  //root finder for LJ+constant force
  double min=x0 > initial ? x0 : initial ;
  double mx=max(2*E/LJ::p0,5.);
  if(Control::Print) dm<<"initial "<<initial<<endl;
  int digits=40;
  if(E-V0<0){
    if(Control::DIST) em<<"E V0 "<<E<<" "<<V0<<" "<<E-V0<<endl;
    assert(E-V0>=0);
  }
  
  double guess=E/LJ::p0;
  if( E-V0 < 1.){
    double stiff = (12*13/pow(LJ::min_r,14) - 6*7*LJ::TAIL/pow(LJ::min_r,8))/2;
    guess = x0 + sqrt((E-V0)/stiff);
    if(Control::DIST) dm <<"left quand\n";
  }
  if(guess<min) guess = min*1.05;
  if(guess>mx) guess = mx*.95;

  if(Control::Print) dm<<"lguess min max\t "<<guess<<"\t "<<min<<"\t "<<mx<<endl;

  unsigned long int it=1000;
  double result = boost::math::tools::halley_iterate(  rootLJ(LJ::TAIL, force, E), guess, min, mx, digits, it);
  if(Control::DIST) em<<"lit "<<it<<" "<<E-V0<<" "<<guess<<endl;
  if(it>10 || result-initial <0 ){
    cout<<"Min "<<min<<endl;
    cout.precision(10);
    cout<<"res initial, R-I "<<result<<" "<<initial<<" "<<result-initial<<endl;
    cout<<"Guess G-x0\t"<<guess<<" "<<guess-x0<<endl;
    cout<<"Res -x0\t"<<result<<" "<<result-x0<<endl;
	
    cout<<"left Iterations "<<it<<endl;
    rootLJ t(LJ::TAIL,force,E);
    tuple<double,double,double> tt=t(result);
    cout<<"lTuple "<< get<0>(tt)<<" "<<get<1>(tt)<<" "<<get<2>(tt)<<endl;
    sleep(1);
    
    tt=t(guess);
    cout<<"guess lTuple "<< get<0>(tt)<<" "<<get<1>(tt)<<" "<<get<2>(tt)<<endl<<endl;
    cout<<endl;
    abort();
  }
  return result;  
}

/** The ECMC says go downhill to minimum then jump up */
double
n_rootlj_r(double V, double deltaE, int aa, int r, vector<double> & z) {//2-factor code combine LJ 
  double v;
  double x = z[aa] - z[r];
  prd_half(x);
  double x2 = x * x;
  if (x2 < LJ::r0*LJ::r0 ) {
    v = V + deltaE;
    if(Control::Print) dm<<"r1 "<<endl;
  } else {
    v = deltaE + LJ::V0;
    if(Control::Print) dm<<"r2 "<<endl;

  }
  if(Control::Print)  dm<<"right V "<<v<<endl;

  return rightroot(v, LJ::r0, LJ::p0,LJ::V0, fabs(x) );
}

/** The ECMC says go downhill to minimum then jump up */
double
n_rootlj_l(double V, double deltaE, int aa, int r,  vector<double> & z) {//2-factor code combine LJ 
  double v;
  double x = z[aa] - z[r];
  prd_half(x);
  double x2 = x * x;

  if (x2 > LJ::r0*LJ::r0 ) {
    v = V + deltaE;
    if(Control::Print) dm<<"l1 "<<endl;

  } else {
    v = deltaE +LJ::V0;
    if(Control::Print) dm<<"l2 "<<endl;

  }
  if(Control::Print)  dm<<"left V "<<v<<endl;
  return leftroot(v,LJ::r0,LJ::p0, LJ::V0, fabs(x));
}

/** Print out souce code. */
void LJ_code(){//print the code to this program
#include "roots.hex"
  ofstream opos ("source.code");
  for(auto i=0U;i<roots_cc_len;i++)
    opos<<roots_cc[i];
}

/** Getting the guesses write for the iterations
* is delicate, print out number of iterations/energy to see where things are slow
*/
void
test_LJ(){
  Control::DIST=false;
  vector<double> z(Control::N);
  
  cout<<yellow<<"p0 "<<LJ::p0<<endl;
  cout<<yellow<<"Delta "<<LJ::Delta<<reset<<endl;
  cout<<yellow<<"L\t"<<LJ::L<<reset<<endl;
  cout<<green<<"beta\t"<<LJ::beta<<reset<<endl;
  cout<<green<<"r0\t"<<LJ::r0<<reset<<endl;
  cout<<green<<"V0\t"<<LJ::V0<<reset<<endl;

  for(int l=0;l<50;l++){
    for(int i=0;i<Control::N;i++) z[i]=i*LJ::Delta +.05* uni(gen);

    int act=Control::N*uni(gen);
    if(Control::Print) dm<<"act "<<act<<endl;
    int right = (act + 1) % Control::N; //reducing separation
    int left = (act - 1 + Control::N) % Control::N; //increasing separation
    
    double x2l = z[left] - z[act];
    prd_half(x2l);
    double x2r = z[right] - z[act];
    prd_half(x2r);

    double Vleft=LJ_energy(left,act,z) + LJ::p0* abs(x2l);
    double Vright = LJ_energy(act, right,z) + LJ::p0* abs(x2r);
 
    double r1= -log(uni(gen)) / LJ::beta;
    double r2= -log(uni(gen)) / LJ::beta;
  
    double r_left = n_rootlj_l(Vleft,  r1 , left, act,z);
    double r_right = n_rootlj_r(Vright,  r2  , act, right ,z);
    if(Control::Print) dm<<"----------------------------"<<endl;
    rootLJ v(LJ::TAIL, LJ::p0, Vleft+r1);
    rootLJ v2(LJ::TAIL, LJ::p0, Vright+r2);
    if(Control::Print){
      auto a1=v(r_left);
      auto a2=v2(r_right);
      dm<<"r_left r_right\t"<< r_left<<" "<<r_right<<endl;
      dm<<"a1 a2\t"<< get<0>(a1)<<" "<<get<0>(a2)<<endl;
      dm<<endl;
    }
  }
  Control::DIST=false;
}

double
LJ_step(vector<double> & x, vector<int> & count, long int & a, long int &b, double & displacement){//single ECMC step "a" is active particle
//  cout<<"Step "<<endl;
  int n= x.size();
  double r_step, l_step;
  long int right = (a+1) % n;
  long int left= (a-1+n) % n;
  double sep_vec[4]={0,0};
  {//collision to right
    double sep = x[right] - x[a];
    sep = sep < 0 ? sep + LJ::L : sep;
    sep_vec[0]=sep;
    double Vright = LJ_energy(a, right,x) + LJ::p0 * abs(sep);
    double root= n_rootlj_r(Vright, poisson(gen) , a,  right, x);
    r_step = sep - root ;
    assert(r_step>=0);
  }
  {//collision left
    double sep=x[a] - x[left];
    if(sep < 0 ) sep += LJ::L;
    sep_vec[1]=-sep;

    double Vleft=LJ_energy(left,a,x) + LJ::p0 * abs(sep);
    double root= n_rootlj_l(Vleft, poisson(gen) ,  a, left, x);
    l_step =   root -sep;
    assert(l_step>0);
  }
  vector<double> candidate={ r_step, l_step};//two possible steps
  auto which = min_element( candidate.begin(), candidate.end() ) - candidate.begin();// find index of smallest candidate displacement
  x[a] += candidate[which];// move the active particle fowards
  displacement += sep_vec[which];
  x[a] =  x[a] > LJ::L ? x[a]- LJ::L : x[a]; //periodic
  if( which == 0 ){
    a = right;
    b++;
    acounts[0]++;
  }
  else if (which ==1){
    a = left;
    b--;
    acounts[1]++;
  }
  count[a]++;
  return candidate[which];
}


/** Write data to disk:  blist, and list are the most important for the curves*/
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


/** This calculates the full virial pressure, including perfect gas pressure term. */
pair<double ,double>
virial(vector<double> &x ){
  int n = x.size();
  double v = 0;
  double force=0;
  for( int i = 0; i<n; i++){
    int right = (i+1) % n;
    double sep= x[right] - x[i];
    sep = sep < 0 ? sep + LJ::L : sep;
    double x2= 1./(sep*sep);
    double x6=x2*x2*x2;
    double x12=x6*x6;
    double increment =  12*x12 - 6* LJ::TAIL *x6;
    v += increment;
    force += increment /fabs(sep);
  }
  force /= Control::N;
  v /= (Control::N);
  v=v+1./LJ::beta;
  v /= LJ::Delta;
  pair<double,double> par(v,force);
  return par;
}

/** fixed number of iterations of the update */
void
sweep(vector<double> & x, vector<int> &count, long int & a, long int &b, double & disp){ //N ECMC steps
  for (int i=0; i < 50 * Control::sample_time; i++){
    LJ_step(x, count, a, b, disp);
  }
  count[a]--;
}

/** Proper ecmc chain with imposed total displacement */
void
chain(vector<double> &x, vector<int> & count, long int & a, long int & b, double & displacement){
  double t_chain=0;
  while(true){
    long int olda=a;
    long int oldb=b;
    double s = LJ_step(x,count,a,b, displacement);
    if( t_chain + s > Control::sample_time ){// we have gone slightly too far, step back
        a=olda;
        b=oldb;
        count[olda]--;
        x[olda] -=  (-Control::sample_time + t_chain + s);
        displacement -=  (-Control::sample_time + t_chain + s);
        return;
    }
    else{
        t_chain += s;
    }
  }
}

/** Take the array of possible Delta and create new initial configuration. */
void
sample_hot(vector<double> &x){
  x[0]=1.e-9;
  LJ::L=0;
  for (auto i=1U;i<x.size();i++){ //fill up the array from matlab samples
    int which = samples.size() * uni(gen);
    LJ::L += samples[which];
    x[i] = x[i-1] + samples[which];
  }
  int which = samples.size() * uni(gen);
  LJ::L += samples[which];
  LJ::Delta = LJ::L/x.size();
//  cout<<green<<"reset L Delta"<<LJ::L<<" "<<LJ::Delta<<reset<<endl;
}

/** Loop over multiple chains to produce nu1, nu2.*/
void
lj_mc(vector<double> &  x, vector<int>  & master, list<int> & lst, list<int>  & blst, list<int>  & zlst, list<double>  & xlst, ofstream & log_out){
  ofstream pout ("drift.dat");
  ofstream vout ("virial.dat");

  for(int i=0;i<Control::N;i++)   x[i] = i * LJ::Delta; 

  for (int l=0; l<Control::iter; l++){ //main loop for running multiple simulations

    if( l % (1024*8) == 0) {
      cout<<"l= "<<l/1024<<" / "<<Control::iter/1024<<endl;
      log_out<<"l "<<l/1024<<" / "<<Control::iter/1024<<endl<<flush;
    }
    if( l % (1024*128) == 0)   save(x, master, lst, blst, zlst, xlst); //so we have something to analyse if we don't finish a long run

    if (Control::load_samples) sample_hot(x);
    else{//reset to ground state
      for(int i=0;i<Control::N;i++)  x[i] = i * LJ::Delta; 
    }

    long int a = 0;
    long int a0 = a; //remember origin
    long int b = 0; // "b" is almost "a" but unwrapped to follow activity
    vector<int> count(Control::N,0);//count site visits

    double displacement=0;
    chain( x , count, a, b, displacement); //chain() or replace by sweep()
    pout<<  Control::N * displacement/Control::sample_time/LJ::L/LJ::beta << endl;
    vout << virial(x).first << endl;

    lst.push_back( count[a] ); // H curve, rho_2
    blst.push_back( b ); //displacement, rho_1
    xlst.push_back( x[a] ); //final position
    zlst.push_back( count[a0] ); //origin visits

    for(int i=0; i<Control::N; i++){//global count curve, not analysed
      int pos = (a0+i) % Control::N;
      master[i] += count[pos];
    }
  }//Control::inter
}


/** This is the main, historic run that calculates \f(\rho_1\f) and \f(\rho_2\f)
* Here create the empty lists which are filled in the simulation with data
* and then call mc()
*/
void
run_distributions(ofstream & log_out, vector<double> &x){
  vector<int> count(Control::N);
  vector<int> master(Control::N);
  list<int> lst; //measures h, rho_2
  list<int> blst; //unwrapped final index, rho_1
  list<int> zlst; // returns to origin, local time
  list<double> xlst; //final position in x of chain

  lj_mc(x, master, lst, blst, zlst, xlst, log_out);
  save(x, master, lst, blst, zlst, xlst);
}

/** Use numerical integration to find mean separation, Delta*/
void init (){
  LJ_code();//prince source code to a file
  LJ::Delta = LJ_integrate();
  LJ::L = LJ::Delta * Control::N;
  LJ::r0 = getstationary();
  LJ::V0 = 1/pow(LJ::r0,12) - LJ::TAIL/pow(LJ::r0,6) + LJ::p0*LJ::r0;
  test_LJ();
  cerr<<"done test\n\n"<<endl;
}

double
entropy(vector<double> &x){
  double q0 = 2 * M_PI/LJ::L;
  double S=0;
  for (int i=1; i<Control::N; i++){
    double q=q0*i;
    double s=0;
    double c=0;
    for(double xi : x){
      double qx = q * xi;
      c += cos(qx);
      s += sin(qx);
    }//particle loop
    double sq= (c*c+s*s)/Control::N;
    //cout<<"Sq" <<sq<<endl;
    S+= (log(sq) - sq + 1);
  }//q-loop
  //cout<<endl<<endl;
  return S/LJ::L ;
}

/** Structure factor lowest mode*/
double Sq(vector<double> &x){
  double s = 0;
  double c = 0;
  double q = 2 * M_PI/LJ::L;
  for(double xi : x){
    double qx = q * xi;
    c += cos(qx);
    s += sin(qx);
  }
  double result = (c*c+s*s)/x.size();
  return result;
}

/** an array sent to the graphics */
void rectify2(vector<double> & z){
        double amin=*min_element(z.begin(), z.end());
        double amax=*max_element(z.begin(), z.end());
        for(int i=0;i<Control::N;i++) {
		z[i] = Control::N*(z[i]- amin)/(amax-amin);
        }
}
void rectify(vector<double> & z, double shift){
        for (auto &y:z) y -= shift/Control::N; 
        for(int i=0;i<Control::N;i++) {
		z[i]=  z[i] <0 ? z[i]+LJ::L : z[i];
		z[i]=  z[i] > LJ::L ? z[i]-LJ::L : z[i];
                z[i] -= LJ::Delta*i;
        }
}

/** Evolution of \f[ S(q) \f] and perhaps entropy during simulation, don't calculate the main distribution functions */
void
long_run(vector<double> & x){
  ofstream sq_out("sq.dat");
  ofstream vout ("virial.dat");
  ofstream dout ("displ.dat");
  ofstream ent ("ent.dat");
  long int a=0;
  long int b=0;
  vector<int> count(Control::N);

  if(!Control::load_samples){
  for(int i=0;i<Control::N;i++)
    x[i] = (i+.5) * LJ::Delta    ; //zero temp chain.
   }
   else{
	sample_hot(x);
   }

  for(int i=0;i<Control::sweeps;i++){
  double displacement=0;
#if DOGRAPH
	vector<double>z=x;
        rectify2(z);
	LJ::g->draw(z,  60 , a , a);
#endif
    chain(x, count,  a, b, displacement);
    ent <<entropy(x)<<endl;
    if( i> Control::sweeps/8){//thow away some data : not equilibrated
        dout << displacement/Control::sample_time<<endl;
        sq_out<<Sq(x)<<endl;
        vout<<virial(x).first<<endl;
    }
    if(i%1000==0) cout<<"i "<<i/1000<<"\t/"<<Control::sweeps/1000 <<endl;
  }
  ofstream xout ("pos.dat");
  for (double i : x) xout <<i<<"\n";
}

int
main(int argc, char * argv[]) {
  cout<<"Usage: simul N p0 beta\n";
  if (argc  > 1)    {
    Control::N = stoi(argv[1]);
    //  Control::sample_time=Control::N;
  }
  if (argc > 2) LJ::p0 = stod(argv[2]);
  if (argc > 3) LJ::beta = stod(argv[3]);
  ofstream log_out ("lj.log");
  cout<<"N "<<Control::N<<endl;
  cout<<"Beta p0"<<LJ::beta<<" "<<LJ::p0<<endl;
  log_out<<"N p0"<<Control::N<<" "<<LJ::p0<<endl;
  log_out<<"Beta "<<LJ::beta<<endl;
  vector<double> x(Control::N);
  log_out<<"git commit string "<<GIT_COMMIT<<endl;

  //samples
  if(Control::load_samples) load_samples();

  //reset variables
    init();

#if DOGRAPH
    double dmin=-.03*Control::N;
    double dmax=Control::N* 1.03 ;
    double diameter=0.2;
    LJ::g = new Graphics(Control::N, 1800, dmin, dmax, diameter);
    Control::sample_time /= 64;
#endif
  cout<<green<<" rho "<< 1./LJ::Delta<<reset<<endl;
  Timer t1;

  if(Control::single_run){
    cout<<"Long run"<<endl;
    log_out<<"Long run"<<endl;
    long_run(x);//structure factor+virial to disk from single chain
}
  else{
    cout<<"Run for rho1, rho2"<<endl;
    log_out<<"Run for rho1, rho2"<<endl;
    run_distributions(log_out, x);//rho_1 and rho_2
}
  double s=acounts[0] + acounts[1];
  cout<<"counts "<<acounts[0]/s<<" "<<acounts[1]/s<<endl;
  log_out<<"counts "<<acounts[0]/s<<" "<<acounts[1]/s<<endl;
  t1.stats(log_out);
  return 0;
}
