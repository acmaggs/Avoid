/** @file 
* **Lambert function for factor fields**
*
* Exponentially interacting particles in one dimension
*/
//#include "lambert.h"
long int acounts[4]={0,0,0,0}; //count events according to kind

namespace Lambert{
  const double a[2]={20., 40.}; //exponential prefactor first and second neighbours
  const double delta= -0.9067 ; //vary this to plot figure 5 , splitting FF between two neigbours
  double p = 0.9067  + delta; //1.0170 +delta;
  double p2 = 0.2528  -delta/2.;
  const double ell=2; //exponential decay length
  const double T=.3;
  double pbar, pbar2;
  double xmin, emin,xav, emin2, xmin2; //positions + values of minimum of potential including pressure
  double L; //system size
  static bool cold=true;  //cold start
  static double scale_factor=1;
  const bool second=true; //second nearest neighbour interactions
  const double BIG=1.e9;
  const bool single_run=false;  //calculate structure factor, or rho1, rho2 if false
}

namespace Control{
  const int N = 1024; //number of particles
  const int iter = 1024*64; //number of samples for rho1, rho2
  const int sweeps=1024*8*8; // for Sq/virial etc
  const double sample_time=N; //length of ecmc chain
}

#include <boost/math/special_functions/lambert_w.hpp> // For lambert_w function.
#include <boost/math/quadrature/gauss_kronrod.hpp> //numerical integration to find mean separation
#include <cmath>
#include <vector>
#include "timer.h"
#include "version.h"
#include <cassert>
#include <iostream>
#include <random>
#include <list>
#include <cmath>

using namespace std;
using boost::math::lambert_w0;
using boost::math::lambert_wm1;

static string reset="\u001b[0m"; //color printing
static string green="\u001B[32m"; 
static string yellow="\u001b[33m";
static string red="\u001b[31m";

mt19937_64 gen(random_device{}()); //seed mersenne random numbers from quality source
exponential_distribution<> poisson(1./Lambert::T);

/** Print source code to file */
void code(){
#include "lambert.hex"
  ofstream scout("source.code");
  for(unsigned int i=0; i<lambert_cc_len; i++)
    scout<<lambert_cc[i];
}

/** Structure factor lowest mode*/
double Sq(vector<double> &x){
  double s = 0;
  double c = 0;
  double q = 2 * M_PI/Lambert::L;
  for(double xi : x){
    double qx = q * xi;
    c += cos(qx);
    s += sin(qx);
  }
  double result = (c*c+s*s)/x.size();
  return result;
}

/** Solution of LambertW for first or second nearest neigbour */
inline double
root_right(double r,  int s){
  double pbar = s == 1 ? Lambert::pbar : Lambert::pbar2;
  double p = s == 1 ? Lambert::p : Lambert::p2;
  if( fabs(p) < 1.e-12 ) {
	double res= -Lambert::ell*log(r/Lambert::a[s-1]);
	return res;
}

  double rbar=r/Lambert::ell/p;
  double rhs = -exp( -rbar ) /pbar;
  assert(rhs > -exp(1.));
  double root = -lambert_wm1(rhs); //this is positive
  assert(root>=0);
  double x = -Lambert::ell*log(root*pbar);
  return x;
}

/** Solution of LambertW for first or second nearest neigbour */
inline double
root_left(double r,  int s){
  double pbar = s == 1 ? Lambert::pbar : Lambert::pbar2;
  double p = s == 1 ? Lambert::p : Lambert::p2;
  if(fabs(p) < 1.e-12) {return 1.e9; }

  double  rbar=r/Lambert::ell/p;
  double  rhs = -exp( -rbar)/pbar;
  assert(rhs > -exp(1.));
  double root = -lambert_w0(rhs);
  double x = -Lambert::ell * log(root*pbar);
  return x;
}

inline double res_energy2(double x, double r){
  double z = Lambert::a[1] * exp( -x / Lambert::ell) + Lambert::p2 * x - r;
  return z;
}

inline double res_energy(double x, double r){
  double z = Lambert::a[0] * exp( -x / Lambert::ell) + Lambert::p * x - r;
  return z;
}

/** Calculate virial and mean force in bonds force, nearest neighbour only */
pair<double ,double>
virial(vector<double> &x,  int s ){
    int n = x.size();
    double v = 0;
    double force=0;
    for( int i = 0; i<n; i++){
        int right = (i+s) % n;
	double sep= x[right] - x[i];
	sep = sep < 0 ? sep + Lambert::L : sep;
        v += sep * (Lambert::a[s-1]/Lambert::ell) * exp(-sep/Lambert::ell);
        force += (Lambert::a[s-1]/Lambert::ell) * exp(-sep/Lambert::ell);
    }
    force /= Control::N;
    v /= (Lambert::L);
    pair<double,double> par(v,force);
    return par;
}

/** Take a single step of ECMC. Active particle moves in +ve (right) direction.  
* We calculate the target energy, and add a Poisson distribued, thermal noise.
* Solve Lambert equation in root_right/root_left. Consider 1'st or 2'nd nn interactions.
*/
double
expon_step(vector<double> & x, vector<int> & count, long int & a, long int &b, double & displacement){//single ECMC step "a" is active particle
  int n= x.size();
  double r_step, l_step, ll_step, rr_step;
  long int right = (a+1) % n;
  long int rright = (a+2) % n;

  long int left= (a-1+n) % n;
  long int lleft= (a-2+n) % n;

  rr_step=Lambert::BIG; //if we only look at first-nearest neigbour this is not overwritten, and so not chosen.
  ll_step=Lambert::BIG;
  double sep_vec[4]={0,0,0,0};
  {//collision to the right: nearest neighbour
    double sep = x[right] - x[a];
    sep = sep < 0 ? sep + Lambert::L : sep;
    sep_vec[0]=sep;

    double energy = sep - Lambert::xmin < 0 ?  res_energy(sep,0) : Lambert::emin;
    energy +=  poisson(gen); 
  
    double root = root_right(energy,1); //smaller root of Lambert
    r_step = sep - root ;
      
    assert(r_step>=0);
  }

  if(Lambert::second==true)   {//collision to the right
    double sep = x[rright] - x[a];
    sep = sep < 0 ? sep + Lambert::L : sep;
    sep_vec[2]=sep;

    double energy = sep - Lambert::xmin2 < 0 ?  res_energy2(sep,0) : Lambert::emin2;
    energy +=  poisson(gen); 
  
    double root = root_right(energy,2); //smaller root of Lambert
    if(root < 0) {
      root=0; //no overtaking
    }
    rr_step = sep - root ;
    assert(rr_step>=0);
  }

  {//collision to the left, nearest neigbour
    double sep=x[a] - x[left];
    if(sep < 0 ) sep += Lambert::L;
    sep_vec[1]=-sep;
    

    double energy= sep - Lambert::xmin > 0 ? res_energy(sep,0) : Lambert::emin;
    energy +=  poisson(gen); 

    double root = root_left(energy,1); //larger root
    l_step =   root -sep;

    assert(l_step>0);
  }
  if(Lambert::second==true)   {//collision to the left
    double sep=x[a] - x[lleft];
    if( sep < 0 ) sep += Lambert::L;
    sep_vec[3]=-sep;

    double energy= sep - Lambert::xmin2 > 0 ? res_energy2(sep,0) : Lambert::emin2;
    energy +=  poisson(gen); 

    double root = root_left(energy,2); //larger root
    ll_step =   root -sep;

    assert(ll_step>0);
  }

  vector<double> candidate={ r_step, l_step, rr_step, ll_step};//four possible steps
  auto which = min_element( candidate.begin(), candidate.end() ) - candidate.begin();// find index of smallest candidate displacement
  x[a] += candidate[which];// move the active particle fowards
  displacement += sep_vec[which];
  x[a] =  x[a] > Lambert::L ? x[a]- Lambert::L : x[a]; //periodic 
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
  else if(which ==2){
    a=rright;
    b+=2;
    acounts[2]++;
  }
  else if(which ==3){
    a=lleft;
    b-=2;
    acounts[3]++;
  }
  count[a]++;
  return candidate[which];
}


/** fixed number of iterations of the update */
void
sweep(vector<double> & x, vector<int> &count, long int & a, long int &b, double & disp){ //N ECMC steps
  for (int i=0; i < 50 * Control::sample_time; i++){
    expon_step(x, count, a, b, disp);
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
    double s = expon_step(x,count,a,b, displacement);
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

/** In order to calculate \f( \rho_1 \f) and  \f( \rho_2 \f) we need multiple runs.
* Prepare the system, and perform Parm::iter runs.
*/
void
exp_mc(vector<double> &  x, vector<int>  & master, list<int> & lst, list<int>  & blst, list<int>  & zlst, list<double>  & xlst, ofstream & log_out){
  ofstream vout ("virial.dat");
  ofstream pout ("drift.dat");

  for(int i=0;i<Control::N;i++)
    x[i] = i * Lambert::xav * Lambert::scale_factor; //squeeze the chain by a factor Lambert::factor from zero pressure

  for (int l=0; l<Control::iter; l++){ //main loop for running multiple simulations

    if( l % (1024*8) == 0) {
      cout<<"l "<<l/1024<<endl;
      log_out<<"l "<<l/1024<<" of "<<Control::iter/1024<<endl<<flush;
    }
    if( l % (1024*128) == 0)   save(x, master, lst, blst, zlst, xlst); //so we have something to analyse if we don't finish a long run

    if (!Lambert::cold) { //equilibrate between each run
      long int a=0;
      long int b=0;
      double d=0;
      vector<int> count(Control::N, 0);
      sweep(x,count,a,b,d);
    }
    else{//reset to ground state
      for(int i=0;i<Control::N;i++)
	x[i] = i * Lambert::xav * Lambert::scale_factor; //squeeze the chain by a factor Lambert::factor from zero pressure
    }

    long int a = 0;
    long int a0 = a; //remember origin
    long int b = 0; // "b" is almost "a" but unwrapped to follow activity
    vector<int> count(Control::N,0);//count site visits

    double displacement=0;
    chain( x , count, a, b, displacement); //chain() or replace by sweep()
    pout<< Lambert::T * Control::N * displacement/Control::sample_time/Lambert::L << endl;
    vout << virial(x,1).first<<" " << virial(x,2).first << endl;

    lst.push_back( count[a] ); // H curve, rho_2
    blst.push_back( b ); //displacement, rho_1
    xlst.push_back( x[a] ); //final position
    zlst.push_back( count[a0] ); //origin visits

    for(int i=0; i<Control::N; i++){//global count curve, not analysed
      int pos = (a0+i) % Control::N;
      master[i] += count[pos];
    }
  }//Lambert::inter
}

void
init(){ /** Rescale parameters for solution of Lambert equation */
//  assert(Lambert::p2>0);
 // assert(Lambert::p>0);
  if ( fabs(Lambert::p) < 1.e-12) Lambert::p = 1.e-14;
  if ( fabs(Lambert::p2) < 1.e-12) Lambert::p2 = 1.e-14;
  cout<<red;
  Lambert::pbar = Lambert::p * Lambert::ell/Lambert::a[0];
  Lambert::pbar2 = Lambert::p2 * Lambert::ell/Lambert::a[1];
  if(Lambert::second) assert(Lambert::p2 != 0);
  Lambert::xmin= - Lambert::ell * log(Lambert::pbar); //position of minimum of function
  Lambert::xmin2= - Lambert::ell * log(Lambert::pbar2); //position of minimum of function
  Lambert::emin= res_energy(Lambert::xmin, 0); //minimum energy
  Lambert::emin2= res_energy2(Lambert::xmin2, 0); //minimum energy
  cout<<"xmin emin pbar\t"<<Lambert::xmin<<" "<<Lambert::emin<<" "<<Lambert::pbar<<endl;
  cout<<"xmin emin pbar\t"<<Lambert::xmin2<<" "<<Lambert::emin2<<" "<<Lambert::pbar2<<endl;
  cout<<reset;
}

/** Find average separation by numerical integration: only valid for nearest neighbour interactions */
double
integrate(){
  using namespace boost::math::quadrature;
  auto f1 = [](double t) { return    exp( (-Lambert::a[0]*exp(-t/Lambert::ell) - Lambert::p*t)/Lambert::T ); };
  auto f2 = [](double t) { return t*exp( (-Lambert::a[0]*exp(-t/Lambert::ell) - Lambert::p*t)/Lambert::T ); };
  double error1, error2;
  double Q1 = gauss_kronrod<double, 61>::integrate(f1, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error1);
  double Q2 = gauss_kronrod<double, 61>::integrate(f2, 0, numeric_limits<double>::infinity(), 5, 1e-12, &error2);
  cout <<yellow<< "average length = " <<Q2/Q1<< " " <<reset<< endl;
  return Q2/Q1;
}

/** This is the main, historic run that calculates \f(\rho_1\f) and \f(\rho_2\f) 
* Here create the empty lists which are filled in the simulation with data
* and then call mc()
*/
void
run_distributions(ofstream & log_out, vector<double> & x){

  vector<int> count(Control::N);
  vector<int> master(Control::N);
  list<int> lst; //measures h, rho_2
  list<int> blst; //unwrapped final index, rho_1 
  list<int> zlst; // returns to origin, local time
  list<double> xlst; //final position in x of chain

  exp_mc(x, master, lst, blst, zlst, xlst, log_out);
  save(x, master, lst, blst, zlst, xlst);
 
}

/** Evolution of \f[ S(q) \f] and perhaps entropy during simulation, don't calculate the main distribution functions */
void
long_run(vector<double> & x){
  ofstream sq_out("sq.dat");
  ofstream vout ("virial.dat");
  ofstream fout ("force.dat");
  long int a=0;
  long int b=0;
  double displacement=0;
  vector<int> count(Control::N);

  for(int i=0;i<Control::N;i++)
    x[i] = i * Lambert::xav * Lambert::scale_factor; //squeeze the chain by a factor Lambert::factor from zero pressure

  for(int i=0;i<Control::sweeps;i++){
    chain(x, count,  a, b, displacement);
    if(i> Control::sweeps/8){
        sq_out<<Sq(x)<<endl;
        double EstimatedPressure=  virial(x,1).first + Lambert::T* Control::N/Lambert::L + virial(x,2).first ;
        vout<< EstimatedPressure<<endl;
        double Estimatedforce=  virial(x,1).second+ 2*virial(x,2).second;
	fout<<Estimatedforce<<endl;
	if(i%1024==0) cout<<"i "<<i/1024<<"/"<<Control::sweeps/1024<<endl;
    }
  }
  ofstream xout ("pos.dat"); 
  for (double i : x) xout <<i<<"\n";
}
void test(){// this checks that the algebra for roots is correct, nenergy should be zero
  for(int i=0;i<5;i++){
    double energy =  Lambert::emin + 5*poisson(gen);
    double root = root_right(energy,1); //smaller root of Lambert
    double nenergy = Lambert::a[0] * exp(-root/Lambert::ell) + Lambert::p* root - energy;
    //	    cout <<"test energy, root nenergy "<<energy<<" "<<root<<" "<<nenergy<<endl;
    assert(fabs(nenergy) < 1.e-14);

    root = root_left(energy,1); //smaller root of Lambert
    nenergy = Lambert::a[0] * exp(-root/Lambert::ell) + Lambert::p* root - energy;
    //	    cout <<"test energy, root nenergy "<<energy<<" "<<root<<" "<<nenergy<<endl;
    assert(fabs(nenergy) < 1.e-14);
  }
  cout<<red<<"Energy checks"<<reset<<endl;
}

/** ECMC simulation of exponential interactions, first and second neigbourg */
int
main(){
  code();
  vector<double> x(Control::N);
  ofstream log_out ("lambert.log");
  log_out<<"git commit string "<<GIT_COMMIT<<endl;
  
  init();
  Lambert::xav=integrate();
  Lambert::xav= 4.69703;
  cout<<"Reset xav"<<endl;
  //test();

  Lambert::L = Lambert::xav * Control::N * Lambert::scale_factor;
  cout<<"Length " <<Lambert::L<<endl;
  cout<<" a= "<< Lambert::a[0] <<" "<<Lambert::a[1]<< " p= " << Lambert::p <<" "<<Lambert::p2<<" T=  "<<Lambert::T<<endl;
  cout<<" rho "<< Control::N/Lambert::L<<endl;
  Timer t1;
  
  if(Lambert::single_run)
    long_run(x);//structure factor
  else 
    run_distributions(log_out, x);//rho_1 and rho_2

  double nn= virial(x,1).first;
  double ns= virial(x,1).second;
  double nnn= virial(x,2).first;
  double nns= virial(x,2).second;
  cout<<yellow<<"v1 v2 "<<nn<<" "<<nnn<<reset<<endl;
  cout<<yellow<<"v1s v2s "<<ns<<" "<<nns<<reset<<endl;
  double EstimatedPressure=  nn + Lambert::T* Control::N/Lambert::L + nnn ;
  log_out<<"VirialPressure "<< EstimatedPressure <<endl ;
  cout<<"VirialPressure "<< EstimatedPressure <<endl ;
  log_out<<"Sum_factor_fields " << Lambert::p+2*Lambert::p2<<endl;
  cout<<"Sum_factor_fields " << Lambert::p+2*Lambert::p2<<endl;
  double s=acounts[0] + acounts[1]+acounts[2]+acounts[3];
  cout<<"counts "<<acounts[0]/s<<" "<<acounts[1]/s<<" "<<acounts[2]/s<<" "<<acounts[3]/s<<endl;
  t1.stats(log_out);
  return 0;
}
