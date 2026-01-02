#include <vector>
#include <random>
#include <set>
#include <iostream>
#include <memory>
#include <cassert>
#include <fstream>
#include "version.h"
#include "graphics.h"
#include "timer.h"
#include <unistd.h>
#include <queue>

using namespace std;

namespace Param {
int N = 1024;             // System size
const double big = 1.e14; // large time for empty sites
const double beta = 0.4;  // for rate constants

int sweep_length = N;

const double reset = 1000000.; // subtract when time is too big, in sweeps

const bool do_multi = false;
int nrun = 80000; // for multisweep

const bool dopeak = false;
int peak = 300;      // multiple on same site
int nw = (int)N / 8; // number of builders with dopeak==false

int ns = 1024; // number of steps in single simul
double tstep = 10;

int nreset = 0;

const bool GRAPH = false;
long int n_event = 0;
} // namespace Param

enum direction { no_event, left_event, right_event };
// the main difficulty in the code is that we need an object that can be in a 1-dimensional physical object, with indexes, while also being
// sorted into a tree using the event time as the sorting variable.
class siteEvent {
public:
  int site;                                       // position on chain
  double time;                                    // time to next event
  multiset<shared_ptr<siteEvent>>::iterator link; // a back link from the lattice to the sorted times
  // set<shared_ptr<siteEvent>>::iterator link; // a back link from the lattice to the sorted times
  siteEvent(int _site, double _time) : site(_site), time(_time) {}; // set site data
  struct compareTime {                                              // sort on pointers to a siteEvent
    bool operator()(const shared_ptr<siteEvent> &d1, const shared_ptr<siteEvent> &d2) const { return (d1->time) < (d2->time); }
  };
  void print() const { cout << site << "\tlink " << &*link << "\tt " << time << endl; }
};

// "shared_ptr" automatically implements memory management
multiset<shared_ptr<siteEvent>, siteEvent::compareTime> events; // time-sorted collection of siteEvent
// set <shared_ptr<siteEvent>,  siteEvent::compareTime> events; //time-sorted collection of siteEvent
// priority_queue<shared_ptr<siteEvent> , vector<shared_ptr<siteEvent>> > pqv;
vector<shared_ptr<siteEvent>> wall; // vector of the same set of siteEvent
vector<int> z(Param::N, 0);         // interface shape
vector<direction> dir(Param::N);    // direction of next event at site
vector<int> rho(Param::N);          // number of builders at site
vector<int> count_plus(Param::N);
vector<int> count_minus(Param::N);

void reset_mem() {
  z.resize(Param::N);
  dir.resize(Param::N);
  rho.resize(Param::N);
  count_plus.resize(Param::N);
  count_minus.resize(Param::N);
  Param::nw = (int)Param::N / 8;
  Param::sweep_length = Param::N;
}

ofstream wout("width.dat");
ofstream four_out("fourier.dat");
ofstream zout("zcor.dat");
ofstream coout("count.dat");
ofstream z2out("z2.dat");
ofstream sout("second.dat");
ofstream cycl("cyclic.dat");
Graphics *g;

mt19937 mt(random_device{}()); // seed mersenne-twister from quality source
// mt19937 mt(0); //seed mersenne-twister to repeat the same run
uniform_real_distribution<> uni(0., 1.);
exponential_distribution<> poisson(1.);

static string reset = "\u001b[0m";
static string yellow = "\u001b[33m";
static string red = "\u001b[31m";

double rate(double zz) { // rate function:  r(z) * r(1-z) == 1 for the continuum equations to be valid
  const double pre = sqrt(1. / exp(Param::beta));
  return pre * exp(Param::beta * zz);
}

void print_code() { // print the code to this program
#include "bricks.hex"
  ofstream opos("source.code");
  opos << "git commit string " << GIT_COMMIT << endl;
  for (auto i = 0U; i < bricks_cc_len; i++)
    opos << bricks_cc[i];
}

void print_all() { // debugging
  cout << endl;
  cout << "sorted " << endl;
  for (const auto &e : events) { // print out sorted contents
    cout << e->site << " " << e->time << endl;
  }
  cout << endl;
  cout << "full " << endl;
  for (int i = 0; i < Param::N; i++) { // full information for each site
    wall[i]->print();
  }
  cout << endl;

  vector<int> h(Param::N);
  h[0] = 0;
  cout << "h ";
  for (int i = 1; i < Param::N; i++) {
    h[i] = h[i - 1] - z[i];
    cout << h[i] << " ";
  }
  cout << endl;

  cout << "z ";
  for (int i = 0; i < Param::N; i++) {
    cout << z[i] << " ";
  }
  cout << endl << endl;
  cout << "rho ";
  for (int i = 0; i < Param::N; i++) {
    cout << rho[i] << " ";
  }
  cout << endl;
} // end debug

void fourier_modes() {
  vector<double> h(Param::N, 0.);
  double q = 2. * M_PI / Param::N;

  double c = 0;
  double s = 0;
  double srho = 0;
  double crho = 0;
  for (int i = 1; i < Param::N; i++) { // pair correlations between nn steps
    h[i] = h[i - 1] - z[i];
  }
  for (int i = 0; i < Param::N; i++) { // pair correlations between nn steps
    c += h[i] * cos(q * i);
    s += h[i] * sin(q * i);
    crho += rho[i] * cos(q * i);
    srho += rho[i] * sin(q * i);
  }
  c /= sqrt(Param::N);
  s /= sqrt(Param::N);
  crho /= sqrt(Param::N);
  srho /= sqrt(Param::N);
  four_out << c << "\t " << s << "\t " << crho << "\t " << srho << endl;
}

void stats() {
  vector<double> h(Param::N, 0.);
  double mean = 0;
  double correlation = z[0] * z[Param::N - 1];
  for (int i = 1; i < Param::N; i++) { // pair correlations between nn steps
    h[i] = h[i - 1] - z[i];
    mean += h[i];
    correlation += z[i] * z[i - 1];
  }
  zout << correlation << endl;
  mean /= Param::N;

  double var = 0;
  double z2 = 0;
  for (int i = 0; i < Param::N; i++) { // distribution of size of steps, often "2" for parity reasons dominates
    z2 += z[i] * z[i];
    if (z[i] >= 0)
      count_plus[z[i]]++;
    if (z[i] < 0)
      count_minus[-z[i]]++;
    var += (h[i] - mean) * (h[i] - mean);
  }
  var /= Param::N;
  var = sqrt(var);
  wout << var << endl;
  z2out << z2 / Param::N << endl;

  double second = 0; // second nearest neighbour correlation
  for (int i = 0; i < Param::N - 2; i++) {
    second += z[i] * z[i + 2];
  }
  sout << second << endl;

  if (false) {
    vector<double> full(Param::N, 0.); // correlations over full range of separations
    for (int i = 0; i < Param::N; i++) {
      for (int j = 0; j < Param::N; j++) {
        int l = (i + j) % Param::N;
        full[j] += z[i] * z[l] / (double(Param::N));
      }
    }
    for (int i = 0; i < Param::N; i++) {
      cycl << full[i] << " ";
    }
    cycl << endl;
  }
} // end stats

#include <boost/type_index.hpp>
void reset_site(const int j, const double time) { // calculate rates after change in occupation
  double next_time;
  direction next_direction;
  if (rho[j] != 0) {
    double time_right = time + poisson(mt) / (rate(z[j]) * rho[j]);
    double time_left = time + poisson(mt) / (rate(-z[j]) * rho[j]);
    next_time = min(time_right, time_left);
    next_direction = time_right < time_left ? right_event : left_event;
  } else {
    next_time = Param::big;
    next_direction = no_event;
  }
  if (0) {
    auto what = wall[j]->link;
    cout << "Type " << endl;
    cout << boost::typeindex::type_id_with_cvr<decltype(what)>().pretty_name() << endl;
    exit(0);
  }
  if (dir[j] != no_event)
    events.erase(wall[j]->link); // throw away last values
  shared_ptr<siteEvent> value(new siteEvent(j, next_time));
  dir[j] = next_direction;
  wall[j] = value;
  // if(true){
  if (next_direction != no_event) {
    auto link = events.insert(value);
    wall[j]->link = link;
  }
}

void build() { // add builders to new system and calculate times of events for each site
  fill(z.begin(), z.end(), 0);
  fill(rho.begin(), rho.end(), 0);

  if (Param::dopeak) { // delta function of many builders in middle
    rho[Param::N / 2] = Param::peak;
  } else {
    for (int i = 0; i < Param::nw; i++) { // random placement
      int site = int(Param::N * uni(mt));
      rho[site] += 2; // this to have constant parity
    }
  }

  if (Param::GRAPH) {
    for (int i = 0; i < Param::N; i++) {
      z[i] = 0 * 4. * sin(i * 2. * M_PI / Param::N); // odd initial condition for visualisation
    }
  }

  wall.clear();
  events.clear();
  for (int i = 0; i < Param::N; i++) { // dummy events  to create the data structures
    shared_ptr<siteEvent> value(new siteEvent(i, poisson(mt)));
    dir[i] = right_event;
    wall.push_back(value); // pointers are inserted into sorted collection and also vector
    auto link = events.insert(value);
    (*link)->link = link;
  }
  for (int i = 0; i < Param::N; i++) { // now calculate real events from initial geometry
    double time = 0;
    reset_site(i, time);
  }
  // print_all();
}

void manage_event(const int l, const double current_time) {
  Param::n_event++;
  if (dir[l] == no_event) {
    wall[l]->print();
    abort();
  }
  int lplus = (l + 1) % Param::N;
  int lminus = (l - 1 + Param::N) % Param::N;

  if (dir[l] == right_event) { // page 233 in springer, page 8 arxiv
    rho[l] -= 1;
    rho[lplus] += 1;
    z[l] -= 1;
    z[lplus] += 1;
    reset_site(lplus, current_time);
  } else {
    rho[l] -= 1;
    z[l] += 1;
    rho[lminus] += 1;
    z[lminus] -= 1;
    reset_site(lminus, current_time);
  }
  reset_site(l, current_time);
}

void reset_time() { // stop times getting too big and overflowing
  for (int i = 0; i < Param::N; i++) {
    if (wall[i]->time != Param::big)
      wall[i]->time -= Param::reset;
  }
  Param::nreset++;
}

void time_step(int l) { // move forwards to the time l*tstep
  // cout<<"l "<<l*Param::tstep<<endl;

  double next_end_time = l * Param::tstep;
  while (true) {
    auto first_event = events.begin();
    //    cout << Param::N <<" "<< (*first_event)->time <<endl;
    // sleep(1);
    double event_time = (*first_event)->time;

    if (event_time > next_end_time) {
      // double delta= event_time - next_end_time;
      // cout<<"delta "<<delta<<endl;
      //       for (int i=0;i<Param::N;i++){
      // if (wall[i]->time != Param::big) wall[i]->time -= delta ;
      // }
      //  cout<<"return \n\n";
      //     cout<<"event_time "<<event_time<<" next_time "<<next_end_time<<endl;
      return;
    } else {
      //      cout<<"Manage "<<endl;
      manage_event((*first_event)->site, event_time);
      //      cout<<"manage event_time "<<event_time<<" next_time "<<next_end_time<<endl;
    }
  }
}

void sweep() { // fixed number of events
  for (int i = 0; i < Param::sweep_length; i++) {
    auto first_event = events.begin();
    double current_time = (*first_event)->time;
    if (current_time > Param::reset)
      reset_time();
    manage_event((*first_event)->site, current_time);
  }
}

void theory() { // equilibrium distribution functions from the paper
  double lambda = 0.25;
  double peven = 0;
  double podd = 0;
  vector<long int> fact(20);
  fact[0] = 1;
  for (int i = 1; i < 20; i++) {
    fact[i] = fact[i - 1] * i;
  }
  for (int i = 0; i < 5; i++) {
    peven += pow(lambda, 2 * i) / fact[2 * i];
    podd += pow(lambda, 2 * i + 1) / fact[2 * i + 1];
  }
  peven *= exp(-lambda);
  podd *= exp(-lambda);
  int L = 10;
  //  for(int i=0;i<5;i++){
  //  cout<<" i ri "<<i <<" "<<rate(i)<<endl;
  //}
  vector<double> R(L);
  R[0] = 1;
  for (int i = 1; i < L; i++) {
    R[i] = R[i - 1] / rate(i);
  }
  for (int i = 0; i < L; i++) {
    if (i % 2 == 0)
      R[i] = R[i] * peven;
    if (i % 2 == 1)
      R[i] = R[i] * podd;
  }
  double normR = R[0];
  for (int i = 1; i < L; i++) {
    normR += 2 * R[i];
  }
  for (int i = 0; i < L; i++) {
    R[i] /= normR;
    cout << "Prob " << i << " " << R[i] << endl;
  }
} // end theory

void single_run() {
  build();
  cout << red << "single run " << reset << endl;
  static int first = 0;
  for (int i = 0; i < Param::ns; i++) {
    if (Param::GRAPH) {
      vector<double> h(Param::N, 0.);
      for (int j = 1; j < Param::N; j++) {
        h[j] = h[j - 1] - z[j] / 10.;
      }

      int FPS = 5;
      g->draw(h, FPS, rho);
      if (first == 0) {
        first = 1;
        sleep(2);
      }
    } // end GRAPH
    time_step(i + 1);
    //    sweep();
    // stats();
    fourier_modes();
    if (i % 100000 == 0)
      cout << "loop " << i << "\n" << flush;
  }
  for (int i = 0; i < Param::N; i++) {
    if (count_plus[i] != 0)
      coout << i << " " << count_plus[i] << endl;
    if (count_minus[i] != 0)
      coout << -i << " " << count_minus[i] << endl;
  }
}

void multi_run() {
  cout << red << "Multi run " << reset << endl;
  vector<int> cumul(Param::N, 0);
  vector<int> builders(Param::N, 0);
  for (int i = 0; i < Param::nrun; i++) {
    if (i % 1000 == 0) {
      cout << "loop " << i << "\n" << flush;
    }
    build();
    //    sweep();
    time_step(1);
    int h = 0;
    for (int j = 0; j < Param::N; j++) {
      h += z[j];
      cumul[j] -= h;         // profile
      builders[j] += rho[j]; // builders
    }
  } // repetitions of simulation
  ofstream mout("multi.dat"); // profile
  ofstream fout("final.dat"); // builders
  double mean_cumul = 0;
  double n_cumul = 0;
  double mean_build = 0;
  double n_build = 0;
  for (int i = 0; i < Param::N; i++) {
    mout << cumul[i] << " "; // profile
    mean_cumul += cumul[i] * i;
    n_cumul += cumul[i];
    if (i != Param::N / 2) {
      mean_build += builders[i] * i;
      n_build += builders[i];
    }
    fout << builders[i] << " "; // builders
  }
  fout << endl;
  mean_cumul /= n_cumul;
  mean_build /= n_build;
  cout << yellow << mean_cumul << " " << mean_build << reset << endl;
  double v_cumul = 0;
  double v_build = 0;
  for (int i = 0; i < Param::N; i++) {
    v_cumul += cumul[i] * (i - mean_cumul) * (i - mean_cumul) / n_cumul;
    if (i != Param::N / 2) {
      v_build += builders[i] * (i - mean_build) * (i - mean_build) / n_build;
    }
  }
  v_cumul = sqrt(v_cumul);
  v_build = sqrt(v_build);
  mout << endl;
  cout << v_cumul << " " << v_build << " " << cumul[Param::N / 2] / double(Param::nrun) << endl;
}

int main(int argc, char *argv[]) {
  Timer tm;
  print_code();
  if (argc > 1) {
    Param::peak = stoi(argv[1]);
    cout << "Peak " << Param::peak << endl;
  }
  if (argc > 2) {
    Param::nrun = 1000 * stoi(argv[2]);
    cout << "nrun " << Param::nrun << endl;
  }
  if (argc > 3) {
    Param::tstep = stod(argv[3]);
    cout << "tstep " << Param::tstep << endl;
  }
  if (argc > 4) {
    Param::N = stoi(argv[4]);
    cout << "N " << Param::N << endl;
  }

  if (argc > 5) {
    Param::ns = stoi(argv[5]);
    cout << "ns " << Param::ns << endl;
    reset_mem();
  }
  // theory();
  if (Param::GRAPH) {
    double dmin = 0.02 * Param::N;
    double dmax = Param::N * 1.05;
    double diameter = 0.2;
    g = new Graphics(Param::N, 1800, dmin, dmax, diameter);
  }
  if (Param::do_multi)
    multi_run();
  else
    single_run();
  auto first_event = events.begin();
  double current_time = (*first_event)->time;
  cout << yellow << "total_time " << Param::nreset * Param::reset + current_time << endl;
  cout << Param::N << "\t Events " << Param::n_event << reset << endl;
  tm.stats();
  return 0;
}
