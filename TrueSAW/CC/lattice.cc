#include <vector>
#include <iostream>
#include <random>
#include <cassert>
#include <fstream>
#include <list>
#include <algorithm> // std::shuffle
#include "timer.h"
#include "graphics.h"
#include "version.h"

namespace Param {
const int p = 2;              // lattice is p*Np sites
int Np = 256;                 // number of particles
const double p0 = 1 - 1. / p; // probability
int miter = 1024 * 128 * 4;   // number of samples
const int do_equilibrate = 0; // hot or cold start
const int GRAPH = 0;          // animation on screen
} // namespace Param

using namespace std;

static string reset = "\u001b[0m";
static string yellow = "\u001b[33m";
static string red = "\u001b[31m";

Graphics *g;

mt19937_64 gen(random_device{}()); // seed random numbers from high quality source
uniform_real_distribution<double> uniform(0., 1.);

void code() { // print the code to this program
#include "lattice.hex"
  for (auto i = 0U; i < lattice_cc_len; i++)
    cout << lattice_cc[i];
  exit(0);
}

inline void step(vector<int> &x, vector<int> &count, int &a, int &b) { // single ECMC step
  // "x[i]" is the position of the particle "i".
  // "n=x.size()" is the number of particles
  // "p*n" is the number of sites in the periodic system
  // "a" is the pointer position, or the number of the active particle x[a] is it's position at the start of the routine.
  // "b" is just "a" but without the periodic bc, does not need to be followed to understand logic
  // in all my runs p0=1/2
  // eq number correspond to arxiv: https://arxiv.org/pdf/2306.13059
  int n = x.size();
  int L = Param::p * n;
  int right = (a + 1) % n;
  int left = (a - 1 + n) % n; // a, right, left are indexes

  if ((x[a] + 1) % L != x[right]) { // there is no particle immediately to the right of the particle "a"
    x[a] = (x[a] + 1) % L;          // jump right
    if (uniform(gen) < Param::p0) { // empty, 9b, alpha'
    } else {                        // transfer motion to left 9b alpha
      b--;
      a = left;
    }
  } else {                          // site to the right is full
    if (uniform(gen) < Param::p0) { // transfer active particle, no motion 10b alpha'
      a = right;
      b++;
    } else { // keep same particle
      // 10b alpha
    }
  }
  count[a]++;
}

void sweep(vector<int> &x, vector<int> &count, int &a, int &b) { // N ECMC steps
  int factor = 1;
  if (Param::GRAPH)
    factor = 100;
  for (auto i = 0U; i < factor * x.size(); i++) { // note 1/2 here, prevents wrap-around in smaller systems
    step(x, count, a, b);
    if (Param::GRAPH)
      g->draw(x, 60, i, a);
  }
}

void save(vector<int> const &x, list<int> const &lst, list<int> const &blst, list<int> &sep, vector<int> &count) {
  ofstream opos("pos.dat");
  ofstream lout("list.dat");
  ofstream bout("blist.dat");
  ofstream sout("sep.dat");
  ofstream ctout("count.dat");

  for (double i : x)
    opos << i << "\n";
  for (int n : lst)
    lout << n << "\n";
  for (int n : blst)
    bout << n << "\n";
  for (int n : sep)
    sout << n << "\n";
  for (int n : count)
    ctout << n << "\n";
}

void sf(vector<int> &x) {
  vector<int> trial(Param::Np * Param::p);
  for (auto j = 0; j < Param::Np * Param::p; j++)
    trial[j] = j % Param::p == 0 ? 1 : 0;
  shuffle(trial.begin(), trial.end(), gen); // random permutation of trial
  int found = 0;
  for (auto i = 0; i < Param::p * Param::Np; i++) {
    if (trial[i] == 1) {
      x[found] = i;
      found++;
    }
  }
  assert(found == Param::Np);
}

int main(int argc, char *argv[]) {
  if (argc == 2 && argv[1] == string("code"))
    code();
  if (argc > 1)
    Param::Np = stoi(argv[1]);
  if (argc > 2)
    Param::miter = 1024 * stoi(argv[2]);
  ofstream log_out("tasep.log");
  log_out << GIT_COMMIT << endl;
  Timer t1;
  cout << yellow << "Param::p=" << Param::p << "\tNp=" << Param::Np << "\tmiter=" << Param::miter
       << "\tdo_equilibrate=" << Param::do_equilibrate << reset << endl;
  log_out << "Param::p=" << Param::p << "\tNp=" << Param::Np << "\tmiter=" << Param::miter << "\tdo_equilibrate=" << Param::do_equilibrate
          << endl;
  if (Param::GRAPH) {
    double dmin = 0.02 * Param::p * Param::Np;
    double dmax = Param::p * Param::Np * 1.05;
    double diameter = 0.2;
    g = new Graphics(Param::Np, 1400, dmin, dmax, diameter);
  }
  vector<int> x(Param::Np);
  vector<int> count(Param::Np);
  vector<int> total(Param::Np);
  list<int> blst;
  list<int> lst;
  list<int> sep;

  for (int i = 0; i < Param::miter; i++) { // main loop over samples
    if (i % (8 * 1024) == 0) {
      cout << "l " << i / 1024 << endl;
      log_out << "l " << i / 1024 << " of " << Param::miter / 1024 << endl << flush;
    }
    int b = 0; // b is a, but unwrapped to follow displacement of activity
    int a = 0;
    for (auto j = 0; j < Param::Np; j++)
      x[j] = Param::p * j; // restart
    if (Param::do_equilibrate)
      sf(x);
    std::fill(count.begin(), count.end(), 0); // reset count of site visits
    sweep(x, count, a, b);
    blst.push_back(b);       // displacement
    lst.push_back(count[a]); // for H
    for (auto m = 0U; m < count.size(); m++)
      total[m] += count[m]; // total site visits
  }

  save(x, lst, blst, sep, total);
  cout << red << "Done" << reset << endl;
  log_out << "Done" << endl;
  t1.stats(log_out);
  return 0;
}
