void
step(vector<double> & x, int & a){
  int n= x.size(); // number of particles
  //"x" array of particle positions
 // "a" is the active particle at the start of the routine
  double r_step, l_step, energy;
  int right = (a+1) % n;
  int left= (a-1+n) % n;

  //collision to the right
  double sep=x[right]-x[a];
  energy =  -Param::T*log( uniform(gen)); // the quasi-metroplis criterion for a single bond
  if ( sep < 0 ) 
    energy += sep*sep/2.; // elastic spring energy
  r_step= sqrt(2*energy)+sep;

  //collision to the left
  sep=x[a] - x[left];
  energy = -Param::T*log( uniform(gen)); // the quasi-metroplis criterion for a single bond
  if( sep > 0 )
    energy += sep*sep/2.; // elastic spring energy
  l_step =  sqrt( energy*2) -sep;

  if (r_step<l_step){ //now choose between two candidates
    x[a] += r_step;
    a=right;
  }
  else{
    x[a] += l_step;
    a=left;
  }
}
