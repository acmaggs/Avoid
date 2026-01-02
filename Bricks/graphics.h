// We use the casio graphics library,2018
// documentation: https://www.cairographics.org/manual/cairo-Paths.html
#ifndef _Graphics
#define _Graphics
#include <unistd.h>
#include <cassert>
#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
//#include <X11/Xlib.h>
//#include <X11/Xutil.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

class Graphics{
 private:
  int Np;
  double lmin, lmax,Dim,diam;
  cairo_surface_t *sfc;
  cairo_t *cr;
  Display *dsp;
  Drawable da;
  int count;
 public:
  Graphics(int N, int Pix, double dmin, double dmax,double diameter);
  //some error messages on trying to duplicate the window
  Graphics & operator=(Graphics &g) = delete; // {// stop the program when copying the window
  void draw(std::vector<double> &, double, std::vector<int> & );//draw the particles
  ~Graphics();
};

#endif
