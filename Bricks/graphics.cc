// https:www.cairographics.org/manual/cairo-Paths.html fix
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <unistd.h>
#include <cassert>
#include "graphics.h"

void Graphics::draw(std::vector<double> &p, double FPS, std::vector<int> &active) {

  //  const int FPS=50;//frames per second
  struct timespec tm = {0, (long int)(1000000000 / FPS)}; // Sleep in nanoseconds between frames
  XEvent event;                                           // check if window closed and finish
  if (XCheckWindowEvent(dsp, da, DestroyNotify, &event)) {
    XCloseDisplay(dsp);
    exit(1);
  }

  double gamma = Dim / (lmax - lmin + diam); // scaling between physical units and pixels
  double alpha = gamma * diam / 2.;
  char a[10];

  cairo_push_group(cr);                   // start drawing
  cairo_set_source_rgb(cr, 0.0, 0.1, .1); // dark green background
  cairo_paint(cr);                        // clear screen with green

  for (int i = 0; i < Np; i++) { // place the particles in the graphics buffer, without drawing
    cairo_new_sub_path(cr);
    if (active[i] == 0) {
      cairo_set_source_rgb(cr, 0.79, 0.29, 0.29); // dark red for particles
    } else {
      cairo_set_source_rgb(cr, 1, 1, 1);
    }
    cairo_arc(cr, (alpha + gamma * (i + lmin)), 250 - 4 * p[i], 3.5 * alpha, 0, 2 * M_PI);
    cairo_fill(cr); // draw all particles with solid color
  }

  cairo_move_to(cr, Dim - 120, Dim - 35); // print frame number to graphics window
  // snprintf(a,4, "%d",l);
  cairo_set_font_size(cr, 20);
  cairo_show_text(cr, a);

  cairo_pop_group_to_source(cr); // finished drawing operations for this set of positions
  cairo_paint(cr);               // send to screen
  cairo_surface_flush(sfc);      // send to x11
  nanosleep(&tm, NULL);          // this sets the animations speed
  XFlush(dsp);                   // sync X11 to cairo
}

Graphics::Graphics(int N, int Pix, double dmin, double dmax, double diameter) {
  Np = N;
  count = 1000;
  assert(dmin < dmax);
  assert(Np > 0);
  lmin = dmin;
  lmax = dmax;
  Dim = Pix;
  diam = diameter;

  if ((dsp = XOpenDisplay(NULL)) == NULL)
    exit(1); // window management X11, and cairo graphics
  int screen = DefaultScreen(dsp);
  da = XCreateSimpleWindow(dsp, DefaultRootWindow(dsp), 0, 0, Dim, 500, 0, 0, 0);
  XMapWindow(dsp, da);
  sfc = cairo_xlib_surface_create(dsp, da, DefaultVisual(dsp, screen), Dim, Dim);
  cairo_xlib_surface_set_size(sfc, Dim, Dim - 400);
  cr = cairo_create(sfc);
} // empty window now on screen

Graphics::~Graphics() {
  cairo_destroy(cr);
  cairo_surface_destroy(sfc);
} // clean up function
