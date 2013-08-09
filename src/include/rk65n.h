#ifndef _RK65N_H_INCLUDED
#define _RK65N_H_INCLUDED

#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sundials/sundials_types.h> /* definition of type realtype */
#include <math.h>
#include "neuronmodel.h"

/* #define min(X,Y) (X<Y ? X:Y) */
/* #define max(X,Y) (X>Y ? X:Y) */

#define SIXTH 0.166666666666667

#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */
class RK65n {
 private:
  double a[9][8];
  double b[9];
  double newdt, dtx, theEps;
  N_Vector Y[9];
  N_Vector F[9];
  double *y5;
  double aF;
  double delta;
  int N;
  double maxdt, eps, abseps, releps; 
  int i, j, k;
  realtype *ff, *yy;
  N_Vector yt;
public:

  /* arguments dim, maxdt, error ineps, abstol, reltol */
  RK65n(int, double, double, double, double);
  ~RK65n();
  /* argument y_vector, ydot_vector, and dt */
  double integrate(double, N_Vector, N_Vector, NeuronModel *, double);
  int dim() { return N; }
};

#endif
