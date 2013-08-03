#ifndef _MM_H_INCLUDED_
#define _MM_H_INCLUDED_

/* /\* cvode libraries *\/ */
/* #include <cvode/cvode.h>             /\* prototypes for CVODE fcts., consts. *\/ */
/* #include <nvector/nvector_serial.h>  /\* serial N_Vector types, fcts., macros *\/ */
/* #include <cvode/cvode_dense.h>       /\* prototype for CVDense *\/ */
/* #include <sundials/sundials_dense.h> /\* definitions DlsMat DENSE_ELEM *\/ */
/* #include <sundials/sundials_types.h> /\* definition of type realtype *\/ */

#include "neuron.h"
#include "synapse.h"
#include "electrode.h"
//#include "simplecell.h"
#include "pdneuron.h"
#include "abneuron.h"
#include "gap.h"
#include "neuronmodel.h"
#include "rk.h"
#include "CVodeSolver.h"

int derivative(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int f(realtype t, N_Vector y, N_Vector dy, void* user_data);
void swap(N_Vector y, N_Vector yout);
static int check_flag(void *flagvalue, char *funcname, int opt);

#define NEQ 17 /* Number of equations */
//#define NEQ 8
//#define NEQ 4
//#define TSTEP  0.095367432    /* stepsize */
#define TSTEP 0.05
#define FREQUENCY 0.75
#define NUMCYCLES 30
//#define TSTOP (NUMCYCLES/FREQUENCY)*1000 // duration for constant
                                         // frequency sine waves with
                                         // integer multiple number of
                                         // cycles
#define T0 0.095367432            /* start time */
#define TSTOP 3000  /* duration */

#define NUM_PARS 15 /* number of parameters */

/* set initial calcium concentration at 10e-4 mM */
double PD_INIVARS[NEQ] = {-60, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0001, 0.0, 0.0, 0.0, -60.0, 0.0, 0.0, 0.0, 0.0};

//double SC_INIVARS[NEQ] = {-60.0, 0.0, 0.0, 0.0, -60.0, 0.0, 0.0, 0.0};//1.66619e-05, 0.000197344, 0.000342989};// 0.000681673, 0.00797977, 0.0139604};
//double CC_INIVARS[NEQ] = {-60, 0.1, 0.1, 0.1
//};

/* run first in VCLAMP and then in ICLAMP */
#define ICLAMP 0
#define VCLAMP 1
#define ZAP 1
#define PULSE 2
#define SINE 3

/*Function f(t, y) should be defined in this manner */
typedef int (*RhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#define T1    RCONST(1.0)      /* first output time      */
#define RTOL  RCONST(1.0e-8)   /* scalar relative tolerance            */

#endif /* _MM_H_INCLUDED_ */
