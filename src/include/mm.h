#ifndef _MM_H_INCLUDED_
#define _MM_H_INCLUDED_

#include "rk65n.h"
#include "neuron.h"
#include "synapse.h"
#include "electrode.h"
#include "abneuron.h"
 #include "pdneuron.h"
#include "simplecell.h"
#include "gap.h"
#include "neuronmodel.h"
#include "CVodeSolver.h"
#include "rk.h"

int f(realtype t, N_Vector y, N_Vector dy, void* user_data);
static int g(realtype t, N_Vector y, realtype *gout, void *user_data);
void swap(N_Vector y, N_Vector yout);
static int check_flag(void *flagvalue, char *funcname, int opt);

#define NEQ 4 /* Number of equations */
#define NEQSYN 4

//#define TSTEP  0.095367432    /* stepsize */
#define TSTEP 0.05
#define FREQUENCY 0.75
#define NUMCYCLES 30
//#define TSTOP (NUMCYCLES/FREQUENCY)*1000 // duration for constant
                                         // frequency sine waves with
                                         // integer multiple number of
                                         // cycles
#define T0 0.095367432            /* start time */
#define TSTOP 100000  /* duration */

#define NUM_PARS 15 /* number of parameters */

/* set initial calcium concentration at 10e-4 mM */
double AB_INIVARS[17] = {-38.6992, 0.0289539, 0.232727, 0.0259491, 0.125504, 0.139297, 0.613146, 0.135393, 1.8862,
                         0.0419444, 0.198714, 0.42555, -68.6553, 0.297358, 0.000221852, 0.982698, 0.000181668 };

double PD_INIVARS[17] = {-32.9142, 0.0292655, 0.223969, 0.0323433, 0.125052, 0.129203, 0.730948,
                         0.126772, 3.15326, 0.0574529, 0.184829, 0.431308, 36.636, 0.423301, 0.999846, 0.0275853, 0.000605812
};

double SC_INIVARS[4] = {-60, 0.0, 0.0, 0.0
};
double ABSTOL[4] = {1.0e-12, 1.0e-12, 1.0e-12, 1.0e-12, 
};

double SYN_INIVARS[NEQSYN] = {-60, 0.2, 0.3, 0.3};

/* run first in VCLAMP and then in ICLAMP */
#define ICLAMP 0
#define VCLAMP 1
#define ZAP 1
#define PULSE 2
#define SINE 3
#define RK65 4
#define CVODE 5

#define T1    RCONST(1.0)      /* first output time      */
#define RTOL  RCONST(1.0e-12)   /* scalar relative tolerance            */

#endif /* _MM_H_INCLUDED_ */
