#ifndef __SIMPLECELL3_H
#define __SIMPLECELL3_H

#include <cstring>
#include "neuron.h"

#define CM 4

#define GL p[0]
#define GM p[1]
#define GMI p[2]
#define EL  -60
#define EM -77
#define EMI -10

/* h-current activation pars */
#define VH_M p[3]
#define K_M p[4]

#define TAU1M p[5]
#define TAU2M p[6]
#define TAU_VH_M p[7]
#define TAU_K_M p[8]

/* Imi activation pars */
#define VH_MI_M p[9]
#define K_MI_M p[10]

#define TAU1_MI p[11]
#define TAU2_MI p[12]

using namespace std;


class SimpleCell3: public Neuron {
    int type;
    int pno;
    double *p;
    double iapp;
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    SimpleCell3(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~SimpleCell3();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, N_Vector, ostream&);
};

#endif
