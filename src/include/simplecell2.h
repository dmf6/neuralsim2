#ifndef __SIMPLECELL2_H
#define __SIMPLECELL2_H

#include <cstring>
#include "neuron.h"

#define CM 4

#define GL p[0]
#define GH p[1]
#define GMI p[2]
#define EL  -60
#define EH -20
#define EMI -10

/* h-current activation pars */
#define VH_H p[3]
#define K_H p[4]

#define TAU1H p[5]
#define TAU2H p[6]
#define TAU_VH_H p[7]
#define TAU_K_H p[8]

/* Imi activation pars */
#define VH_MI_M p[9]
#define K_MI_M p[10]

#define TAU1_MI p[11]
#define TAU2_MI p[12]

using namespace std;


class SimpleCell2: public Neuron {
    int type;
    int pno;
    double *p;
    double iapp;
    
  public:
        /* Make another constructor to pass in an array of parameter values *argv*/
    SimpleCell2(const string name, double *inp, int inpno, int iniVarNo, int mode);
    virtual ~SimpleCell2();
    void set_p(double *);
    int derivative(realtype , N_Vector *, N_Vector *, void *);
    void currents(realtype, N_Vector, N_Vector, ostream&);
};

#endif
